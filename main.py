# -*- coding: utf-8 -*-
"""
Created on Sat Sep 28 2024

@author: chen

Copyright 2024 Edgar-Mendonca/SHPB-Analysis. Github URL: https://github.com/Edgar-Mendonca/SHPB-Analysis
"""


import tkinter as tk
from tkinter import filedialog

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import SpanSelector



from scipy.signal import savgol_filter
from scipy.integrate import cumulative_trapezoid


time_data = []
select_incident = []
select_reflected = []
select_transmitted = []

#主体
def main():
    file_path = filedialog.askopenfilename(filetypes=[("CSV Files", "*.csv")])
    
    if file_path:
        data = pd.read_csv(file_path)
        crop(data)


#裁剪
def crop(data):
    global time_data
    time_data = data['Time'].values
    incident_voltage = data['Incident'].values
    transmitted_voltage = data['Transmitted'].values

    fig, ((ax1, ax2), (ax3, ax4), (ax5, ax6)) = plt.subplots(3, 2, figsize=(12, 8))

    ax1.plot(time_data, incident_voltage)
    ax1.set_title('Press left mouse button and drag to select a region of Incident Voltage')
    line1, = ax2.plot([], [])
    
    ax3.plot(time_data, incident_voltage)
    ax3.set_title('Press left mouse button and drag to select a region of Reflected Voltage')
    line2, = ax4.plot([], [])

    ax5.plot(time_data, transmitted_voltage)
    ax5.set_title('Press left mouse button and drag to select a region of Tranmitted Voltage')
    line3, = ax6.plot([], [])
    
    #入射
    def onselect_incident(xmin, xmax):
        indmin, indmax = np.searchsorted(time_data, (xmin, xmax))
        indmax = min(len(time_data) - 1, indmax)

        region_x = time_data[indmin:indmax]
        region_y = incident_voltage[indmin:indmax]

        if len(region_x) >= 2:
            line1.set_data(region_x, region_y)
            ax2.set_xlim(region_x[0], region_x[-1])
            ax2.set_ylim(region_y.min(), region_y.max())
            fig.canvas.draw_idle()

            #保存数据
            global select_incident
            select_incident = region_y
    
    global span_incident
    span_incident = SpanSelector(
        ax1,
        onselect_incident,
        "horizontal",
        useblit=True,
        props=dict(alpha=0.5, facecolor="tab:red"),
        interactive=True,
        drag_from_anywhere=True
    )

    #反射
    def onselect_reflected(xmin, xmax):
        indmin, indmax = np.searchsorted(time_data, (xmin, xmax))
        indmax = min(len(time_data) - 1, indmax)

        region_x = time_data[indmin:indmax]
        region_y = incident_voltage[indmin:indmax]  # Change to reflected voltage data

        if len(region_x) >= 2:
            line2.set_data(region_x, region_y)
            ax4.set_xlim(region_x[0], region_x[-1])
            ax4.set_ylim(region_y.min(), region_y.max())
            fig.canvas.draw_idle()

            #保存数据
            global select_reflected
            select_reflected = region_y

    global span_reflected
    span_reflected = SpanSelector(
        ax3,
        onselect_reflected,
        "horizontal",
        useblit=True,
        props=dict(alpha=0.5, facecolor="tab:green"),
        interactive=True,
        drag_from_anywhere=True
    )


    #透射
    def onselect_transmitted(xmin, xmax):
        indmin, indmax = np.searchsorted(time_data, (xmin, xmax))
        indmax = min(len(time_data) - 1, indmax)

        region_x = time_data[indmin:indmax]
        region_y = transmitted_voltage[indmin:indmax]  

        if len(region_x) >= 2:
            line3.set_data(region_x, region_y)
            ax6.set_xlim(region_x[0], region_x[-1])
            ax6.set_ylim(region_y.min(), region_y.max())
            fig.canvas.draw_idle()

            #保存数据
            global select_transmitted
            select_transmitted = region_y

            analysis()

    global span_transmitted
    span_transmitted = SpanSelector(
        ax5,
        onselect_transmitted,
        "horizontal",
        useblit=True,
        props=dict(alpha=0.5, facecolor="tab:blue"),
        interactive=True,
        drag_from_anywhere=True
    )
    
    plt.tight_layout()
    plt.show()


#分析
def analysis():
    #统一值的数量
    min_length = min(len(select_incident), len(select_reflected), len(select_transmitted))

    selected_time = time_data[:min_length]
    selected_incident = select_incident[:min_length]
    selected_reflected = select_reflected[:min_length]
    selected_transmitted = select_transmitted[:min_length]
    
    global entries
    if entries[0].get() == "":
        #与实验装置和试件有关的参数 (自动代入)
        Eb = 2.06e5  #杆的弹性模量 (MPa)
        rho_b = 7850  #杆的密度 (kg/m^3)
        Ab = 7850  #杆的横截面积 (mm^2)
        lo = 70.7  #试件的初始长度 (mm)
        As = 4998.49  #试件的横截面积(mm^2)
        Voltage_to_Strain = 1e-3  #电压转化为应变的系数
        
        #滤波器参数
        window_length = 51  #取值为奇数且不能超过len(x)。越大，则平滑效果越明显；越小，则更贴近原始曲线。
        poly_order = 3  #多项式拟合的阶数。越小，则平滑效果越明显；越大，则更贴近原始曲线。
        
    else:
        #与实验装置和试件有关的参数 (手动输入)
        Eb = float(entries[0].get())  #2.06e5  #杆的弹性模量 (MPa)
        rho_b = float(entries[1].get())  #7850  #杆的密度 (kg/m^3)
        Ab = float(entries[2].get())  #7850  #杆的横截面积 (mm^2)
        lo = float(entries[3].get())  #70.7  #试件的初始长度 (mm)
        As = float(entries[4].get())  #4998.49  #试件的横截面积(mm^2)
        Voltage_to_Strain = float(entries[5].get())  #1e-3  #电压转化为应变的系数
        
        #滤波器参数
        window_length = int(entries[6].get())  #51  #取值为奇数且不能超过len(x)。越大，则平滑效果越明显；越小，则更贴近原始曲线。
        poly_order = int(entries[7].get())  #3  #多项式拟合的阶数。越小，则平滑效果越明显；越大，则更贴近原始曲线。
    
    #滤波
    selected_incident_filtered = -1 * savgol_filter(selected_incident, window_length, poly_order)
    selected_reflected_filtered = savgol_filter(selected_reflected, window_length, poly_order)
    selected_transmitted_filtered = -1 * savgol_filter(selected_transmitted, window_length, poly_order)
    
    #将电压转化为应变
    strain_incident = Voltage_to_Strain * selected_incident_filtered
    strain_reflected = Voltage_to_Strain * selected_reflected_filtered
    strain_transmitted = Voltage_to_Strain * selected_transmitted_filtered
    
    #杆波速 (mm/s)
    cb = np.sqrt(Eb / rho_b) * 1e6
    
    #计算试件应变率
    es_dot = cb /lo * (strain_incident + strain_reflected - strain_transmitted)
    
    #利用scipy.integrate中的culative_梯形积分应变率得到应变
    es = cumulative_trapezoid(es_dot, selected_time, initial=0)

    #计算试样的应力
    Ss = (Eb * Ab / As /2) * (strain_incident - strain_reflected + strain_transmitted)

    #用εs(t) = -ln(1 - es(t))计算真应变:
    epsilon_tolerance = 1e-10  #小常数防止被零除
    es_safe = np.clip(es, epsilon_tolerance, 1 - epsilon_tolerance)  #裁剪es以确保其保持在有效范围内
    true_strain = -np.log(1 - es_safe)

    #用σs(t) = Ss(t) * (1 - es(t))计算真应力
    true_stress = Ss * (1 - es)

    #结果曲线图
    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(12, 8))

    #滤波后的电压时间
    ax1.plot(selected_time, selected_incident_filtered, label="Incident Voltage (Filtered)", color="red")
    ax1.plot(selected_time, selected_reflected_filtered, label="Reflected Voltage (Filtered)", color="green")
    ax1.plot(selected_time, selected_transmitted_filtered, label="Transmitted Voltage (Filtered)", color="blue")
    ax1.legend()
    ax1.set_title("Filtered Voltage vs Time")

    #应力-应变率曲线
    ax2.plot(selected_time, es_dot, color='red')
    ax2.set_title("Strain Rate vs Time")

    #应力-应变曲线
    ax3.plot(es, Ss, color='green')
    ax3.set_title("Stress vs Strain")

    #真应力-应变曲线
    ax4.plot(true_strain, true_stress, color='blue')
    ax4.set_title("True Stress vs True Strain")

    plt.tight_layout()
    plt.show()


root = tk.Tk()
root.title("Open CSV File")
root.geometry("800x600")

labels = ["    Young's modulus of the bars (MPa)", 
          "    Density of the bars (kg/m^3)", 
          "    Cross-sectional area of the bars (mm^2)", 
          "    Initial length of the specimen (mm)", 
          "    Cross-sectional area of the specimen (mm^2)", 
          "    Conversion of Voltage to Strain", 
          "    Adjust window size", 
          "    Adjust polynomial order"]
entries = []

for label in labels:
    frame = tk.Frame(root)
    frame.pack(fill='x', pady=10)
    
    lbl = tk.Label(frame, text=label, width=40, anchor='w')
    lbl.pack(side=tk.LEFT)
    
    entry = tk.Entry(frame, width=20)
    entry.pack(side=tk.LEFT)
    entries.append(entry)

button = tk.Button(root, text="Open File", command=main)
button.pack(expand=True)

root.mainloop()