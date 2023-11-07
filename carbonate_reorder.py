# -*- coding: utf-8 -*-
# ------------------------------------------------------------------------------------ #
# Author: Keran Li, Nanjing University, keranli98@outlook.com
# This module is mainly designed to input must information
# Use add parse to run code on the terminal
# ------------------------------------------------------------------------------------ #

# 在 carbonate_thermal_simulation.py 文件中定义以下函数
import numpy as np
from scipy.integrate import odeint

def carbonate_reorder_simulation(simulated_temprature, old_time, Δ47eq):
    # 参数
    T = simulated_temprature  # 温度
    t = old_time  # 时间区间
    Δ47 = [0.58]  # 初始化Δ47数组

    # 定义微分方程
    def reorder(Δ47_ode, t, T, Δ47eq):
        Δ47 = Δ47_ode[0]

        for i in range(len(t)):
            # 使用T[i]取温度值
            exp_term = -23000 / (T[i] * (1e+6) * 365 * 24 * 60 * 60 + 1) + 20
            exp_value = np.exp(exp_term)
            Δ47_new = Δ47eq + (Δ47 - Δ47eq) * np.exp(-exp_value * t[i])

            # 添加判断
            if Δ47 > 0.6:
                Δ47_new *= 0.25

            Δ47_ode[0] = Δ47_new

        return Δ47_ode

    # 求解Δ47t
    Δ47t = odeint(reorder, Δ47, t, args=(T, Δ47eq), rtol=1e-5, atol=1e-8)[:, -1] / 100
    # 第一个值特殊处理
    Δ47t[0] = Δ47t[0] * (100 if t[0] == 0 else 1)

    return Δ47t