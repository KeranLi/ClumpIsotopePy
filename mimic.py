# -*- coding: utf-8 -*-
# ------------------------------------------------------------------------------------ #
# Author: Keran Li, Nanjing University, keranli98@outlook.com
# This module is mainly designed to input must information
# Use add parse to run code on the terminal
# ------------------------------------------------------------------------------------ #

from carbonate_reorder import carbonate_reorder_simulation_newton
from visualization import plot_delta47, plot_temp_and_delta47

# Like Liu et al., 2023 "碳酸盐岩团簇同位素约束下的川东地区二叠系热演化", "地质学报", if we can know the strata temprature varations
old_time = list(range(275, 0, -5))
# 1->11
simulated_temprature = [
    10, 26, 31, 46, 61, 76, 91, 106, 121, 136, 151,
    166, 181, 196, 211, 226, 241, 256, 271, 286, 276, 266,
    256, 246, 236, 226, 216, 206, 196, 186, 176, 166, 156,
    146, 136, 126, 116, 106, 91, 76, 84, 92, 100, 108, 104,
    102, 106, 112, 92, 72, 52, 47, 40, 33, 20
]
delta47_init = 0.9

# 调用函数
result = carbonate_reorder_simulation_newton(
    simulated_temprature,
    old_time,
    delta47_init
)

print(result)

fig1 = plot_delta47(old_time, result)
fig2 = plot_temp_and_delta47(old_time, simulated_temprature, result)
