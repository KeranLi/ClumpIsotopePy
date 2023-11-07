# -*- coding: utf-8 -*-
# ------------------------------------------------------------------------------------ #
# Author: Keran Li, Nanjing University, keranli98@outlook.com
# This module is mainly designed to input must information
# Use add parse to run code on the terminal
# ------------------------------------------------------------------------------------ #

from carbonate_reorder import carbonate_reorder_simulation

# Like Liu et al., 2023 "碳酸盐岩团簇同位素约束下的川东地区二叠系热演化", "地质学报", if we can know the strata temprature varations
old_time = list(range(0, 275, 25))
simulated_temprature = [25, 75, 130, 140, 150, 200, 225, 140, 100, 50, 25]
Δ47eq = 0.6  # 平衡Δ47值

# 调用函数
result = carbonate_reorder_simulation(simulated_temprature, old_time, Δ47eq)
print(result)