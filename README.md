## The carbonate clumped isotope calulcation by Python

#### Author: Keran Li (Nanjing University)

#### This repository is licensed under the ***MIT License***.

### Introduction

This repository contains the code files for the carbonate clumped isotope calculation by Python. 

### Dependencies

* Python 3.10
* Jupyter Notebook
* Numpy
* Scipy
* Matplotlib
* Pandas
  
### What has been done?

1. ~~Realize paper "Hemingway, J. D., and Henkes, G. A disordered kinetic model for clumped isotope bond reordering in carbonates, 2021, EPSL."~~
2. ~~郭炀锐, 邓文峰, 韦刚健. 碳酸盐成岩作用中的团簇同位素地球化学研究进展. 2022, 矿物岩石地球化学通报.~~
3. ~~The carbonate clumped isotope reordering calculation by Python (Exchange/diffusion model from Stolper et al., 2015| Paper "Stolper, D. A., Eiler, J. M., THE KINETICS OF SOLID-STATE ISOTOPE-EXCHANGE REACTIONS FOR CLUMPED ISOTOPES: A STUDY OF INORGANIC CALCITES AND APATITES FROM NATURAL AND EXPERIMENTAL SAMPLES. 2015．American Journal of Science. ").~~

### Results

1. Single initial Δ47 input
    (1) ours (time bar is converted):
    <div align="center">
    <img width="500" alt="image" src="https://user-images.githubusercontent.com/66153455/280744777-30c49076-eb05-42c7-b12f-0e7e169b6bd9.png">
    </div>
    (2) 刘鑫, 邱楠生, 冯乾乾. 碳酸盐岩团簇同位素约束下的川东地区二叠系热演化. 2023, 地质学报. results:
    <div align="center">
    <img width="500"  alt="image" src="https://user-images.githubusercontent.com/66153455/280745312-d4974462-7839-429e-9a60-23fb31db4722.png">
    </div>

### To do

1. Add Monte-Carlo simulation in defferent Δ47 input
