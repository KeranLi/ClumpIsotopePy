import matplotlib.pyplot as plt

def plot_temp_and_delta47(old_time, simulated_temprature, delta47_t_list):
    # 绘图
    fig, ax1 = plt.subplots()

    color = 'tab:red'
    ax1.set_xlabel('From Perianm (Ma)')
    ax1.set_ylabel('simulated temprature', color=color)
    ax1.plot(old_time, simulated_temprature, color=color)
    ax1.tick_params(axis='y', labelcolor=color)

    ax2 = ax1.twinx()  

    color = 'tab:blue'
    ax2.set_ylabel('Δ47', color=color)  
    ax2.plot(old_time, delta47_t_list, color=color)
    ax2.tick_params(axis='y', labelcolor=color)

    plt.gca().invert_xaxis()  # 反转x轴
    plt.show()

def plot_delta47(old_time, delta47_t_list):
    # 绘图
    fig, ax2 = plt.subplots()

    color = 'tab:blue'
    ax2.set_ylabel('Δ47', color=color)
    ax2.set_xlabel('From Perianm (Ma)', color=color)  
    ax2.plot(old_time, delta47_t_list, color=color)
    ax2.tick_params(axis='y', labelcolor=color)

    plt.gca().invert_xaxis()  # 反转x轴
    plt.show()