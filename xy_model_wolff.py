import numpy as np
from matplotlib import pyplot as plt

def get_neighbors(i, j, n):
    '''获取点(i, j)的最近邻'''
    return [(i, (j+1)%n), (i, (j-1)%n), ((i+1)%n, j), ((i-1)%n, j)]

def spin_correlation(th_array, n):
    '''计算自旋关联函数'''
    # 周期性边界条件下, n*n的格点有2*n**2个自旋对
    # 以下循环法在数学上正确但是效率低下
    # s = 0
    # for i in range(n):
    #     for j in range(n):
    #         s += np.cos(th_array[i, j] - th_array[(i+1)%n, j]) + np.cos(th_array[i, j] - th_array[i, (j+1)%n])
    # return s/(2*n**2)
    s_c = np.sum(np.cos(th_array - np.roll(th_array, 1, axis=0)) + np.cos(th_array - np.roll(th_array, 1, axis=1)))/(2*n**2)
    return s_c

def wolff_flip(th_array, n, T):
    '''使用Wolff算法生成簇并翻转一次'''

    # 先生成簇

    # 待处理的自旋们
    to_process = [] # 这里用了python的原生方法, 可能拖累效率

    # 是否加入簇
    cluster_array = np.zeros_like(th_array, dtype=bool)

    # 生成随机的初始点
    (i_initial, j_initial) = np.random.randint(0, n, 2)

    to_process.append((i_initial, j_initial))

    while to_process: # Aieee! 原生循环! 原生循环为何......
        (i, j) = to_process.pop()
        if cluster_array[i, j] == False:
            cluster_array[i, j] = True
            for (i_, j_) in get_neighbors(i, j, n):
                # 下面的if语句利用了逻辑短路
                if cluster_array[i_, j_] == False and 1 - np.exp(-2 * np.cos(th_array[i, j] - th_array[i_, j_])/T) >= np.random.rand():
                    to_process.append((i_, j_))

    # 再翻转簇内的自旋

    # 生成随机的参考方向
    theta_0 = np.random.rand() * 2 * np.pi

    # for i in range(n):
    #     for j in range(n):
    #         if cluster_array[i, j]:
    #             th_array[i, j] = 2 * theta_0 - th_array[i, j]
    # 以上循环法效率太低, 重写之
    th_array[cluster_array] = 2 * theta_0 - th_array[cluster_array]

    return th_array

def spin_correlation_calculation(th_array, n, T):
    '''计算不同温度下的自旋关联函数'''
    s_cs = np.zeros(150)
    for i in range(1500):
        th_array = wolff_flip(th_array, n, T)
        if i % 10 == 0:
            s_cs[i//10] = spin_correlation(th_array, n)
        if i % 100 == 0:
            print('step %d, T=%.2f' % (i, T))
    return s_cs

def plot_spins(th_array, n, T):
    '''绘制XY模型系统的自旋状态'''
    fig = plt.figure(figsize=(10.8, 10.8))
    plt.quiver(np.cos(th_array), np.sin(th_array), pivot='mid', headwidth=2, headlength=6, headaxislength=6)
    plt.title('XY model at T=%.2f' % T)
    plt.xlabel('x')
    plt.ylabel('y')
    plt.axis('equal')
    plt.grid()
    fig.tight_layout()
    fig.savefig('xy_model_wolff_T%.2f.png' % T)

def main():

    # 预设kb=J=S=1, h=0, 采用周期性边界条件

    # 模型规模
    n = 50

    # 初始化位矢
    th_array_initial = np.random.rand(n, n) * 2 * np.pi
    while abs(spin_correlation(th_array_initial, n)) >= 0.05:
        th_array_initial = np.random.rand(n, n) * 2 * np.pi

    # # 绘制一个温度下自旋关联函数的随蒙特卡洛步的变化
    
    # # 温度设定
    # T_list = [0.5, 1.0, 1.5, 2.0]

    # s_c_Ts = []
    # steps = np.arange(0, 1500, 10)
    # for T in T_list:
    #     th_array = th_array_initial.copy()
    #     s_c_Ts.append(spin_correlation_calculation(th_array, n, T))

    # colour_list = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b', '#e377c2', '#7f7f7f', '#bcbd22', '#17becf']
    # fig = plt.figure(figsize=(19.2, 10.8))
    # for i in range(len(T_list)):
    #     plt.plot(steps, s_c_Ts[i], label='T=%.1f' % T_list[i], color=colour_list[i], linewidth=2)
    #     plt.axhline(np.mean(s_c_Ts[i][-20:]), linestyle='--', color = colour_list[i], linewidth=2)
    # plt.legend()
    # plt.xlabel('steps')
    # plt.ylabel('spin correlation')
    # plt.title('XY model Wolff algorithm spin correlation with Monte Carlo steps')
    # plt.grid()
    # fig.tight_layout()
    # fig.savefig('xy_model_wolff_spin_correlation.png')

    # 绘制不同温度下稳态时的自旋关联函数
    T_list = np.array([0.4, 0.6, 0.7, 0.8, 5/6, 13/15, 0.9, 14/15, 29/30, 1, 1.1, 1.2, 1.4])
    s_c_means = np.zeros_like(T_list)
    step_target = [909, 919, 929, 939, 949, 959, 969, 979, 989, 999]
    for i in range(len(T_list)):
        th_array = th_array_initial.copy()
        s_c = []
        for step in range(1000):
            th_array = wolff_flip(th_array, n, T_list[i])
            if step in step_target:
                s_c.append(spin_correlation(th_array, n))
            if step % 100 == 0:
                print('step %d, T=%.2f' % (step, T_list[i]))
        s_c_means[i] = np.mean(s_c)
    fig = plt.figure(figsize=(19.2, 10.8))
    plt.scatter(T_list, s_c_means, marker='o', s=200)
    plt.xlabel('T')
    plt.ylabel('spin correlation')
    plt.title('XY model Wolff algorithm spin correlation at different $T$')
    plt.grid()
    fig.tight_layout()
    fig.savefig('xy_model_wolff_spin_correlation_T.png')

    # # 绘制自旋状况
    # th_array = th_array_initial.copy()
    # T = 1.0
    # for i in range(1000):
    #     th_array = wolff_flip(th_array, n, T)
    #     if i % 100 == 0:
    #         print('step %d, T=%.2f' % (i, T))
    # plot_spins(th_array, n, T)

if __name__ == '__main__':
    main()