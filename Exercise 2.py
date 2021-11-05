import numpy as np
import random
import math
import matplotlib.pyplot as plt

def single_sweep(initial_config,N,J,T,h):
    initial_config_energy = energy(initial_config, N, J,h)
    magnetization = np.array([])
    Energy=np.array([])
    for sweep in range(N**2):
        row,column=random.randint(0,N-1),random.randint(0,N-1)
        new_config=np.copy(initial_config)
        new_config[row][column] = -(new_config[row][column])
        new_config_energy=energy(new_config,N,J,h)
        delta_s=(new_config_energy/T) - (initial_config_energy/T)
        if  delta_s < 0:
            magnetization=np.append(magnetization,np.sum(new_config))
            Energy=np.append(Energy,math.exp(-new_config_energy/T))
            initial_config=np.copy(new_config)
        else:
            y=random.uniform(0,1)
            if y<= math.exp(-delta_s):
                magnetization = np.append(magnetization, np.sum(new_config))
                Energy = np.append(Energy, math.exp(-new_config_energy/T))
                initial_config = np.copy(new_config)
    return magnetization,Energy


def energy(config,N,J,h):
    sum=0
    for i in range(N):
        for j in range(N):
            sum += config[i][j] * (config[(i + 1) % N][j] + config[(i - 1)][j] + config[i][(j + 1) % N] + config[i][(j - 1)])
    return (-J/2)*sum - (h*np.sum(config))

def average_magnetization(Energy,magnetization,N):
    avg_m=0
    for i in range(Energy.size):
        avg_m += magnetization[i]*Energy[i]
    return (avg_m / np.sum(Energy))/(N**2)

def main():
    grid_size = 5
    J=0.25
    T=1
    iterations=1000
    h_range = np.arange(-1, 1.1, 0.1)
    fixxed_N_m_list=[]
    initial_config = np.empty([grid_size, grid_size])
    for i in range(grid_size):
        for j in range(grid_size):
            choice = random.uniform(0, 1)
            if choice <= 0.5:
                initial_config[i][j] = 1
            elif 0.5 < choice <= 1:
                initial_config[i][j] = -1

    for h in h_range:
        m_list = []
        for k in range(iterations):
            magnetization, Energy = single_sweep(initial_config,grid_size,J,T,h)
            avg_m=average_magnetization(Energy,magnetization,grid_size)
            m_list.append(avg_m)
        fixxed_N_m_list.append(sum(m_list)/iterations)
    plt.plot(h_range,fixxed_N_m_list,'b')
    plt.xlabel("h")
    plt.ylabel("<m>")
    plt.title("Plot for fixed grid size (J=0.25,grid size= 10*10,T=1, iterations = 1000 )")
    plt.show()

main()

