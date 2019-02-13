import matplotlib.pyplot as plt
import numpy as np
import os

import oned_consolidation as oc

fig = plt.figure()
plot1 = fig.subplots(1, 1)
plot1.set_title('Pore Pressure - Time Relation')
plot1.set_xlabel('Time')
plot1.set_ylabel('Pore Pressure')
plot1.set_xlim([0.0, 1.0])

x_data = [];
y_data = [];

with open("res.txt", 'r') as res_f:
    while True:
        res_line = res_f.readline()
        if not res_line:
            break
        data = list(map(lambda x: float(x), res_line.split(',')[:-1]))
        #print(data)
        
        if data[0] >=1.0 and data[0] <=2.0:
            x_data.append(data[0]-1.0)
            y_data.append(data[22])
    
    line1, = plot1.plot(x_data, y_data)

Es = 40.0e6
kv = 1.0e-5
miu = 1.0 # dynamic viscosity
Cv = kv * Es / miu
u0 = 40.0e3
H = 10.0
con_res = oc.OneDConsolidation(Cv, Es, u0, H)

data_num = 100
t_list = np.zeros(data_num)
p_list = np.zeros(data_num)
u_list = np.zeros(data_num)
for i in range(data_num):
    t_list[i] = 0.01 * float(i)
    p_list[i] = con_res.calPorePressure(t_list[i], 10.0)
    u_list[i] = con_res.calSettlement(t_list[i])


line2, = plot1.plot(t_list, p_list, 'k--')
#plot1.plot(t_list, u_list, 'r--')

plt.legend(handles=[line1, line2], labels=['FEM du-form', 'Analytical Solution'])
plt.show()
