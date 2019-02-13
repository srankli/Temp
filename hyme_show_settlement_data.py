import matplotlib.pyplot as plt
import numpy as np
import os

import oned_consolidation as oc

fig = plt.figure()
plot1 = fig.subplots(1, 1)
plot1.set_title('Settlement - Time Relation')
plot1.set_xlabel('Time')
plot1.set_ylabel('Settlement')
plot1.set_xlim([0.0, 0.5])

x_data_FEM_a = [];
y_data_FEM_a = [];
x_data_FEM_du = [];
y_data_FEM_du = [];
x_data_MPM_a = [];
y_data_MPM_a = [];

file_FEM_a = open("res_FEM_a.txt", 'r')
while True:
    res_line = file_FEM_a.readline()
    if not res_line:
        break
    data = list(map(lambda x: float(x), res_line.split(',')[:-1]))
    if data[0] >=1.0 and data[0] <=2.0:
        x_data_FEM_a.append(data[0]-1.0)
        y_data_FEM_a.append(data[21])
file_FEM_a.close()

file_FEM_du = open("res_FEM_du.txt", 'r')
while True:
    res_line = file_FEM_du.readline()
    if not res_line:
        break
    data = list(map(lambda x: float(x), res_line.split(',')[:-1]))
    if data[0] >=1.0 and data[0] <=2.0:
        x_data_FEM_du.append(data[0]-1.0)
        y_data_FEM_du.append(data[21])
file_FEM_du.close()

file_MPM_a =  open("res_MPM_a.txt", 'r')
while True:
    res_line = file_MPM_a.readline()
    if not res_line:
        break
    data = list(map(lambda x: float(x), res_line.split(',')[:-1]))
    if data[0] >=1.0 and data[0] <=2.0:
        x_data_MPM_a.append(data[0]-1.0)
        y_data_MPM_a.append(data[50]-9.9)

file_MPM_a.close()

line1, = plot1.plot(x_data_FEM_a, y_data_FEM_a, "b--")
line2, = plot1.plot(x_data_FEM_du, y_data_FEM_du, "g--")
line3, = plot1.plot(x_data_MPM_a, y_data_MPM_a, "r--")

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


#plot1.plot(t_list, p_list, 'k--')
line4, = plot1.plot(t_list, u_list, 'k')

plt.legend(handles=[line1,line2, line3, line4], labels=['FEM a-form','FEM du-form', 'MPM a-form', 'Analytical Solution'])

plt.show()
