import matplotlib.pyplot as plt
import numpy as np
import os

import oned_consolidation as oc

fig = plt.figure()
plot1 = fig.subplots(1, 1)
plot1.set_title('Settlement - Time relation')
plot1.set_xlabel('Time')
plot1.set_ylabel('Settlement')
#plot1.set_xlim([0.0, 1.0])

x_data = [];
y_data = [];

with open("res.txt", 'r') as res_f:
    while True:
        res_line = res_f.readline()
        if not res_line:
            break
        data = list(map(lambda x: float(x), res_line.split(',')[:-1]))
        #print(data)
        if data[0] >=0.0 and data[0] <=3.0:
            x_data.append(data[0])
            y_data.append(data[22])
    
    plot1.plot(x_data, y_data)

plt.show()
