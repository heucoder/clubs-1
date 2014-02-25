import matplotlib.pyplot as plt
from pylab import *

import matplotlib.cm as cm

import sys
import numpy as np

close()
close()

# read the data
f = open('data.txt', 'r')
lines = f.readlines()
file_len = len(lines)
no_clusters = file_len/2

print no_clusters
fig2 = plt.figure()

colors = cm.rainbow(np.linspace(0, 1, no_clusters))
# colors = ["r","g","b","m","c"]

for i in range(0,no_clusters):
    x=[]
    y=[]
    line_x = lines[2*i]
    line_y = lines[2*i+1]
    line_arr_x = line_x.split(' ')
    line_arr_y = line_y.split(' ')
    line_arr_len = len(line_arr_x)-1
    for j in range(0,line_arr_len):
        x.append(float(line_arr_x[j]))
        y.append(float(line_arr_y[j]))
    plt.scatter(x,y,color=colors[i])

plt.show()

print x
print y

f.close()

# read the true values
f = open('TrainingData.csv')
lines = f.readlines()
file_len = len(lines)
no_clusters = 5

x1=[]
y1=[]
x2=[]
y2=[]
x3=[]
y3=[]
x4=[]
y4=[]
x5=[]
y5=[]
count=0
for i in range(0,file_len):
    line = lines[i]
    line_arr = line.split(',')
    if (int(line_arr[0])==6733):
        count = count+1
        if (str(line_arr[5][0:2]) == str("FM")):
            x1.append(float(line_arr[2])*0.001)
            y1.append(float(line_arr[3])*0.001)
        if (str(line_arr[5][0:2]) == str("TT")):
            x2.append(float(line_arr[2])*0.001)
            y2.append(float(line_arr[3])*0.001)
        if (str(line_arr[5][0:4]) == str("Both")):
            x3.append(float(line_arr[2])*0.001)
            y3.append(float(line_arr[3])*0.001)
        if (str(line_arr[5][0:12]) == str("Undetermined")):
            x4.append(float(line_arr[2])*0.001)
            y4.append(float(line_arr[3])*0.001)
        if (str(line_arr[5][0:8]) == str("Aberrant")):
            x5.append(float(line_arr[2])*0.001)
            y5.append(float(line_arr[3])*0.001)

print count
fig1 = plt.figure()
plt.scatter(x1,y1,color="g")
plt.scatter(x2,y2,color="r")
plt.scatter(x3,y3,color="b")
plt.scatter(x4,y4,color="m")
plt.scatter(x5,y5,color="c")

plt.show()

