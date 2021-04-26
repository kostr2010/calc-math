from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import pandas
import sys

points = pandas.read_csv("res_explicit.csv")
fig1 = plt.figure(figsize=(16, 16))
ax1 = fig1.add_subplot(111, projection='3d')

x = points['x'].values
y = points['y'].values
z = points['z'].values

ax1.scatter(x, y, z, c='b', s=1)
plt.savefig("res_explicit.png", dpi=500)


points = pandas.read_csv("res_implicit.csv")
fig2 = plt.figure(figsize=(16, 16))
ax2 = fig2.add_subplot(111)

x = points['y1'].values
y = points['y2'].values

ax2.scatter(x, y, c='b', s=1)
plt.savefig("res_implicit.png", dpi=500)
