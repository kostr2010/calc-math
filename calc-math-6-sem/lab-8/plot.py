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

x1 = points['y1_11'].values
y1 = points['y2_11'].values

ax2.scatter(x1, y1, c='b', s=1)

x2 = points['y1_12'].values
y2 = points['y2_12'].values

ax2.scatter(x2, y2, c='r', s=1)

# plt.plot(x2, y2)

# x3 = points['y1_21'].values
# y3 = points['y2_21'].values

# plt.plot(x3, y3)

# x4 = points['y1_22'].values
# y4 = points['y2_22'].values

# plt.plot(x4, y4)

plt.savefig("res_implicit.png", dpi=500)
