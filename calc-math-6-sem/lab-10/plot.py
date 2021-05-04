from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import pandas
import sys

# points = pandas.read_csv("res_explicit.csv")
# fig1 = plt.figure(figsize=(16, 16))
# ax1 = fig1.add_subplot(111, projection='3d')

# x = points['x'].values
# y = points['y'].values
# z = points['z'].values

# ax1.scatter(x, y, z, c='b', s=1)
# plt.savefig("res_explicit.png", dpi=500)


points = pandas.read_csv("res_shooting.csv")
fig = plt.figure(figsize=(8, 8))

ax = fig.add_subplot(111)

for i in range(1, len(points.axes[1])):
    name = 'y' + str(i - 1)

    if (i == len(points.axes[1]) - 1):
        ax.scatter(points['x'].values, points[name].values, c='r', s=5)
    else:
        ax.scatter(points['x'].values, points[name].values, c='b', s=1)

# plt.legend()
plt.title('evolution of shooting metod\'s solution')
plt.xlabel('x')
plt.ylabel('y, red is the final solution')
plt.grid()
plt.savefig("res_shooting.png", dpi=300)
