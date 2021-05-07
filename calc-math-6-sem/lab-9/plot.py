from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import pandas
import sys

points = pandas.read_csv("res.csv")
fig1 = plt.figure(figsize=(8, 8))
ax1 = fig1.add_subplot(111)

x = points['x'].values
y = points['y'].values

ax1.scatter(x, y, c='b', s=1)
plt.savefig("res.png", dpi=100)
