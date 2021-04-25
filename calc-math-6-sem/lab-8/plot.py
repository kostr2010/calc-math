from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import pandas
import sys

points = pandas.read_csv(sys.argv[1])

fig = plt.figure(figsize=(16, 16))
ax = fig.add_subplot(111, projection='3d')

x = points['x'].values
y = points['y'].values
z = points['z'].values

ax.scatter(x, y, z, c='b', s=1)

plt.savefig(sys.argv[1][:-4] + '.png', dpi=500)
