import matplotlib.pyplot as plt
import matplotlib.cbook as cbook

import numpy as np
import pandas as pd

filenames = ["cos(sin(x))", "e^sin(cos(x))",
             "ln(x + 3)", "sin^2(x)", "sqrt(x + 3)"]
colors = ['b', 'r', 'g', 'c', 'm']

for filename in filenames:
    msft = pd.read_csv("output/calculator/" + filename +
                       ".csv")

    ax = plt.gca()

    for i in range(1, 6):
        msft.plot(x=0, y=i, logy=True, logx=True,
                  subplots=True, sharex=True, ax=ax, color=colors[i - 1], linewidth=.5)

    plt.ylabel("error")
    plt.xlabel("h")
    plt.grid(True)

    plt.savefig("output/drawer/" +
                filename + ".pdf")
    plt.close()
