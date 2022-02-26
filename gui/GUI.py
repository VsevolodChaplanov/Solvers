import numpy as np
from matplotlib import pyplot as plt
import pandas as pd
from mpl_toolkits import mplot3d
import seaborn as sns

# ТУТ МЕНЯЕМ В СООТВЕСТВИИ С С++ файлом
Nx = 101

Ny = 101

df = pd.read_csv("2D.csv", dtype = np.double, names = ["heat"], index_col = False)
Field = np.array(df["heat"].values).reshape(Nx, Ny)
x = np.arange(0, 1.0 + 1 / (Nx), 1/(Nx-1))
y = np.arange(0, 1.0 + 1 / (Ny), 1/(Ny-1))
X, Y = np.meshgrid(x, y)
fig, ax = plt.subplots()
CS = ax.contour(X, Y, Field, levels = 20)
ax.clabel(CS, inline=True, fontsize=10)
fig.set_figwidth(9)    #  ширина и
fig.set_figheight(9)    #  высота "Figure"
plt.savefig("Solution " + str((len(df["heat"])-1)) + ".png")