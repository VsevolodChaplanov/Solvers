from matplotlib import pyplot as plt
import pandas as pd
import numpy as np

df = pd.read_csv("Residuals.csv", names = ["Residual"])
fig, ax = plt.subplots()
df2 = pd.read_csv("2D.csv", names = ["U(x) Numerical solution"]) 
df["N"] = np.arange(0, len(df["Residual"]) * 10, step = 10)

ax.plot("N" , "Residual", data = df, marker = "*", linestyle = "-")
ax.set_title(r'$||r||_2$ versus iteration number',
						fontsize = 20)
ax.set_xlabel("Iteration")
ax.set_ylabel(r"${\rm ln}(||r||_2)$")
ax.set_yscale("log")
ax.legend()
ax.grid()
ax.minorticks_on()
fig.set_figwidth(16)    #  ширина и
fig.set_figheight(9)    #  высота "Figure"

plt.savefig("Problems/GD/Residuals " + str(len(df2["U(x) Numerical solution"])) + ".png")
plt.cla()
