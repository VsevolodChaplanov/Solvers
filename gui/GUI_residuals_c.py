from matplotlib import pyplot as plt
import pandas as pd
import numpy as np

df_conv = 		pd.DataFrame()
CGP_convergence = 	pd.read_csv("Residuals_P.csv", dtype = np.double, names = ["CGP_Residue_vs_Iteration"], index_col = False)
CG_convergence = 	pd.read_csv("Residuals.csv", dtype = np.double, names = ["CG_Residue_vs_Iteration"], index_col = False)
df2 = pd.read_csv("2D.csv", dtype = np.double, names = ["heat"], index_col = False)
#
# Сходимость
#

N = max(len(CGP_convergence["CGP_Residue_vs_Iteration"]), len(CG_convergence["CG_Residue_vs_Iteration"]))

df_conv = pd.concat([CG_convergence, CGP_convergence] , axis = 1)
df_conv["N"] = np.arange(0, N * 10, step = 10)
fig, ax = plt.subplots()
ax.plot("N", "CGP_Residue_vs_Iteration", data = df_conv, marker = "*", linestyle = "-", label = "Preconditioned CG")
ax.plot("N", "CG_Residue_vs_Iteration", data = df_conv, marker = "o", linestyle = "-", label = "CG")
ax.set_title(r'$-\dfrac{d}{dx}(k(x)\cdot\dfrac{d}{dx}(u(x)) = f(x)$',
														fontsize = 20)
ax.set_xlabel("Iterations")
ax.set_ylabel(r"${\rm ln}(||r||_2)$")
ax.set_yscale("log")
ax.legend()
ax.grid()
ax.minorticks_on()
fig.set_figwidth(16)    #  ширина и
fig.set_figheight(9)    #  высота "Figure"

plt.savefig("Residuals Scaling " + str(len(df2["heat"])) + ".png")
plt.cla()