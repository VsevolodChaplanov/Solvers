from matplotlib import pyplot as plt
import pandas as pd
import numpy as np


def exact_solution(x):

	return np.array(x - np.sin(16 * np.pi * x) / 16)

df = pd.read_csv("Solution.csv", names = ["U(x) Numerical solution"]) # Считать решение из файла
df["X"] = np.arange(0,1 + 1/(len(df["U(x) Numerical solution"])-1), 1/(len(df["U(x) Numerical solution"])-1)) # Построить равномерную одномерную сетку
df["Exact solution"] = exact_solution(df["X"])
# Построение графиков
fig, ax = plt.subplots()
ax.plot("X" , "Exact solution", data = df, marker = "o", linestyle = None)
ax.plot("X" , "U(x) Numerical solution", data = df, marker = "*", linestyle = "-")
ax.set_title(r'$-\dfrac{d}{dx}(k(x)\cdot\dfrac{d}{dx}(u(x)) = f(x)$',
						fontsize = 20)
ax.set_xlabel("x")
ax.set_ylabel("u(x)")
ax.legend()
ax.grid()
ax.minorticks_on()
fig.set_figwidth(16)    #  ширина и
fig.set_figheight(9)    #  высота "Figure"


plt.savefig("Solution " + str((len(df["U(x) Numerical solution"])-1)) + ".png")
plt.cla()