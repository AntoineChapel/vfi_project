import numpy as np
import matplotlib.pyplot as plt


kmin = 1
kmax = 25
precision=200

kgrid = np.linspace(kmin, kmax, precision)
gk = np.linspace(kmin, kmax, precision)
Vk0 = np.ones(precision)

A=10
α=0.5
β=0.9

norm = 1e5
tol = 1e-6
maxiter=1000
n_iter=0


Vk = Vk0
while n_iter < maxiter and norm > tol:
    value_array = np.empty((precision, precision))
    
    for iprim, kprim in enumerate(kgrid):
        for i, k in enumerate(kgrid):
            c = A*(k**α) - kprim
            if c > 0:
                value_array[i, iprim] = np.log(c) + β*Vk[iprim] 
            else:
                value_array[i, iprim] = -np.inf
                
    Vkprim = np.empty(precision)
    for row in range(value_array.shape[1]):
        gk[row] = kgrid[int(np.argmax(value_array[row, :]))]
        Vkprim[row] = np.max(value_array[row, :])
    
    norm = np.max(np.abs(Vkprim - Vk))
    Vk = Vkprim
    
    n_iter += 1
    print("iteration: ", n_iter, " norm: ", norm)


kstar = kgrid[np.argmin(np.abs(gk - kgrid))]
print(f"The steady-state value of capital is {kstar}")


plt.rcParams["figure.figsize"] = (10, 8)

plt.plot(kgrid, gk, label="Optimal k'")
plt.plot(kgrid, kgrid, label="45° line")
plt.vlines(kstar, ymin=0, ymax=30, label="Equilibrium value k*")

plt.title("Value Function Iteration: ")
plt.legend()
plt.savefig("vfipy.png")
