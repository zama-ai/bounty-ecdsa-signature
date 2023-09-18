import matplotlib.pyplot as plt
import numpy as np

## plot this equation t(w)=T_S\cdot(2^w-1)\cdot\frac{B}{w}+T_K\cdot\lceil\frac{B}{w}\rceil
def fn(t_s, t_k, b, w):
    return t_s*(2**w-1)*b//w + t_k*np.ceil(b/w)

b = 256
t_s = 2.8
t_k = 1000

result = []
for i in range(1,12):
    result.append(fn(t_s, t_k, b, i))

plt.figure(figsize=(8,6))
plt.plot(range(1,12), result)
plt.xlabel('window size')
plt.ylabel('estimated time')
plt.title('Estimated time for different window size')
plt.savefig("window.png")