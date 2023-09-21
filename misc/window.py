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
# add annotation on the plot
for i in range(4,7):
    plt.annotate('(%d, %d)' % (i, result[i-1]), xy=(i, result[i-1]), xytext=(i, result[i-1]+5000))
i = 7
plt.annotate('(%d, %d)' % (i, result[i-1]), xy=(i, result[i-1]), xytext=(i, result[i-1]+15000))
i = 8
plt.annotate('(%d, %d)' % (i, result[i-1]), xy=(i, result[i-1]), xytext=(i, result[i-1]+20000))


plt.xlabel('window size (bits)')
plt.ylabel('estimated time (seconds)')
plt.title('Estimated time vs window size for 256-bit on n2d-standard-64')
plt.savefig("window.png")