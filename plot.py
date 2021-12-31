import numpy as np
import matplotlib.pyplot as plt

r = np.loadtxt('fflt.txt')
rr = np.loadtxt('pyflt.txt')
plt.plot(r)
plt.plot(rr, '--')
# plt.xlim(1000, 4000)
plt.savefig('cmp.png')