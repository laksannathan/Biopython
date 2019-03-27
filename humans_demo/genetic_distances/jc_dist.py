import scipy.io
import matplotlib.pyplot as plt

d = scipy.io.loadmat('D.mat')

for key, value in d.items() :
    mat = value

print(mat)

plt.imshow(mat, cmap='hot', interpolation='nearest')
plt.show()
