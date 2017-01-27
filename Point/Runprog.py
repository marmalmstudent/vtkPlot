import numpy as np
import os


xSize = 1001
ySize = 1001
filename = "/home/marcus/git/vtkPlot/Point/filename.txt"
A, B = np.meshgrid(np.linspace(-2*np.pi, 2*np.pi, xSize), np.linspace(-2*np.pi, 2*np.pi, ySize))
C = np.sin(A**2)*np.sin(B**2)
X = np.reshape(A, (-1, np.size(A)))
Y = np.reshape(B, (-1, np.size(B)))
Z = np.reshape(C, (-1, np.size(C)))
np.savetxt("X.txt", X)
np.savetxt("Y.txt", Y)
np.savetxt("Z.txt", Z)
os.system("build/Point "+filename+" "+str(xSize)+" "+str(ySize))
