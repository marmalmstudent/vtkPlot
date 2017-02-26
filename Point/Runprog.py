import numpy as np
import os


xSize = 501
ySize = 3001
filename = "/home/marcus/git/vtkPlot/Point/filename.txt"
A, B = np.meshgrid(np.linspace(-2*np.pi, 2*np.pi, xSize),
                   np.linspace(-2*np.pi, 2*np.pi, ySize))
R = np.sin(A)*np.cos(B)
C = 5*np.sin(R)**2
X = np.reshape(A, (-1, np.size(A)))
Y = np.reshape(B, (-1, np.size(B)))
Z = np.reshape(C, (-1, np.size(C)))
fh = open("X.dat", "wb")
fh.write(X.tobytes())
fh.close()
fh = open("Y.dat", "wb")
fh.write(Y.tobytes())
fh.close()
fh = open("Z.dat", "wb")
fh.write(Z.tobytes())
fh.close()
os.system("build/Point $HOME/git/vtkPlot/Point/X.dat "
          +"$HOME/git/vtkPlot/Point/Y.dat "
          +"$HOME/git/vtkPlot/Point/Z.dat "
          +str(xSize)+" "+str(ySize))
