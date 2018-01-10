from numpy import pi, sin, cos, array, linspace, meshgrid, reshape, size, float32
from os import system


xSize = 1001
ySize = 1001
A, B = meshgrid(linspace(-2*pi, 2*pi, xSize), linspace(-2*pi, 2*pi, ySize))
C = sin(A)*cos(B)
X = array(reshape(A, (-1, size(A))), dtype=float32)
Y = array(reshape(B, (-1, size(B))), dtype=float32)
Z = array(reshape(C, (-1, size(C))), dtype=float32)
for a, s in zip((X, Y, Z), ("X", "Y", "Z")):
    with open("./res/%s.dat" % s, "wb") as f:
        f.write(a.tobytes())
system("build/vtkPlot $HOME/git/vtkPlot/res/X.dat $HOME/git/vtkPlot/res/Y.dat $HOME/git/vtkPlot/res/Z.dat %d %d" % (xSize, ySize))
