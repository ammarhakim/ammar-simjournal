from pylab import *
import tables

def plotLines(X, q0, q1, c):
    for i in range(q0.shape[0]):
        plot([X[i], X[i+1]], [q0[i], q1[i]], c)
        plot([0.5*(X[i]+X[i+1])], [0.5*(q0[i]+q1[i])], c[0]+'.')

figure(1)

fh = tables.open_file("a1-two-stream_distf_0.h5")
q = fh.root.StructGridField.read()
q0 = q[:,0,0]
q1 = q[:,0,1]
X = linspace(0, 1, q0.shape[0]+1)
plotLines(X, q0, q1, 'r-')
        
fh = tables.open_file("a1-two-stream_distf_1.h5")
q = fh.root.StructGridField.read()
q0 = q[:,0,0]
q1 = q[:,0,1]
X = linspace(0, 1, q0.shape[0]+1)
plotLines(X, q0, q1, 'k-')
grid()

show()

