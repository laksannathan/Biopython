import scipy.io
%matplotlib inline
import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import string

d = scipy.io.loadmat('D.mat')

for key, value in d.items() :
    mat = value

#plt.imshow(mat, cmap='hot', interpolation='nearest')
#plt.show()

'''
We can calculate also the mean distance among any two H. sapiens and between Neanderthal and any modern human.

D_human=reshape((D(1:n,1:n)),1,n^2);
human_dist=mean(D_human(find(D_human>0)))

D_mix=[D(1:n,n+1); D(1:n,n+2)];
mixed_dist=mean(D_mix)
'''

#dist between any two human
human_dist = 206*[0]

for i in range(206):
  for j in range(206):
    if i!=j:
        human_dist[i] += mat[i][j]

print([h/205 for h in human_dist])

#dist between any two human and neas
nea1 = 0
nea2 = 0

for j in range(206):
      nea1 += mat[206][j]
      nea2 += mat[207][j]

print(nea1/206)
print(nea2/206)


A = np.array(mat)*10

dt = [('len', float)]
A = A.view(dt)

G = nx.from_numpy_matrix(A)
G = nx.relabel_nodes(G, dict(zip(range(len(G.nodes())),string.ascii_uppercase)))    

G = nx.drawing.nx_agraph.to_agraph(G)

#G.node_attr.update(color="red", style="filled")
G.node_attr.update(color="red")
G.edge_attr.update(color="blue", width="2.0")

G.draw('graph.png', format='png', prog='neato')
