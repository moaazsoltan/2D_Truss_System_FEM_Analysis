#System input and definition begins at line 70
import numpy as np
from objects import nodes, TrussNode, Truss_Element,elements

np.set_printoptions(suppress=True)
np.set_printoptions(precision=4)

#Node Definition
TrussNode(id=1, X=0, Y=0)
TrussNode(id=2, X=3, Y=0)
TrussNode(id=3, X=3, Y=2)

#Conenctivity = Nodes [ith, jth]
Truss_Element(id=1, conn=[1,2], E=2000, A=1)
Truss_Element(id=2, conn=[2,3], E=4000, A=1)
Truss_Element(id=3, conn=[1,3], E=6000, A=1)

#Restraint Definition [X,Y]
nodes[1].rest = [1, 1]
nodes[2].rest = [0, 1]

#Forces Definition [X, Y]
nodes[3].force = [50, -40]


N = 0 #Unknown Displacements
M = 0 #Total DOFs

for id, node in nodes.items():
  if node.rest[0] == 0: # If the node has not restrained  in X direction then assign DOF
    node.dof[0] = M
    M += 1
  if node.rest[1] == 0: # If the node has not restrained  in Y direction then assign DOF
    node.dof[1] = M
    M += 1

N = M

for id, node in nodes.items():
  if node.rest[0] != 0:  # if restrained, then the dof should be equal the last value of M, then add +1 to M for the next DOF
    node.dof[0] = M
    M += 1
  if node.rest[1] != 0: # if restrained, then the dof should be equal the last value of M, then add +1 to M for the next DOF
    node.dof[1] = M
    M += 1

KSystem = np.zeros((M, M)) # K Sysmte, = M x M matrix size
U_Matrix = np.zeros(M)     # U Matrix is a Column Matrix of M size

for id, elm in elements.items(): #Forming the K System by superposing K Elements using Code Vectors
  KE = elm.K()
  cv = elm.codevector()
  for i in [0, 1, 2, 3]:
    for j in [0, 1, 2, 3]:
      ii = cv[i]
      jj = cv[j]
      KSystem[ii, jj] =  KSystem[ii, jj] + KE[i, j]

for id, node in nodes.items():                  #Forming P_System (P1 + P2)
  ii = node.dof[0] 
  U_Matrix[ii] = U_Matrix[ii] + node.force[0]
  ii = node.dof[1] 
  U_Matrix[ii] = U_Matrix[ii] + node.force[1]


"""Solving Algorithm"""
K11 = KSystem[0:N, 0:N]
P1 = U_Matrix[0:N]

Inverse_K11 = np.linalg.inv(K11)              #Inverse[K11]
U1 = Inverse_K11 @ P1                         # U1 = Inverse[K11] . P1
USystem = np.asarray(U1.tolist() + [0]*(M-N)) # Combining U1 and U2 to obtain U System. U2 is zero displacements due to restrained conditions. List Comprehension
PSystem = KSystem @ USystem                   #Reaction Forces



"""Printing Solution Values"""

#Print out Node & DOF data
print("N (Unknown Displacements):", N)
print("M (Total Dofs):", M)

print("Node DOFS:")
for id, node in nodes.items():
  print(id, node.dof)

print("Element Code Vectors")
for id, elm in elements.items():
  print(id, elm.codevector())

#Printing out K11, P1, U1, Usystem, Psystem
print("K11: \n", K11)
print("P1 \n", P1)
print("U1: \n", U1)
print("u_system \n:", USystem)
print("P_system: \n", PSystem)