import numpy as np
nodes = dict()

class TrussNode:
  def __init__(node, id, X, Y):
    node.id = id               # Node Id
    node.X = X                 # node x coordinates
    node.Y = Y                 # node y coordinates
    node.rest = [0, 0]         # Node restraints [X, Y]
    node.force = [0, 0]        # Node forces positive direction Up, Right [X, Y]
    node.dof = [-1, -1] 
    nodes[id] = node           # including node Trussnode objecting in nodes dictionary using ID as refrence
 

elements = dict()

class Truss_Element:
  def __init__(elm, id, conn, E, A):
    elm.id = id                # Element ID
    elm.conn = conn            # Element conectivity data, ith note and jth node
    elm.EA = E*A 
    elements[id] = elm         # Including  Truss_element Object in Elements dictionary

  def Lx(element):
    ni = element.conn[0]        # first entry of the conn array
    nj = element.conn[1]        # second entry in the conn array
    x1 = nodes[ni].X            # Fetching the x coordinate of node i of the element
    x2 = nodes[nj].X            # Fetching the x coordinate of node j of the element
    return x2 - x1              # Returning x projection of the element
        
  def Ly(element):       
    ni = element.conn[0]        # first entry of the conn array
    nj = element.conn[1]        # second entry in the conn array
    y1 = nodes[ni].Y            # Fetching the x coordinate of node i of the element
    y2 = nodes[nj].Y            # Fetching the x coordinate of node j of the element
    return y2 - y1              # Returning y projection of the element

  def L(element):               # Element Length Function Definition
    Lx = element.Lx()           # Fetching X projection
    Ly = element.Ly()           # Fetching Y projection
    return (Lx**2 + Ly**2)**0.5 # Calculation Element Length

  def codevector(elm):          # Defining the Code Vector of an element
    ni = nodes[elm.conn[0]]     # Fetching Element's ith node DOFs
    nj = nodes[elm.conn[1]]     # Fetching Element's jth node DOFs
    return ni.dof + nj.dof      # Combining both nodes DOFs to form the elements code vector | Listing them

  def K(row):                   # Consutrction K element Matrix
    EA = row.EA
    Lx = row.Lx()
    Ly = row.Ly()
    L = row.L()

    c = Lx / L                 # Cosine 
    s = Ly / L                 # Sine
    cc = c * c                 # Consine squared
    ss = s * s                 # Sine Squared
    cs = c * s                 # Cosine * sine

    return EA/L * np.asarray( [ [cc,   cs,  -cc,  -cs ],
                                [cs,   ss,  -cs,  -ss ],
                                [-cc, -cs,   cc,   cs ],
                                [-cs, -ss,   cs,   ss ] ] )
