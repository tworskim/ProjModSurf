import numpy as np
from cvxopt import solvers, matrix
from scipy import sparse
import math
np.set_printoptions(threshold=np.inf)



filepath = "testnex.off"
print("Lecture Fichier")
(verts, faces) = read_off(open(filepath))
len_edges = len_edges(verts,faces)
print("Transformation mesh")
edges = find_edges(faces, len_edges)


longueur = longueur(edges)
gamma = 1
beta = 0.001
lambd = 0.02 * longueur**2 * gamma
alpha = 0.1 * gamma



####################
#FILE READER
def read_off(file):
    if 'OFF' != file.readline().strip():
        raise ValueError('Not a valid OFF header')
    n_verts, n_faces, n_dontknow = tuple([int(s) for s in file.readline().strip().split(' ')])
    verts = []
    for i_vert in range(n_verts):
        verts.append([float(s) for s in file.readline().strip().split(' ')])
    faces = []
    for i_face in range(n_faces):
        faces.append([int(s) for s in file.readline().strip().split(' ')][1:])
    return verts, faces

def write_off(verts, faces, filepath):
    file = open(filepath,"w") 
    file.write("OFF\n")
    file.write(str(len(verts)) + " " + str(len(faces)) +" " + "0\n")
    for i in range(len(verts)):
        for j in range(3):
            file.write(str(verts[i][j]) + " ")
        file.write("\n")
    #file.writelines(verts)
    for i in range(len(faces)):
        file.write("3 ")
        for j in range(3):
            file.write(str(faces[i][j]) + " ")
        file.write("\n")
    file.close()
    return 

####################
#MESH TRANSFORM
def common_edge(face1, face2):
    intersect = np.intersect1d(face1[:3], face2[:3])
    if len(intersect) >= 2 :
        p1 = intersect[0]
        p3 = intersect[1]
        for i in range(3):
            if (face1[i%3] == p1):
                if face1[(i+1)%3] == p3:
                    p4= np.setxor1d(face1[:3], intersect)[0]
                    p2= np.setxor1d(face2[:3], intersect)[0]
                elif face1[(i-1)%3] == p3:
                    p2= np.setxor1d(face1[:3], intersect)[0]
                    p4= np.setxor1d(face2[:3], intersect)[0]
                else:
                    print("Orientation couldn't be solved")
        edge = np.array([p1,p2,p3,p4])
        if (len(edge) < 4):
            print("Erreur longueur edge:")
            #print(edge)
            #print(face1[:3])
            #print(face2[:3])
            return(True, edge)
        else:
            return(True, edge)
    else:
        return(False, np.empty)

def len_edges(verts, faces):
    #wtf sort this mess 
    return(len(faces) + len(verts)+2)



def param(edges):
    l = 0
    gam = 0
    for i in range (len_edges):
        coords = coord(edges[i])
        p12 = coords[1] - coords[0]
        p14 = coords[3] - coords[0]
        p32 = coords[1] - coords[2]
        p34 = coords[3] - coords[2]
        p13 = coords[2] - coords[0]
        l += norm(edges[3]-edges[1])
        gam += math.acos(np.dot(p13,p12)/(norm(p13)*norm(p12)))
    l = l/len_edges
    print("Mean length")
    print("Mean dihedral angle")
    print(l)
    return(l,gam)
    
def find_edges(faces, len_edges):

    faces = np.array(faces)
    z = np.zeros((len(faces),1), dtype=faces.dtype)
    faces = np.append(faces, z, axis=1)
    
    edges = np.zeros((len_edges,4), dtype=np.uint32)
    k = 0
    for i, currentFace in enumerate(faces):
        j = i + 1
        while (currentFace[3] <3) & (j < len(faces)):
            if faces[j][3] < 3:
                (bool, edge) = common_edge(currentFace, faces[j])
                if bool:
                    if len(edge) < 4:
                        print(i,j)
                        print(currentFace, faces[j])
                        print("\n")
                    else:
                        edges[k] = edge
                        currentFace[3] +=  1
                        faces[j][3] += 1
                        k +=1
                elif faces[j][3] >3:
                    print("This face has common edges with more than 3 faces")
                    print(j, faces[j])
            j += 1;
        
    return(edges)
 
def coord(edge):
    coords = np.zeros((4,3), dtype=np.float)
    coords[0] = verts[edge[0]]
    coords[1] = verts[edge[1]]
    coords[2] = verts[edge[2]]
    coords[3] = verts[edge[3]]
    return(coords)
    
def norm(point):
    return(math.sqrt(np.dot(point,point)))
    
    
#########################
#OPERATOR CALCULATION
    
def calcops(coords):
    thetas = np.zeros((len_edges,4), dtype=np.float)
    deltas = np.zeros((len_edges,2), dtype=np.float)
    
    D = np.zeros((len_edges,3), dtype=np.float)
    Dmatrix = np.zeros((3*len_edges,len_edges * 12), dtype=np.float)
    R = np.zeros((len_edges,3), dtype=np.float)
    
    for i in range(len_edges):
        coords = coord(edges[i])
        p12 = coords[1] - coords[0]
        p14 = coords[3] - coords[0]
        p32 = coords[1] - coords[2]
        p34 = coords[3] - coords[2]
        p13 = coords[2] - coords[0]
        p31 = -p13
        
        #calcul thetas
        thetas[i][0] = math.acos(np.dot(p13,p12)/(norm(p13)*norm(p12)))
        thetas[i][1] = math.acos(np.dot(p14,p12)/(norm(p14)*norm(p12)))
        thetas[i][2] = math.acos(np.dot(p32,p31)/(norm(p32)*norm(p31)))
        thetas[i][3] = math.acos(np.dot(p34,p31)/(norm(p34)*norm(p31)))
        
        #Calcul deltas
        deltas[i][0] = 0.5 * norm(p12) * norm(p13) * math.sin(thetas[i][0])
        deltas[i][1] = 0.5 * norm(p14) * norm(p13) * math.sin(thetas[i][1])
        
        #Calc D
        a0 = (deltas[i][0]* np.dot(p34,p13) + deltas[i][1]* np.dot(p31,-p32))/(norm(p13)**2*(deltas[i][0] + deltas[i][1]))
        a1 = deltas[i][1]/(deltas[i][0] + deltas[i][1])
        a2 = (deltas[i][0]* np.dot(-p14,p13) + deltas[i][1]* np.dot(p31,p12))/(norm(p13)**2*(deltas[i][0] + deltas[i][1]))
        a3 = deltas[i][0]/(deltas[i][0] + deltas[i][1])
        
        Dmatrix[i][12*i:12*(i+1)] = np.array([a0, 0, 0, a1, 0, 0, a2, 0, 0, a3, 0, 0])
        Dmatrix[i+1][12*i:12*(i+1)] = np.array([0, a0, 0, 0, a1, 0, 0, a2, 0, 0, a3, 0])
        Dmatrix[i+2][12*i:12*(i+1)] = np.array([0, 0, a0, 0, 0, a1, 0, 0, a2, 0, 0, a3])
        D[i] = a0*coords[0] + a1 * coords[1] + a2*coords[2] + a3*coords[3]
        
        r_aux = (p12 + p34)
        R[i] = np.dot(r_aux,r_aux)
    return(thetas, deltas, D, R, Dmatrix)

###########################
# PARAMETERS OPTIMIZATION
def optdelt(D):
    delt = np.zeros((len_edges,3), dtype=np.float)
    count = 0
    for i in range(len_edges):
        if norm(beta * D[i]) > lambd:
            delt[i] = D[i]  
            count += 1 
    print("Norm L0 of delt:" + str(count))
    return(delt)

def verttoedge (edges):
    verttoedge = np.zeros((12*len_edges,3*len(verts)), dtype=np.float)
    for i in range (len_edges):
        for j in range(4):
            pos = edges[i][j]
            verttoedge[12*i + 3*j][3*pos:3*pos+3] = np.array([1,0,0])
            verttoedge[12*i + 3*j+1][3*pos:3*pos+3] = np.array([0,1,0])
            verttoedge[12*i + 3*j+2][3*pos:3*pos+3] = np.array([0,0,1])
    return(verttoedge)
    
def Rmatrix ():
    Rmatrix = np.zeros((3*len_edges, 12*len_edges), dtype=np.float)
    for i in range (len_edges):
        Rmatrix[3*i][12*i:12*(i+1)]   = np.array([1, 0, 0,-1, 0, 0, 1, 0, 0,-1, 0, 0])
        Rmatrix[3*i+1][12*i:12*(i+1)] = np.array([0, 1, 0, 0,-1, 0, 0, 1, 0, 0,-1, 0])
        Rmatrix[3*i+2][12*i:12*(i+1)] = np.array([0, 0, 1, 0, 0,-1, 0, 0, 1, 0, 0,-1])
    return(Rmatrix)


def transform1d(verts):
    pstar = np.zeros((3 *len(verts)), dtype = np.float)
    for i in range(len(verts)):
        pstar[3*i:3*i+3] = np.array([verts[i][0], verts[i][1], verts[i][2]])
    return(pstar)

def transformvect(p, verts):
    #verts = np.zeros((len(verts),3), dtype = np.float)
    for i in range(len(verts)):
        verts[i][0] =  float(p[3*i])
        verts[i][1] =  float(p[3*i + 1])
        verts[i][2] =  float(p[3*i + 2])
    return(verts)
    


print("Calcul Opérateurs")
(thetas,delta, D, R, Dmatrix) = calcops(edges)
delt = optdelt(D)

print("Calcul Matrices")
print("Calcul unique R et P")
r= Rmatrix()
e= verttoedge(edges)
r = np.mat(r)
a = sparse.coo.coo_matrix(np.mat(e))
Rmat = sparse.coo.coo_matrix(r) * a

Rmatsym = Rmat.transpose() * Rmat
Pconstpart = sparse.coo.coo_matrix(np.eye(3*len(verts)))
pstar = transform1d(verts)
pstar = np.mat(pstar)
len(pstar)
pstar = pstar.T
len(pstar)

step = 1
while (beta < 1000):

    print("Calcul Opérateurs")
    (thetas,delta, D, R, Dmatrix) = calcops(edges)
    deltvect = optdelt(D)
    print("Calcul récurrent D")
    Dmatrix = sparse.coo.coo_matrix(np.mat(Dmatrix))
    #Dmatrix = sparse.csc_matrix(Dmatrix)
    DM = Dmatrix * a
    DMT = DM.transpose()
    DMsym =  DMT * DM
    #DM = DM.todense()
    
    print("Calcul récurrent P")
    P = Pconstpart + beta * (DMsym) + alpha * Rmatsym
    P = P.todense()
    delt = matrix(transform1d(deltvect))
    print("Calcul récurrent Q")
    Qaux = beta *(DMT.todense().dot(delt))
    Q = - matrix(pstar) - Qaux
    print("opt")
    Psolv = matrix(P)
    Qsolv = matrix(Q)
    print(len(Qsolv))
    print(len(Qsolv.T))
    optres = solvers.qp(Psolv, Qsolv)
    p = optres['x']
    p = np.mat(p)
    pstar = p
    verts = transformvect(p, verts)
    alpha = 0.5*alpha
    beta = math.sqrt(2) * beta
    print(beta)
    print(step)
    write_off(verts, faces, "result"+ str(step)+ ".off")
    step +=1

print("fin")
print("Creation fichier off result")
write_off(verts, faces, "result.off")
write_off(transformvect(verts, pstar), faces, "testytest.off")






test  = np.matrix([[0,0,1], [0,0,1], [0,1,0]])
test2 = np.matrix([[0,0,1], [0,0,1], [0,1,0]])
a = sparse.coo.coo_matrix(test) * sparse.coo.coo_matrix(test2)

print(a.todense())
print(test* test2)







