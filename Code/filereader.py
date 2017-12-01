import numpy as np
from scipy import sparse
from joblib import Parallel, delayed
import math
import time
np.set_printoptions(threshold=np.inf)

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
        p13 = coords[2] - coords[0]
        l += norm(p13)
        gam += math.acos(np.dot(p14,p12)/(norm(p14)*norm(p12)))
    l = l/len_edges
    gam = gam/len_edges
    print("Mean length")
    print("Mean dihedral angle")
    print(l)
    print(gam)
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
    return(math.sqrt(np.dot(point,point.T)))
    
    
#########################
#OPERATOR CALCULATION
#timer = 0
def calcops():

    dcol = []
    drow = []
    ddata = []
    nrow = 3*len_edges
    ncol = len_edges * 12
    Parallel(n_jobs = 5) (delayed(subcal)(i) for i in range(len_edges))
        
    Dmatrix = sparse.coo_matrix((ddata, (drow,dcol)), shape=(nrow, ncol), dtype=np.float)
    return(Dmatrix)


def subcal(i):
    coords = coord(edges[i])
        #start = time.time()
    c1 = coords[0]
    c2 = coords[1]
    c3 = coords[2]
    c4 = coords[2]
    #timer += time.time() - start
    p12 = c2 - c1
    p14 = c4 - c1
    p32 = c3 - c2
    p34 = c4 - c3
    p13 = c3 - c1
#       p12 = p12.T
#       p14 = p14.T
#       p32 = p32.T
#       p34 = p34.T
#       p13 = p13.T
    p31 = -p13

#calcul thetas
    theta0 = math.acos(np.dot(p13,p12.T)/(norm(p13)*norm(p12)))
    theta1 = math.acos(np.dot(p14,p12.T)/(norm(p14)*norm(p12)))
        #thetas[i][2] = math.acos(np.dot(p32,p31)/(norm(p32)*norm(p31)))
        #thetas[i][3] = math.acos(np.dot(p34,p31)/(norm(p34)*norm(p31)))
        
        #Calcul deltas
    delta0 = 0.5 * norm(p12) * norm(p13) * math.sin(theta0)
    delta1 = 0.5 * norm(p14) * norm(p13) * math.sin(theta1)
        
        #Calc D
    a0 = (delta0* np.dot(p34,p13.T) + delta1* np.dot(p31,-p32.T))/(norm(p13)**2*(delta0 + delta1))
    a1 = delta1/(delta0 + delta1)
    a2 = (delta0 * np.dot(-p14,p13.T) + delta1* np.dot(p31,p12.T))/(norm(p13)**2*(delta0 + delta1))
    a3 = delta0/(delta0 + delta1)

    drow += 12* [i] + 12* [i+1] + 12* [i+2]
    dcol +=  3 * list(range(12*i, 12*(i + 1)))
    ddata += [a0, 0, 0, a1, 0, 0, a2, 0, 0, a3, 0, 0] + [0, a0, 0, 0, a1, 0, 0, a2, 0, 0, a3, 0] + [0, 0, a0, 0, 0, a1, 0, 0, a2, 0, 0, a3]
    return()
###########################
# PARAMETERS OPTIMIZATION
def optdelt(D):
    delt = np.zeros((3*len_edges,1), dtype=np.float)
    #count = 0
    for i in range(len_edges):
        dp = D[3*i:3*(i+1)]
        if beta * norm(dp.T)  > lambd:
            delt[3*i:3*(i+1)] = dp
            #count += 1 
    #print("Norm L0 of delt:" + str(count))
    return(delt)

def verttoedge (edges):
    nrow= 12*len_edges
    ncol=3*len(verts)
    col = []
    row = []
    val = []
    for i in range (len_edges):
        for j in range(4):
            pos = edges[i][j]
            col += [3*pos, 3*pos +1, 3*pos +2]
            row += [12*i + 3*j, 12*i + 3*j+1, 12*i + 3*j+2]
            val += [1,1,1]
    verttoedge = sparse.coo_matrix((val, (row,col)), shape=(nrow, ncol), dtype=np.uint32)
    return(verttoedge)
    

def Rmatrix ():
    col = []
    row = []
    val = []
    nrow= 3*len_edges
    ncol=12*len_edges
    for i in range (len_edges):
        row += 12 * [3*i] +  12 * [3*i+1] + 12 * [3*i+2]
        col += 3 * list(range(12*i, 12*(i + 1)))
        val += [1, 0, 0,-1, 0, 0, 1, 0, 0,-1, 0, 0] + [0, 1, 0, 0,-1, 0, 0, 1, 0, 0,-1, 0] + [0, 0, 1, 0, 0,-1, 0, 0, 1, 0, 0,-1]
    Rmatrix = sparse.coo_matrix((val, (row,col)), shape=(nrow, ncol), dtype=np.float)
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
    
filepath = "input.off"
print("Lecture Fichier")
(verts, faces) = read_off(open(filepath))
len_edges = len_edges(verts,faces)

print("Transformation mesh")
edges = find_edges(faces, len_edges)
(longueur, gamma) = param(edges)
beta = 0.001
lambd = 0.02 * longueur**2 * gamma
alpha = 0.1 * gamma

print("Calcul Opérateurs")
#Dmatrix = calcops(edges)

print("Calcul Matrices")
print("Calcul unique R et P")
r= Rmatrix()
e= verttoedge(edges)
Rmat = r * e
Rmatsym = Rmat.transpose() * Rmat
Pconstpart = sparse.coo.coo_matrix(np.eye(3*len(verts)))
pstar = transform1d(verts)
pstar = np.mat(pstar)
pstar = pstar.T

step = 1

###
#Time mesurement variable 
t0 = 0
t1 = 0
t2 = 0
t3 = 0
t4 = 0
t5 = 0
t6 = 0
t01 = 0
#TODO factoriser sparse.coo.coo_matrix(pstar)
while (beta < 1000):

    #print("Calcul Opérateurs")
    start =time.time()
    coords = e * sparse.coo.coo_matrix(pstar)
    coords = coords.todense()
    t0 += time.time() - start
    start =time.time()
    Dmatrix  = calcops()
    t01 += time.time() - start
    #print("Calcul récurrent D")
    start =time.time()
    
    DM = Dmatrix * e
    DMT = DM.transpose()
    DMsym =  DMT * DM
    
    t1 += time.time() - start
                 
    #print("Optim delta")
    start = time.time()    
    DMP = DM * sparse.coo.coo_matrix(pstar)
    DMP = DMP.todense()
    delt= optdelt(DMP)
    t2 += time.time() - start
    #print("Calcul")
    start = time.time()
    P = Pconstpart + beta * (DMsym) + alpha * Rmatsym
    t3 += time.time() - start
    #♠delt = matrix(transform1d(deltvect))
    start = time.time()
    Qaux = beta *(DMT * delt)
    Q =  sparse.coo.coo_matrix(pstar) + Qaux
    t4 += time.time() - start
    #print("opt")
    start = time.time()
    p = sparse.linalg.spsolve(P, Q)
    t5 += time.time() - start
    #optres = solvers.qp(Psolv, Qsolv)
    #p = optres['x']
    start = time.time() 
    p = np.mat(p).T
    tmp1 = DMP - delt
    #print("Norme de beta *  D(p) - delt: " + str(beta * tmp1.T*tmp1))
    tmp2 = pstar - p
    #print("Norme de p-p*: " + str(tmp2.T * tmp2))
    pstar = p
    verts = transformvect(p, verts)
    #alpha = 0.5*alpha
    beta = math.sqrt(2) * beta
    
    #print(beta)
    print(step)
    #write_off(verts, faces, "result"+ str(step)+ ".off")
    step +=1
    t6 += time.time() - start

print("fin")
print("Creation fichier off result")
write_off(verts, faces, "result.off")
print("Temps")
print(t0)
print(t01)
print(t1)
print(t2)
print(t3)
print(t4)
print(t5)
print(t6)
#
#len_egdes








