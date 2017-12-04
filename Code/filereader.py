import numpy as np
from scipy import sparse
#from joblib import Parallel, delayed
import math
import time
import os
import sys
np.set_printoptions(threshold=np.inf)


filename = "chat"


dirpath = "../Meshs/" + filename + "/" +  "input/"
dirpathout = "../Meshs/" + filename
if not os.path.exists(dirpathout):
    os.makedirs(dirpathout)
files = os.listdir(dirpath)

print("Select input std:")
for i,s in enumerate(files):
    print(str(i)+ ": " + s)
files = [files[int(input("Enter file id: "))]]
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
        rl = file.readline().strip().split(' ')
        if(len(rl) == 5):
            n_faces += 1
            rlbis = [rl[k] for k in [3,4,1]]
            faces.append([int(s) for s in rl[1:4]])
            faces.append([int(s) for s in rlbis])
        else:
            faces.append([int(s) for s in rl][1:4])
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

def write_edges(edges, filepath):
    file = open(filepath,"w")
    file.write(str(len_edges)+ "\n")
    for i in range(len(edges)):
        for j in range(4):
            file.write(str(edges[i][j]) + " ")
        file.write("\n")
    #file.writelines(verts)
    file.close()
    return 

def read_edges(filepath):
    file = open(filepath) 
    #len_edges = int(file.readline())
    edges = np.zeros((len_edges, 4), dtype = np.uint32)
    for i in range(len_edges):
        s = file.readline().strip().split(' ')
        for j in range(4):
            edges[i][j] = s[j]
    file.close()
    return edges

####################
#MESH TRANSFORM 

def lenedges(verts, faces): 
    return(len(faces) + len(verts)+ 2* (filename=="eight") - 2*(filename!="eight"))
    
def common_edge(face1, face2):
    intersect = np.intersect1d(face1, face2)
    if len(intersect) >= 2 :
        p1 = intersect[0]
        p3 = intersect[1]
#        for i in range(3):
#            if (face1[i%3] == p1):
#                if face1[(i+1)%3] == p3:
#                    p4= np.setxor1d(face1, intersect)[0]
#                    p2= np.setxor1d(face2, intersect)[0]
#                elif face1[(i-1)%3] == p3:
#                    p2= np.setxor1d(face1, intersect)[0]
#                    p4= np.setxor1d(face2, intersect)[0]
#                else:
#                    print("Orientation couldn't be solved")
        p4= np.setxor1d(face1, intersect)[0]
        p2= np.setxor1d(face2, intersect)[0]
        
        edge = np.array([p1,p2,p3,p4])
        return(True, edge)
    else:
        return(False, np.empty)

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
    #print("Mean length: " + str(l))
    #print("Mean dihedral angle: " + str(gam))
    return(l,gam)
    
    
def find_edges2(faces, len_edges): 
    edgelist = []
    edges = []
    for i, currentFace in enumerate(faces):
        edgeface = list_edges(currentFace[:3])
        edgelist, edges = matchedges(edgelist,edgeface, edges)
        sys.stdout.write("\r" + str(100 *int(len(edges)/len_edges)) + "%")
        sys.stdout.flush()
    return(edges)
    
def list_edges(currentFace):
    edge1 = currentFace[0:]
    edge2 = currentFace[1:] + [currentFace[0]]
    edge3 = [currentFace[2]] + currentFace[0:2]
    return([edge1,edge2,edge3])
    len(edges[2303])
def matchedges(edgelist, edgeface, edges):
    crctf = 0
    for i in range(3):
        edge = edgeface[i-crctf]
        (boolean, edgelist, edges) = compedges(edge, edgelist, edges)
        if boolean:
            edgeface.pop(i-crctf)
            crctf +=1
    edgelist += edgeface
    return(edgelist,edges)
    
def compedges(edge, edgelist, edges):
    i = len(edgelist) -1
    boolean = False
    while  (i >=0) & (boolean == False):
        #print(len(np.intersect1d(edge[:2], edgelist[i][:2])))
        if len(np.intersect1d(edge[:2], edgelist[i][:2]))>=2:
            edges.append([edge + [edgelist[i][2]]][0])
            edgelist.pop(i)
            boolean = True
        i += -1
    return(boolean, edgelist, edges)
    
def find_edges(faces, edgelist):
    

    faces = np.array(faces)
    z = np.zeros((len(faces),1), dtype=faces.dtype)
    faces = np.append(faces, z, axis=1)
    
    edges = np.zeros((len_edges,4), dtype=np.uint32)
    k = 0
    for i, currentFace in enumerate(faces):
        j = i + 1
        while (currentFace[3] < 3) & (j < len(faces)):
            if faces[j][3] < 3:
                face2 = faces[j]
                (bool, edge) = common_edge(currentFace[:3], face2[:3])
                if bool:
                    edges[k] = edge
                    currentFace[3] +=  1
                    face2[3] += 1
                    k +=1
                    faces[j][3] = face2[3]
            j += 1;
        sys.stdout.write("\r" + str(int(100*k/len_edges))+ "% " + str(k) + "/" + str(len_edges))
        sys.stdout.flush()
    if k < len_edges:
        edges = edges[:k]
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

def calcops():
    dcol = []
    drow = []
    ddata = []
    nrow = 3*len_edges
    ncol = len_edges * 12
    for i in range(len_edges) : 
        drow += 4* [i] + 4* [i+1] + 4* [i+2]
        for k in range(3):
            dcol +=  [12 * i + k, 12 * i +3 + k, 12 * i + 6 + k, 12 * i + 9 + k]
        dat = subcal(i)
        ddata += dat
    Dmatrix = sparse.coo_matrix((ddata, (drow,dcol)), shape=(nrow, ncol), dtype=np.float)
    return(Dmatrix)

def subcal(i):
    coords = coord(edges[i])
    
    c1 = coords[0]
    c2 = coords[1]
    c3 = coords[2]
    c4 = coords[2]
    
    p12 = c2 - c1
    p14 = c4 - c1
    p32 = c3 - c2
    p34 = c4 - c3
    p13 = c3 - c1
    
#calcul thetas
    
    np13 = norm(p13)
    np14 = norm(p14)
    np12 = norm(p12)
        #Calcul deltas
    
    delta0 = 0.5 * norm(p12) * np13 * math.sqrt(1 - (np.dot(p13,p12)/(np13*np12))**2)
    delta1 = 0.5 * norm(p14) * np13 * math.sqrt(1 - (np.dot(p14,p12)/(np14*np12))**2)
    sumd = delta0 + delta1
    
        #Calc D
    a0 = (delta0* np.dot(p34,p13) + delta1* np.dot(-p13,-p32))/(np13**2*sumd)
    a1 = delta1/(delta0 + delta1)
    a2 = (delta0 * np.dot(-p14,p13) + delta1* np.dot(-p13,p12))/(np13**2*sumd)
    a3 = delta0/sumd
    ddata = 3 * [a0, a1, a2, a3]
    return(ddata)
###########################
# PARAMETERS OPTIMIZATION   

def optdelt(D):
    delt = np.zeros((3*len_edges,1), dtype=np.float)
    count = 0
    for i in range(len_edges):
        dp = D[3*i:3*(i+1)]
        if beta * norm(dp.T)  >= lambd:
            delt[3*i:3*(i+1)] = dp
            count += 1 
    #sys.stdout.write("\r" + "Sharp conserved proportion:" +str(int(100*count/len_edges))+ "%")
    #sys.stdout.flush()
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
    
    
    
(verts, faces) = read_off(open(dirpath + "/" + "0.001.off"))
len_edges = lenedges(verts,faces)


######################################

#Gets edges if saved, calculate it if not
######################################
print("Edges recuperation for mesh: " + filename)
if os.path.isfile(dirpathout + "/" + "edges.off"):
    edges = read_edges(dirpathout + "/" + "edges.off")
    len_edges = len(edges)
    print("Edges were saved")
else:
    start = time.time()
    edges = find_edges(faces, len_edges)
    len_edges = len(edges)
    print("\nEdges saved for future use")
    write_edges(edges, dirpathout + "/" + "edges.off")
    timer = time.time() - start
    print("\nEdges recuperation time = " + str(int(timer/60)) + " m " + str(int(timer - 60 * int(timer/60))) + " s")


######################################

#Perform smoothing on all files from input folder
######################################
print("\nSharpness option:\n1: Paper behavior\n<1: Less sharp edges\n>1: more sharp edges")
sharpness = float(input("Enter sharpness: "))
print("\nSmoothness option:\n Float parameter\n0: No smoothing\n1: Paper behavior\nn: Smoothen n times harder \n2: Really smooth")
smoothness = float(input("Enter smoothness: "))
print("\nSpeed option:\n0: Paper behavior\n1: 2x faster\n2: 4x faster. Be careful\n3: 8x faster. Try if you want")
speed = int(input("Enter speed: "))

a = 0.5
alphaprev = 0.1 * smoothness
if(smoothness == 314159):
    a =1
    alphaprev = 0.1
if (speed == 0):
    b = math.sqrt(2)
elif (speed == 1):
    b = 2
    a = a**2
elif (speed == 2):
    b = 4
    a = a**4
elif (speed == 3):
    b = 16
    a = a**8
            
for i,file in enumerate(files):
    start = time.time()
    print("\nSmoothing mesh " + str(i+1) + "/" + str(len(files)))
    filepath = dirpath + "/" + file
    (verts, faces) = read_off(open(filepath))
    #len_edges = lenedges(verts,faces)
    (longueur, gamma) = param(edges)
    beta = 0.001
    if(sharpness <= 0):
        print("Sharpness must be > 0")
    if(sharpness > 0):
        lambd = 0.02 * longueur**2 * gamma / (sharpness)
    if(smoothness < 0):
        print("Smoothness must be >= 0")
    alpha = alphaprev *  gamma
    r= Rmatrix()
    e= verttoedge(edges)
    Rmat = r * e
    Rmatsym = Rmat.transpose() * Rmat
    Pconstpart = sparse.identity(3*len(verts))
    pstar = transform1d(verts)
    pstar = np.mat(pstar)
    pstar = pstar.T
    
    step = 1
    #TODO factoriser sparse.coo.coo_matrix(pstar)
    while (beta < 1000):
    
        coords = e * sparse.coo.coo_matrix(pstar)
        coords = coords.todense()
        Dmatrix= calcops()
        DM = Dmatrix * e
        DMT = DM.transpose()
        DMsym =  DMT * DM 
        DMP = DM * sparse.coo.coo_matrix(pstar)
        DMP = DMP.todense()
        delt= optdelt(DMP)
        P = Pconstpart + beta * (DMsym) + alpha * Rmatsym
        Qaux = beta *(DMT * delt)
        Q =  sparse.coo.coo_matrix(pstar) + Qaux
        (p,info) = sparse.linalg.cg(P, Q)
        if(info>0):
            print("Tolerance not achieved")
        p = np.mat(p).T
        pstar = p
        verts = transformvect(p, verts)

        alpha = a*alpha
        beta = b * beta
        
        sys.stdout.write("\r" + str(int(step*(2.5*2**(speed))))+ "%")
        sys.stdout.flush()
        #write_off(verts, faces, "result"+ str(step)+ ".off")
        step +=1
    
    
    filepathout = dirpathout + "/" + "result/" + "sh"+ str(sharpness) + "sm" + str(smoothness) + "sp" + str(speed) + "std" +file
    if not os.path.exists(dirpathout + "/" + "result/"):
        os.makedirs(dirpathout + "/" + "result")
    #print("\nCreation" + filepathout)
    write_off(verts, faces, filepathout)
    timer = time.time()-start
    print("\nSmoothing time = " + str(int(timer/60)) + " m " + str(int(timer - 60* int(timer/60))) + " s")
    write_off(verts, faces, "test.off")