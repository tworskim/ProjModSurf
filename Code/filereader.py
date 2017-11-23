import numpy as np
np.set_printoptions(threshold=np.inf)


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

def common_edge(face1, face2):
    intersect = np.intersect1d(face1[:3], face2[:3])
    if len(intersect) >= 2 :
        #print(np.union1d(face1[:3], face2[:3]))
        edge = np.union1d(face1[:3], face2[:3])
        if (len(edge) < 4):
            print("Erreur longueur edge\n")
            print(edge)
            print(face1[:3])
            print(face2[:3])
            return(True, edge)
        else:
            return(True, edge)
    else:
        return(False, np.empty)

def len_edges(verts, faces):
    return(len(faces) + len(verts) -2)

def find_edges(faces, len_edges):

    faces = np.array(faces)
    z = np.zeros((len(faces),1), dtype=faces.dtype)
    faces = np.append(faces, z, axis=1)
    
    edges = np.zeros((len_edges,4), dtype=np.float32)
    k = 0
    for i, currentFace in enumerate(faces):
        j = i + 1
        while (currentFace[3] <3) & (j < len(faces)):
            (bool, edge) = common_edge(currentFace, faces[j])
            if bool:
                print(edge)
                print(i,j)
                edges[k] = edge
                currentFace[3] +=  1
                faces[j] += 1
            j += 1;
    return(edges)
    
    
filepath = "../Data/homer.off"
(verts, faces) = read_off(open(filepath))
len_edges = len_edges(verts,faces)

edges = find_edges(faces, len_edges)
print(edges)