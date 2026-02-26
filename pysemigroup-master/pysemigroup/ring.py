import numpy as np

class BuchiMatrix(object):
    def __repr__(self):
        s = ""
        for x in self.dimension[0]:
            s += "\n"
            for y in self.dimension[1]:
                s+= str(self.data[(x,y)])+" "
        return s
    def __hash__(self):
        if not hasattr(self, "_hash"):
            self._hash= hash(str(self.data))
        return self._hash

    def is_idempotent(self):
        return (self == (self*self))
    def diagonal(self):
        data = {}
        if not(self.dimension[0] == self.dimension[1]):
            raise ValueError("Diagonal are only for square matrix")
        for x in self.dimension[0]:
            data[(x,0)] = self.data[(x,x)]
        return BuchiMatrix((self.dimension[0],(0,)),data)
    def projection(self,proj):
        for x in self.data:
            if self.data[x] in proj:
                self.data[x] = proj[self.data[x]]

    def __init__(self,dimension,data):
        self.dimension = dimension
        self.data = data
        self._hash_data = {}

    def __mul__(self,other):
        if not (self.dimension[1] == other.dimension[0]):
            raise ValueError("Dimensions mismatch for matrix product")
        dimension=(self.dimension[0],other.dimension[1])
        data = {}
        for x in dimension[0]:
            for y in dimension[1]:
                data[(x,y)] = "-oo"
                for k in self.dimension[1]:
                    k = buchiMul[(self.data[(x,k)],other.data[(k,y)])]
                    data[(x,y)] = buchiAdd[(data[(x,y)],k)]
        return BuchiMatrix(dimension,data)
    def __eq__(self,other):
        return hash(self)==hash(other)

buchiAdd = {
    ("-oo","-oo"):"-oo",
    ("-oo",0):0,
    ("-oo",1):1,
    (0,"-oo"):0,
    (0,0):0,
    (0,1):1,
    (1,"-oo"):1,
    (1,0):1,
    (1,1):1
}
buchiMul = {
    ("-oo","-oo"):"-oo",
    ("-oo",0):"-oo",
    ("-oo",1):"-oo",
    (0,"-oo"):"-oo",
    (0,0):0,
    (0,1):1,
    (1,"-oo"):"-oo",
    (1,0):1,
    (1,1):1
}



class hash_matrix(np.matrixlib.defmatrix.matrix):
    def __hash__(self):
        if not hasattr(self, "_hash"):
            self._hash= hash(self.tostring())
        return self._hash
    def __eq__(self,other):
        return (hash(self)==hash(other))
    def is_idempotent(self):
        return (self == (self*self))
    def inversible(self):
        if np.linalg.det(self) == 0:
            return False
        else:
            return True
    def inverse(self):
        N =  np.linalg.inv(self).astype(self.dtype)
        return hash_matrix(N)

