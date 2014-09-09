from numpy import matrix
from numbers import Number

dotprod = lambda x,y : float(x.T * y)
distance = lambda x,y : (x-y).norm

class vector(matrix):
    def __new__(S, *args, **kwargs):
        try:
            x= matrix.__new__(S, *args, **kwargs)
        except TypeError:
            x= matrix.__new__(S, [[i] for i in args])

        if isinstance(x, vector): return x
        else: return vector(x.A)

    def __getattr__(self, name):
        if name == 'T':
            self.T=self.getT()#self.__dict__['T'] = T
        elif name in {'norm', 'norm_squared'}:
            #db print self
            try:
                self.norm_squared = dotprod(self,self)
            except ValueError:
                self.norm_squared = dotprod(self.T, self.T)
            self.norm = self.norm_squared **.5
        elif name in {'x','y','z'}:
            self.__dict__[name] = self.item( 'xyz'.find(name))
        else:
            raise AttributeError, 'vector object has no attribute '+name
                
        return self.__dict__[name]
        
    def __init__(self, *args, **kwargs):
        matrix.__init__(self, *args, **kwargs)

        if 1 not in self.shape:
            raise ValueError, 'input does not describe a vector.'
        
    def __repr__(self):
        return matrix.__repr__(self).replace('matrix','vector')

    def __set__(self, *args, **kwargs):
        raise TypeError, "'vector' object is immutable"

    def __mul__(self, right):
        x = super(vector, self).__mul__(right)
        if x.size == 1:
            return x.item()
        elif 1 in x.shape:
            return vector(x)
        else: return x

    def __rmul__(self, left):
        x = super(vector, self).__rmul__( left)
        if x.size == 1:
            return x.item()
        elif 1 in x.shape:
            return vector(x)
        else: return x
    

def xy_from_lengths(p0, p1, A, B):
    '''p0,p1: vectors (x,y) specifying two known points
       A: ||p2 - p0||
       B: ||p2 - p1||

       returns p2 as a vector (x,y).'''

    Asq = A**2
    Bsq = B**2
    v = p1 - p0
    C = v.norm
    Csq = v.norm_squared

    y = .5*(Asq + Csq - Bsq)/C
    xplus = (Asq - y**2)**.5
    xminus = -xplus

    T = matrix( (( v.y/C, v.x/C),
                 (-v.x/C, v.y/C)) )
    
    pplus = T*vector(xplus, y) + p0
    pminus = T*vector(xminus, y) + p0

    if pplus.x < 0: return pminus
    if pminus.x < 0: return pplus
    return pminus, pplus

