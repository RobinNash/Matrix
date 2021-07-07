## matrix ##
## June, 2021 ##
## By Robin Nash ##
'''
This module contains Matrix, Vector, and RowOp classes.
Matrix objects store entries as fractions and implement matrix operations.
Matrix also does more like RREF function implements Gaussian elimination/row reduction to return
a matrix in reduced row echelon form.

Vector is a subclass of Matrix that implements vector operations.

RowOp is a class that breaks down a row operation string to make performing row operations
simpler. For example, if the user wants to add 3 * row 2 to row 1 of a 4x5 Matrix M,
they can simply pass "R1 + 3*R2" into M.rowop() and the string will be given meaning by the RowOp class.
'''
from fractions import Fraction as Frac
from math import *

# Some math functions for returning values in degrees
def dcos(x): return cos(radians(x))
def dtan(x): return tan(radians(x))
def dsin(x): return sin(radians(x))
def dacos(x): return degrees(acos(x))
def dasin(x): return degrees(asin(x))
def datan(x): return degrees(atan(x))



class Matrix(list):
    '''2D matrix object
    simplements matrix operations and functions
    Matrix can be initialized by passing a 2D array
    or a string of form "a b c\nd e f" where a-c are entries of
    first row and d-f are entries of second row
    '''
    def __init__(self,matrix):
        super().__init__(matrix)
        # make each row a Matrix
        if self and type(self[0]) == list:
            for i in range(len(self)):
                self[i] = Matrix([Frac(a) for a in self[i]])
        # initialize n,m if matrix is 2D
        if self and isinstance(self[0],list):
            self.m,self.n = self.mn()
    
    def __repr__(self):
        return str([[str(c) for c in row] for row in self]).replace("'",'')

    def __add__(self, matrix):
        '''Return self + matrix'''
        if not self.is_same_size(matrix): raise ValueError("Matrices not of compatable size")
        return Matrix([Matrix([self[i][j]+matrix[i][j] for j in range(len(self[0]))]) for i in range(len(self))])
    def __radd__(self, matrix):
        '''Return matrix + self'''
        return matrix.__add__(self)

    def __neg__(self):
        '''Return self where each element is negated'''
        return Matrix([Matrix([-x for x in row]) for row in self])
    def __sub__(self, matrix):
        '''Return self - matrix'''
        return self + -matrix

    def __rsub__(self,matrix):
        '''Return matrix - self'''
        return matrix - self
    def __mul__(self,value):
        '''return self*value (value can be matrix or constant)'''
        m,r = self.mn()
        if isinstance(value,(int,Frac)):
            return Matrix([Matrix([x*value for x in row]) for row in self])
        r2,n = value.mn()
        if r != r2: raise ValueError("Matrices of incompatable sizes")
        return Matrix([Matrix([sum([self[i][t]*value[t][j] for t in range(r)]) for j in range(n)]) for i in range(m)])
    def __rmul__(self, value):
        '''return value*self (value can be matrix or constant)'''
        if isinstance(value,(int,Frac)):
            return Matrix([Matrix([x*value for x in row]) for row in self])
        return value.__mul__(self)

    def __floordiv__(self, value):
        '''Return self where each element x is x//value, value is int or Fraction'''
        return Matrix([Matrix([x//value for x in row]) for row in self])
    def __div__(self, value):
        '''Return self where each element x is x/value, value is a constant'''
        return Matrix([Matrix([x/value for x in row]) for row in self])

    def __pow__(self, value):
        '''Return self**value'''
        # if value is less than 0, we have to invert first, but we'll worry about this later
        if value > 0:
            M = self.copy()
            for i in range(value-1):
                M = M*self
        return M
	
    def print(self):
        '''display formatted matrix to console'''
        # holds format specifier for each column of self
        form = [max([len(str(self[i][j])) for i in range(self.m)]) for j in range(self.n)]
        M = [[f"{str(self[i][j]):>{form[j]}s}" for j in range(self.n)] for i in range(self.m)]
        print(str(M).replace('], [',']\n [').replace("'",''))#.replace(',',''))

    def copy(self):
        ''' Return a 2 level copy of self'''
        return Matrix([Matrix([x for x in row]) for row in self])
    def is_same_size(self, matrix):
	'''return if self has the same number of rows and columns as matrix'''
        return self.mn() == matrix.mn()
    def mn(self):
        '''Return (row,columns) of self'''
        return len(self),len(self[0])
    
    def remove_col(self,c):
        '''return self with column c removed'''
        return Matrix([[self[r][i] for i in range(self.n) if i != c] for r in range(self.m)])
    def remove_row(self,r):
        '''return self with row r removed'''
        return Matrix([self[i] for i in range(self.m) if i != r])
    
    # Row operations
    def swap(self,r1,r2):
        '''r1 <-> r2'''
        M = self.copy()
        M[r1],M[r2] = M[r2],M[r1]
        return M
    def scale(self,r1,c):
        '''r1 <- c*r1'''
        M = self.copy()
        M[r1] = c*M[r1]
        return M
    def pivot(self,r1,r2,c = 1):
        '''r1 <- r1 + c*r2'''
        m,n = self.mn()
        M = self.copy()
        M[r1] = M[r1]+c*M[r2]
        return M
    def row_op(self,opst):
        '''return matrix with row operation object or string opst applied to self'''
        opst = RowOp(str(opst))
        if opst.op == 0: return self.swap(opst.r1,opst.r2)
        if opst.op == 1: return self.scale(opst.r1,opst.c)
        if opst.op == 2: return self.pivot(opst.r1,opst.r2,opst.c)
        
    def T(self):
        '''Return transpose of self'''
        return Matrix([[self[j][i] for j in range(len(self))] for i in range(len(self[0]))])


    def REF(self, get_row_ops = False):
        '''Return self in a row echelon form, and the list of row operations that reduce it
		if get_row_ops = True.'''
        # Sort rows by least amount of leading zeros
        def leading_zeros(row,n):
            '''return the number of leading zeros in a list/Vector'''
            return n if row==[] or row[0]!=0 else leading_zeros(row[1:],n+1)
        def get_sort(M,start=0, ops = []):
            '''return M (with rows sorted by number of leading zeros) and row ops'''
            if start == M.m:
                return (M, ops)
            M = M.copy()
            leads = [leading_zeros(row,0) for row in M]
            r2 = leads.index(min(leads[start+1:]+[leads[start]]),start)
            if r2 != start:
                M[start],M[r2] = M[r2],M[start]
                ops.append(RowOp(0,start,r2))
            return get_sort(M, start+1, ops)
        
##            return Matrix(M[:start] + sorted(M[start:],key = lambda row: leading_zeros(row,0))) # if row_ops not involved
        
        M, row_ops = get_sort(self)

        for r in range(M.m):
            lead = leading_zeros(M[r],0) #where the current row's leading 1 will be
            # if zero row, no ops necessary
            if lead == M.n:
                break
            # Transform row so lead is 1
            if M[r][lead] != 1:
                row_ops.append(RowOp(1,r,Frac(1,M[r][lead])))
                M = M.scale(r,Frac(1,M[r][lead]))
            # Remove entries below
            for r2 in range(r+1,M.m):
                if M[r2][lead] == 0:
                    break
                lead2 = leading_zeros(M[r2],0)
                row_ops.append(RowOp(2,r2,r,-M[r2][lead2]))
                M = M.pivot(r2,r,-M[r2][lead2])
                
            # Sort the below by leading zeros again
            M,row_ops = get_sort(M,r+1, row_ops)

        return M if not get_row_ops else (M, row_ops)

    def RREF(self, get_row_ops = False):
        '''return self in reduced row echelon form, and the list of row operations that reduce it
		if get_row_ops = True'''
        def leading_zeros(row,n):
            return n if row==[] or row[0]!=0 else leading_zeros(row[1:],n+1)
        # put it in REF
        M, row_ops = self.REF(True)
        
        leads = [leading_zeros(row,0) for row in M]
        for r in range(M.m):
            for c in range(leads[r]+1,M.n):
                if c in leads:
                    r2 = leads.index(c)
                    row_ops.append(RowOp(2,r,r2, Frac(-M[r][c],M[r2][c])))
                    M = M.pivot(r,r2, Frac(-M[r][c],M[r2][c]))

        return M if not get_row_ops else (M, row_ops)
        
        
        
    def tr(self):
        '''return trace of self'''
        return sum([self[i][i] for i in range(n)])
    
    def det(self):
        '''Return determinant of self if self is square'''
        m,n = self.mn()
        if n!=m:
            raise ValueError("This Matrix is not sqaure")
        if n == 1:
            return self[0][0]
        if n == 2:
            return self[0][0]*self[1][1] - self[0][1]*self[1][0]
        # row expansion
        return sum([self[0][j]*self.C(0,j) for j in range(n)])
    def M(self,i,j):
        '''return the Minor of self at i,j(i.e. det of the matrix of self with row i and col j removed)'''
        return Matrix([[self[ii][jj] for jj in range(self.n) if jj != j] for ii in range(self.m) if ii!=i]).det()
    def C(self,i,j):
        '''return the cofactor of self at i,j'''
        return (-1)**(i+j)*self.M(i,j)

    def adj(self):
        '''return the adjoint matrix of self'''
        return Matrix([[self.C(j,i) for j in range(self.n)] for i in range(self.m)])

    def inverse(self):
        '''return the inverse matrix of self if it exists'''
        return Frac(1,self.det()) * self.adj()

    def TA(self, x):
        '''return Matrix transformation of self*x where x is a Vector'''
        return self*x.col()

def I(n):
	'''Return an n x n identity matrix'''
    return Matrix([[(1 if i==j else 0) for j in range(n)  ] for i in range(n)])
    

def Elementary(n, op, *args):
    ''' Return elementary matrix where a row operation
	is performed on the identity matrix of size n.
	row is row number (from 1) for op to be performed on.
    op is op number 0-2 or 's','m','p'
    args following op contains either:
    0/'s' (r1,r2)   : r1 <-> r2
    1/'m' (r1,c)    : r1 <- r1*c
    2/'p' (r1,r2,c) : r1 <- r1 + c*r2'''
    if str(op) in 'smp': op = 'smp'.find(op)
    if op == 0:
        self = I(n).swap(*args[:2]) # b is constant not row in this case
    if op == 1:
        self = I(n).scale(*args[:2])
    if op == 2:
        self = I(n).pivot(*args)
    return self

def ElementaryOpst(n, opst):
    '''Return elementary matrix where a row operation
	is performed on the identity matrix of size n.
	opst is row op string. ex "R2*-3", "R2sR3", "R2 - 3/2R3". no spaces necessary'''
    opst = RowOp(str(opst))
    return Elementary2(n,opst.op,*opst.tuple())


class RowOp:
    '''Holds details about an elementary row operation to be performed on a matrix.
    These are descriptions corresponding to op numbers:
    0 - Swap: two row numbers to indicate these rows to be interchanged  (r1,r2) r1 <-> r2
    1 - Scale: a target row, then a constant to multiply that row by (r1,c) r1 <- r1*c
    2 - Pivot: a row number, a another row number, then a constant to
    add constant * second row to first (r1,r2,c)  r1 <- r1 + c*r2'''
    
    def __init__(self, *args):
        '''args can be opst which is row op string. 
        Examples: 
        "R2*-3" -> multiply each entry in row to by constant -3 
        "R2sR3" -> switch rows 2 and 3
        "R2 - 3/2R3" -> add -3/2 of each entry in row 3 to row 2
        spaces in format are optional
        args can be op number (0-2), then r1,c or r1,r2 or r1,r2,c based on the number'''
        if len(args) == 1:
            if type(args[0]) == str:
                self.init_opst(args[0])
            elif type(args[0]) == RowOp:
                self.init_opst(args[0].opst)
        else:
            args += (0,)
            self.op = op = args[0]
            self.r1 = args[1]
            self.r2,self.c = [(args[2],1), (None,args[2]), (args[2],args[3])][op]
            # assign self.opst
            self.reconstruct()

    def init_opst(self,opst):		
        self.opst = opst
        self.op = op = ['s' in opst, '*' in opst, True].index(True)
        
        opst = opst.replace(' ','')
        r1,r2,c = None,None,0

        if op == 0:
            r1,r2 = map(int,opst.replace('R','').split('s'))
        if op == 1:
            r1,c = opst[1:].split('*')
            r1 = int(r1)
            if '/' in c:
                a,b = map(int,c.split("/"))
                c = Frac(a,b)
            else:
                c = int(c)
        if op == 2:
            pm = '+-'[int('-' in opst)]
            r1 = int(opst[1:opst.find(pm)])
            r2 = int(opst[opst.rfind("R")+1:])
            c = opst[opst.find(pm):opst.rfind("R")]
            if '/' in c:
                a,b = map(int,c.split("/"))
                c = Frac(a,b)
            else:
                c = int(c+('1' if len(c)==1 else ''))
        self.r1 = r1 - 1
        self.r2 = (None if not r2 else r2 - 1)
        self.c = c
        self.reconstruct()
        
    def __repr__(self):
        return self.opst
    
    def reconstruct(self):
        '''sets self.opst based on op,r1,r2,c values'''
        pm = "+-"[int(self.c < 0)]
        r1,r2 = self.r1+1, (None if self.r2 == None else self.r2+1)
        self.opst = [f"R{r1}sR{r2}", f"R{r1} * {self.c}", f"R{r1} {pm} {abs(self.c)}R{r2}"][self.op]

    def tuple(self):
        '''return op as tuple of form (r1,c,None), (r1,r2,None), or (r1,r2,c) based on self.op'''
        return [(self.r1,self.r2,None), (self.r1,self.c,None), (self.r1,self.r2,self.c)][self.op]

    def invert(self):
        '''Return the inverse row operation string of self'''
        opst = RowOp(self.opst)
        if opst.op == 1: opst.c = Frac(1,opst.c)
        if opst.op == 2: opst.c = -opst.c
        opst.reconstruct()
        return opst

class Vector(list):
    def __init__(self, vector):
        # cast Vector to str
        if type(vector) == str:
            vector = list(map(Frac,vector.split()))
        # convert matrix of 1 col to vector
        if isinstance(vector[0], list):
            if len(vector[0]) == 1: vector = [row[0] for row in vector]
            else: vector = vector[0] # would be a matrix with one row
        super().__init__(vector)

    def __repr__(self):
        return str([str(c) for c in self]).replace("'",'')
    def __neg__(self):
        '''return -self'''
        return -1*self
    def __mul__(self, value):
        '''return self*value. value is a constant'''
        return Vector([a*value for a in self])
    def __rmul__(self, value):
        '''return value*self. value is a constant'''
        return Vector([a*value for a in self])
    def __add__(self, vector):
        '''return self+vector'''
        return Vector([self[i]+vector[i] for i in range(len(self))])
    def __sub__(self, vector):
        '''return self - vector'''
        return self + -1*vector
##    def __setitem__(self, key, value):
##        '''set self[key] to value'''
##        self[key] = Frac(value)
    def norm(self):
        '''return the norm (length) of self'''
        return sqrt(self.normsq())
    def normsq(self):
        '''return the norm^2 of self'''
        return sum([v**2 for v in self])
    def unit(self):
        '''return unit vector of self'''
        return (1/self.norm())*self

    def dot(self,vector):
        '''return dot product of self and vector'''
        return sum([self[i]*vector[i] for i in range(len(self))])	
    def angle(self,vector):
        '''return angle between two vectors in radians'''
        return acos( self.dot(vector) / (self.norm()*vector.norm()) )
    def dangle(self,vector):
        '''return angle between self and vector in degrees'''
        return degrees(self.angle(vector))
	
    def cross(self,vector):
        '''return self x vector'''
        M = Matrix([self,vector])
        return Vector([M.remove_col(0).det(),-M.remove_col(1).det(),M.remove_col(2).det()])
    def lagrange(self,vector):
        '''return length of self cross vector using lagrange identity'''
        return sqrt( self.norm()**2 * vector.norm()**2 - self.dot(vector)**2 )
	
    def proj_len(self, vector):
        ''''return the length of the projection of self onto vector; proj vector self'''
        return (self.dot(vector)/vector.norm())#Frac(self.dot(a),a.norm())
    def proj(self,vector):
        '''return projection of self onto vector; proj vector self'''
        return Vector([Frac(self.dot(vector),vector.normsq())*c for c in vector])
    
    def col(self):
        '''return self as a column matrix'''
        return Matrix([[a] for a in self])
    def row(self):
        '''return self as a row matrix'''
        return Matrix([self])
	
def mul_list(matrices):
    '''multiply each matrix in order'''
    M = matrices[0].copy()
    for E in matrices[1:]:
        M = M*E
    return M

def apply_ops(M, row_ops):
    '''return a matrix where a list of RowOps or opst are applied to M in order'''
    M = M.copy()
    for op in row_ops:
        M = M.row_op(op)
    return M
  
  
if __name__ == "__main__":


    #RREF
    A = Matrix([[0,3,-1,2,-5], [3,6,9,-3,15], [3,9,8,-1,10]])
    A.REF().print()
    A.RREF().print()
    print()
    A = Matrix([[0,0,0,2,1,9],[0,-2,-6,2,0,2],[0,2,6,-2,2,0],[0,3,9,2,2,19]])
    A.REF().print()
    A.RREF().print()
    
##    x = Vector([[0,2,0]])
##    y = Vector([[0,0,3]])
##    u = Vector(x-2*y)
##    v = Vector(2*x+3*y)
##    print(v.norm())
##    print(u.norm())
##    angle = v.angle(u)
##    print(angle)


##
##    optests = ["R2-3/2R3","R2sR3","R2*-3/4"]
##    for op in optests:
##        op = RowOp(op)
##        print(op,op.invert())


    
