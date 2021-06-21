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
        print(str(self).replace('], [',']\n ['))#.replace(',',''))

    def copy(self):
        ''' Return a 2 level copy of self'''
        return Matrix([Matrix([x for x in row]) for row in self])
    def is_same_size(self, matrix):
	'''return if self has the same number of rows and columns as matrix'''
        return self.mn() == matrix.mn()
    def mn(self):
        '''Return (row,columns) of self'''
        return len(self),len(self[0])
    
    # Row operations
    def mul(self,target,c):
        '''target <- c*target'''
        M = self.copy()
        M[target] = [c*a for a in M[target]]
        return M
    def switch(self,r1,r2):
        '''r1 <-> r2'''
        M = self.copy()
        M[r1],M[r2] = M[r2],M[r1]
        return M
    def add(self,target,r1,r2,c = 1):
        '''target <- r1 + c*r2'''
        m,n = self.mn()
        return Matrix([[self[r][col] if r != target else self[r1][col]+c*self[r2][col] for col in range(n) ] for r in range(m)])
    def row_op(self,opst):
        '''return matrix with row operation object or string opst applied to self'''
        opst = RowOp(str(opst))
        if opst.op == 0: return self.mul(opst.r1,opst.c)
        if opst.op == 1: return self.switch(opst.r1,opst.r2)
        if opst.op == 2: return self.add(opst.r1,opst.r1,opst.r2,opst.c)
        
    def T(self):
        '''Return transpose of self'''
        return Matrix([Matrix([self[j][i] for j in range(len(self))]) for i in range(len(self[0]))])

    def REF(self):
        '''Return self in row echelon form'''
        # Sort rows by least amount of leading zeros
        def leading_zeros(row,n):
            return n if row==[] or row[0]!=0 else leading_zeros(row[1:],n+1)
        def get_sort(M,start=0):
            return Matrix(M[:start] + sorted(M[start:],key = lambda row: leading_zeros(row,0)))
        M = get_sort(self)

        for r in range(M.m):
            lead = leading_zeros(M[r],0) #where the current row's leading 1 will be
            # if zero row, no ops necessary
            if lead == M.n:
                break
            # Transform row so lead is 1
            if M[r][lead] != 1:
                M = M.mul(r,Frac(1,M[r][lead]))
            # Remove entries below
            for r2 in range(r+1,M.m):
                if M[r2][lead] == 0:
                    break
                lead2 = leading_zeros(M[r2],0)
                M[r2][lead2]
                M = M.add(r2,r2,r,-M[r2][lead2])
            # Sort the below by leading zeros again
            M = get_sort(M,r+1)

        return M

    def RREF(self):
        '''return self in reduced row echelon form'''
        def leading_zeros(row,n):
            return n if row==[] or row[0]!=0 else leading_zeros(row[1:],n+1)
        # put it in REF
        M = self.REF()
        leads = [leading_zeros(row,0) for row in M]
        for r in range(M.m):
            for c in range(leads[r]+1,M.n):
                if c in leads:
                    r2 = leads.index(c)
                    M = M.add(r,r,r2, Frac(-M[r][c],M[r2][c]))
        return M
        
        
    def tr(self):
        '''return trace of self'''
        return sum([self[i][i] for i in range(n)])
    
    def det(self):
        '''Return determinant of self if self is square'''
        m,n = self.mn()
        if n!=m:
            raise ValueError("Matrix not sqaure")
        if n == 1:
            return self[0][0]
        if n == 2:
            return self[0][0]*self[1][1] - self[0][1]*self[1][0]
        # row expansion
        return sum([self[0][j]*self.C(0,j) for j in range(n)])
    def M(self,i,j):
        '''return the det of the matrix of self with row i and col j removed'''
        m,n = self.mn()
        return Matrix([[self[ii][jj] for jj in range(n) if jj != j] for ii in range(m) if ii!=i]).det()
    def C(self,i,j):
        '''return the cofactor of self at i,j'''
        return (-1)**(i+j)*self.M(i,j)

    def adj(self):
        '''return the adjoint matrix of self'''
        m,n = self.mn()
        return Matrix([[self.C(j,i) for j in range(n)] for i in range(m)])

    def inverse(self):
		'''return the inverse matrix of self if it exists'''
        return Frac(1,self.det()) * self.adj()

    def TA(self, x):
        '''return Matrix transformation of self*x'''
        return self*x


def I(n):
	'''Return an n x n identity matrix'''
    return Matrix([[(1 if i==j else 0) for j in range(n)  ] for i in range(n)])
    
def Elementary(n, op, row_op_tuple):
    '''row is row number (from 0) for op to be performed on.
    op is op number 0-2 or 'm','s','a'
    row op tuple contains either:
    - 0: a target row, then an int or frac to multiply that row by (r1,c),
    - 1: a row number then another row number to indicate these
    rows to be interchanged  (r1,r2)
    - 2: a row number, a constant, then another row number to
    add constant * second row to first (r1,c,r2)  r1 <- r1 + c*r2'''
    if str(op) in 'msa': op = 'msa'.find(op)
    a,b,c = (row_op_tuple+(0,))[:3]
    if op == 0:
        self = I(n).mul(a,b)
    if op == 1:
        self = I(n).switch(a,b)
    if op == 2:
        self = I(n).add(a,a,c,b)
    return self
          
El = Elementary

def Elementary2(n, op, a, b, c = 1):
    '''row is row number (from 1) for op to be performed on.
    op is op number 0-2 or 'm','s','a'
    row op tuple contains either: (r1,c), (r1,r2), (r1,r2,c)  r1 <- r1 + c*r2'''
    if str(op) in 'msa': op = 'msa'.find(op)
    a-=1
    b-=1
    if op == 0:
        self = I(n).mul(a,b+1) # b is constant not row in this case
    if op == 1:
        self = I(n).switch(a,b)
    if op == 2:
        self = I(n).add(a,a,b,c)
    return self
E2 = Elementary2

def Elementary3(n, opst):
    '''opst is row op string. ex "R2*-3", "R2sR3", "R2 - 3/2R3". no spaces necessary'''
    opst = RowOp(str(opst))
    a,b,c = opst.tuple()
    return Elementary2(n,opst.op,a,b,c)

E3 = Elementary3

class RowOp:
    def __init__(self,opst):
        '''opst is row op string. 
		Examples: 
		"R2*-3" -> multiply each entry in row to by constant -3 
		"R2sR3" -> switch rows 2 and 3
		"R2 - 3/2R3" -> add -3/2 of each entry in row 3 to row 2
		spaces in format are optional'''
		
        self.opst = opst
        self.op = op = ['*' in opst, 's' in opst, '*' not in opst and 's' not in opst].index(True)
        
        opst = opst.replace(' ','')
        r1,r2,c = None,None,0
        if op == 0:
            r1,c = opst[1:].split('*')
            r1 = int(r1)
            if '/' in c:
                a,b = map(int,c.split("/"))
                c = Frac(a,b)
            else:
                c = int(c)
        if op == 1:
            r1,r2 = map(int,opst.replace('R','').split('s'))
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
        self.r1 = r1
        self.r2 = r2
        self.c = c
        self.reconstruct()
    def __repr__(self):
        return self.opst
    
    def reconstruct(self):
        '''sets self.opst based on op,r1,r2,c values'''
        pm = "+-"[int(self.c < 0)]
        self.opst = [f"R{self.r1} * {self.c}", f"R{self.r1}sR{self.r2}", f"R{self.r1} {pm} {abs(self.c)}R{self.r2}"][self.op]

    def tuple(self):
        '''return op as tuple of form (r1,c,None), (r1,r2,None), or (r1,r2,c) based on self.op'''
        return [(self.r1,self.c,None), (self.r1,self.r2,None), (self.r1,self.r2,self.c)][self.op]

    def invert(self):
        '''Return the inverse row operation string of self'''
        opst = RowOp(self.opst)
        if opst.op == 0: opst.c = Frac(1,opst.c)
        if opst.op == 2: opst.c = -opst.c
        opst.reconstruct()
        return opst

class Vector(Matrix):
    def norm(self):
        '''return the norm (length) of self'''
        return sum([v**2 for v in self[0]])**0.5
    def dot(self,vector):
        '''return dot product of self and vector'''
        return (self*vector.T())[0][0]
    def angle(self,vector):
        '''return angle between two vectors in radians'''
        return acos( self.dot(vector) / (self.norm()*vector.norm()) )
    def dangle(self,vector):
        '''return angle between self and vector in degrees'''
        return degrees(self.angle(vector))

def mul_list(matrices):
    '''multiply each matrix in order'''
    M = matrices[0]
    for E in matrices[1:]:
        M = M*E
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


    
