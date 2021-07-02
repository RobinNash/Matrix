# Testing Vector Object of matrix.py module #
# By Robin Nash
from matrix import *

def allFrac(vector):
    '''returns if all entries in vector are of type Fraction'''
    return not int in [type(v) for v in Vector(vector)]

### init ###

# vector as list
a = Vector([2,2,0])
# row vector as matrix
b = Vector(Matrix([[1,0,0]]))
# column vector as matrix
c = Vector(Matrix([[1],[0],[1]]))
# str
d = Vector('9 8 1/2')
# Vector to Vector
e = Vector(a)
print(a,b,c,d,e)

# row/col functions
print(a.col())
print(a.row())

#norm
print(a.norm())
#unit
print(a.unit())

# proj_len
print(c.proj_len(b)) # 1
print(a.proj_len(d))

# proj
print(c.proj(b)) # should be b [1,0,0]
print(a.proj(b)) # [2,0,0]
print(a.proj(d))

