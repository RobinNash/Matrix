# Testing RowOp object of matrix.py module #
# and Matrix functions implementing RowOp  #
# By Robin Nash
from matrix import *

A = Matrix('''1 0 1 1
1 1 1 4
2 1 1 2''')

# init
ops = [RowOp(0,1,2),RowOp(1,1,2),RowOp(2,1,2,3)]
for op in ops:
    print(op)
print()

# testing invert function of RowOp
optests = ["R2-3/2R3","R2sR3","R2*-3/4"]
for op in optests:
    op = RowOp(op)
    print(op,op.invert(),sep=', ')
print()

# RowOp on matrix
optests = ["R2+R1","R1sR3","R3*-2"]
for opst in optests:
    A.row_op(opst).print()
    print()

    

# test if Row Ops from RREF work
# and test apply_ops function

A.RREF().print()

M, row_ops = A.RREF(True)
print(row_ops)
apply_ops(A, row_ops).print()

