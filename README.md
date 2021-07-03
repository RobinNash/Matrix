# Matrix
This repository holds a Python modules matrix.py. This module allows users to create matrix objects, perform matrix and vector operations, and easily perform elementary row operations on a matrix. This module is great for applications in linear algebra.

The Matrix class holds data of a matrix and can perform matrix operations.  
The Vector class holds data of a Vector and can perform vector operations.  
The RowOp class exists to store data on a elementary row operation. It can be represented by a readable string or data values.  

## How to Use
To create a Matrix object, cast Matrix to a 2D list or a string.
For example, the 2x3 matix:
```
1 2 3
4 5 6
```
can be initialized in the following ways
```python
A = Matrix([[1, 2, 3], [4, 5, 6]])
A = Matrix("1 2 3\n4 5 6")
A = Matrix('''1 2 3
4 5 6''')

# to print the matrix:
A.print()
```
A Vector object can be initialized by a 1D list or a string in similar ways:
```python
v = Vector([1,2,3])
v = Vector('1 2 3')
```

## Features
- Matrix
  - Supports scalar multiplication of a matrix, matrix multiplication, matrix addition/subtration, and matrix powers using arithmetic operators `* + - **`
  - REF, RREF functions to get a row echelon form or the reduced row echelon form of a matrix
  - Many transformative functions including `.T()`, `.inverse()`, `.adj()` to get the transpose, inverse, or adjoint of a matrix
  - Supports elementary row operations on a matrix using `.swap()`, `.scale()`, `.pivot()`, or `.row_op()` functions
  - scalar functions `.det()` and `.tr()` to get the determinant and trace of a square matrix
- Vector
  - Supports scalar/vector multiplication, vector addition/subtraction using arithmetic operators `* + -`
  - row, col functions to convert Vector into a row Matrix or column Matrix object
  - Includes lots of vector operations like cross product, dot product, norm
  - Includes `.proj(vector)` which returns the projection of the current vector on another vector
- RowOp
  - This class converts an elementary row operation string like `"R1 * 3"`, `"R2sR3"`, or `"R1 - 3/2R3"` into values that can be used by Matrix to perform the row operation
  - `.invert()` function returns the inverted version of the operation string. Ex. `"R1 * 3" -> "R1 * 1/3"`  
## License
[MIT](https://choosealicense.com/licenses/mit/)
