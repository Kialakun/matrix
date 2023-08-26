# matrix

Linear Algebra for humans.

## Example
```
a := matrix.M{
    {1, 2},
    {3, 4},
}

// get transpose
aT := a.T()

// get Inverse
aI := a.Inv()

// matrix multiplication
b := a.X(a.T()).X(a.Inv())

// LU Decomposition
l, u := a.LU()
```
PS: This was initially developed to help me code Kalman Filters in Golang. 
