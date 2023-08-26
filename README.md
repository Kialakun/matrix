# matrix

Linear Algebra for humans.

## Example
'''
a = matrix.M{
    {1, 2},
    {3, 4},
}

// get transpose
aT := a.T()

// get Inverse
aI := a.Inv()

// cross product
b := a.X(a.T()).X(a.Inv())
'''
PS: This was initially developed to help me code Kalman Filters in Golang. 
