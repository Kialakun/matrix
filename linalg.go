package matrix

import (
	"log"
	"math"
)

// LU (Lower Upper Decomposition)
// ------------------------------
// This function decomposes a matrix to its LU components.
// It takes in two matrices, A (the coeffecient matrix of the augmented matrix)
// and B (the result matrix of the augmented matrix).
func LU(A, B M) (L, U M) {
	dimx, dimy := A.GetDimensions()
	if dimx != dimy {
		panic("A is not a square matrix")
	}
	L = zeroes(dimx, dimy)
	for i := 0; i < dimy-1; i += 1 {
		if A[i][i] == 0 { // nth (ie:first) row element
			panic("math error, first row matrix has zero value, cannot divide by zero")
		}
		L[i][i] = 1                        // set U element
		for j := i + 1; j < dimy; j += 1 { // j is the (n + 1)th (ie:second) row element
			ratio := A[j][i] / A[i][i]
			L[j][i] = ratio
			for k := 0; k < dimx; k += 1 {
				A[j][k] = A[j][k] - ratio*A[i][k] // calculate upper trianglar matrix element
			}
		}
	}
	return L, U
}

// Calculate the Inverse of a Matrix usig the Gauss Jordan Method
func (a M) Inv() M {
	dimx, dimy := a.GetDimensions()
	if dimx == dimy && dimx == 1 {
		// a is a scalar. The inverse is just the reciprocal of the scalar
		return M{{1 / a[0][0]}}
	}
	if dimx != dimy {
		panic("error: non square matrix detected")
	}
	n := len(a)
	// augment the  matrix with "a" diagonal identity matrix the same dimensions as "a"
	for i := 0; i < n; i++ {
		for j := 0; j < n; j++ {
			if i == j {
				a[i] = append(a[i], 1)
			} else {
				a[i] = append(a[i], 0)
			}
		}
	}
	// Apply Gauss Jordan Elimination on Augmented Matrix (A):
	for i := 0; i < n; i++ {
		if a[i][i] == 0 {
			panic("mathemtical error, ai,i division by zero")
		}
		for j := 0; j < n; j++ {
			if i != j {
				ratio := a[j][i] / a[i][i]
				for k := 0; k < 2*n; k++ {
					a[j][k] = a[j][k] - ratio*a[i][k]
				}
			}
		}
	}
	// Row Operation to Convert Principal Diagonal to 1.
	for i := 0; i < n; i++ {
		for j := n; j < 2*n; j++ {
			a[i][j] = a[i][j] / a[i][i]
		}
	}
	// slice matrix
	for i := 0; i < n; i++ {
		a[i] = a[i][n:]
	}
	return a
}

// Cholesky Decomposition algorithm, takes in A and returns U where A = U x U.T()
func (a M) CholeskyDecomposition() M {
	L := len(a)
	u := zeroes(L, L)
	for i := 0; i < L; i++ {
		for j := 0; j < i+1; j++ {
			sum := 0.
			if i == j {
				for k := 0; k < j; k++ {
					sum += u[j][k] * u[j][k]
				}
				u[j][j] = math.Sqrt(a[j][j] - sum)
				if math.IsNaN(u[j][j]) {
					panic("NaN detected")
				}
			} else {
				for k := 0; k < j; k++ {
					sum += u[i][k] * u[j][k]
				}
				if u[j][j] > 0 {
					u[i][j] = (a[i][j] - sum) / u[j][j]
				} else {
					log.Println(u[j][j])
					log.Println(a[i][j])
					log.Println(sum)
					log.Println(a)
					log.Println("WARNING!! Dvision by zero. matrix.linalg.CholeskyDecomposition")
				}
			}
		}
	}
	return u
}

func (a M) HadamardProduct(b M) (R M) {
	dx, dy := a.GetDimensions()
	bx, by := b.GetDimensions()
	if dx != bx || dy != by {
		log.Fatal("matrix mismatch, cannot take hammard product of matrices of different sizes")
		return
	}
	R = Zeroes(dx, dy)
	for y := 0; y < dy; y++ {
		for x := 0; x < dx; x++ {
			R[y][x] = a[y][x] * b[y][x]
		}
	}
	return R
}
