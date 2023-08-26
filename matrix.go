// Matrix operations library for computation.

package matrix

import (
	"log"
	"math"
)

type Matrice interface {
	X(M) M                    // cross product
	T() M                     // Transpose
	Add(M) M                  // add a matrix
	Sub(M) M                  // subtract the matrix from a matrix
	Det2x2() float64          // calculate determinant
	Adj2x2() M                // returns the adjoint of the matrix
	Inv2x2() M                // returns the inverse of a 2x2 matrix
	Inv() M                   // returns the inverse of an arbitrary size square matrix
	Xscalar(float64) M        // returns the product of the matrix with a scalar
	CholeskyDecomposition() M // returns the cholesky factor of the matrix
}

type M [][]float64

// function returns a diagonal identity matrix of specified dimensions
func NewDiagonalMatrix(dimx, dimy int) M {
	result := zeroes(dimx, dimy)
	for i := 0; i < dimy; i++ {
		if i < dimx {
			result[i][i] = 1
		}
	}
	return result
}

// check two matrices if they have the same dimensions
func hasSameDimensions(a, b M) bool {
	adimx, adimy := a.GetDimensions()
	bdimx, bdimy := b.GetDimensions()
	if adimx != bdimx || adimy != bdimy {
		return false
	}
	return true
}

// Add a matrix
func (a M) Add(b M) M {
	if len(a) == 1 && len(a[0]) == 1 {
		result := make(M, len(b))
		// a is a scalar
		for i, row := range b {
			for j := range row {
				result[i] = append(result[i], a[0][0]+b[i][j])
			}
		}
		return result
	}
	result := make(M, len(a))
	if len(b) == 1 && len(b[0]) == 1 {
		// b is a scalar
		for i, row := range a {
			for j := range row {
				result[i] = append(result[i], a[i][j]+b[0][0])
			}
		}
		return result
	}
	if !hasSameDimensions(a, b) {
		panic("matrix dimension mismatch, trying to add matrices with different dimensions")
	}
	for i, row := range a {
		for j := range row {
			result[i] = append(result[i], a[i][j]+b[i][j])
		}
	}
	return result
}

func (a M) GetDimensions() (dimx, dimy int) {
	return len(a[0]), len(a)
}

// Subtract a matrix
func (a M) Sub(b M) M {
	if len(a) == 1 && len(a[0]) == 1 {
		result := make(M, len(b))
		// a is a scalar
		for i, row := range b {
			for j := range row {
				result[i] = append(result[i], a[0][0]-b[i][j])
			}
		}
		return result
	}
	result := make(M, len(a))
	if len(b) == 1 && len(b[0]) == 1 {
		// b is a scalar
		for i, row := range a {
			for j := range row {
				result[i] = append(result[i], a[i][j]-b[0][0])
			}
		}
		return result
	}
	if !hasSameDimensions(a, b) {
		panic("dimension mismatch, trying to subtract matrices of different dimensions")
	}
	for i, row := range a {
		for j := range row {
			result[i] = append(result[i], a[i][j]-b[i][j])
		}
	}
	return result
}

// perform a matrix cross product
func (a M) X(b M) M {
	// log.Println(len(b[0]), len(a))
	result := zeroes(len(b[0]), len(a))
	if len(a) == 1 && len(a[0]) == 1 {
		// a is a scalar
		return b.Xscalar(a[0][0])
	}
	if len(b) == 1 && len(b[0]) == 1 {
		// b is a scalar
		return a.Xscalar(b[0][0])
	}
	if len(a[0]) != len(b) {
		log.Println("a column:", len(a[0]), ", b row:", len(b))
		panic("row column mismatch")
	}
	for i := 0; i < len(a); i++ {
		for j := 0; j < len(b[0]); j++ {
			for k := 0; k < len(b); k++ {
				// make sure j doesnt exceed the number of columns in b
				if j < len(b[0]) {
					result[i][j] += a[i][k] * b[k][j]
				}
			}
		}
	}
	// log.Println(result)
	return result
}

// product with a scalar
func (a M) Xscalar(b float64) M {
	result := zeroes(len(a[0]), len(a))
	for i, row := range a {
		for j, col := range row {
			result[i][j] = b * col
		}
	}
	return result
}

func Zeroes(dimx, dimy int) M {
	return zeroes(dimx, dimy)
}

// create a zero matrix from supplied dimensions
func zeroes(dimx, dimy int) M {
	m := make(M, dimy)
	for i := range m {
		for j := 0; j < dimx; j++ {
			m[i] = append(m[i], float64(0))
		}
	}
	return m
}

// check if matrix is a square matrix
func (a M) IsSquareMatrix() bool {
	return len(a) == len(a[0])
}

// Transpose of the matrix
func (a M) T() M {
	dimy := len(a)
	dimx := len(a[0])
	transposeMat := zeroes(dimy, dimx)
	for i, rows := range a {
		for j := range rows {
			transposeMat[j][i] = a[i][j]
		}
	}
	return transposeMat
}

// calculate determinant of 2x2 matrix
func (a M) Det2x2() float64 {
	return a[0][0]*a[1][1] - a[0][1]*a[1][0]
}

// calculate adjoint of 2x2 matrix
func (a M) Adj2x2() M {
	result := zeroes(2, 2)
	result[0][0] = a[1][1]
	result[0][1] = -1 * a[0][1]
	result[1][0] = -1 * a[1][0]
	result[1][1] = a[0][0]
	return result
}

// calculate inverse of 2x2 matrix
func (a M) Inv2x2() M {
	det := a.Det2x2()
	inverse := a.Adj2x2().Xscalar(1 / det)
	return inverse
}

// round off float64 to nearest integer
func (a M) Round() M {
	for i, row := range a {
		for j := range row {
			a[i][j] = math.Round(a[i][j])
		}
	}
	return a
}
