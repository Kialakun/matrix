package matrix

import "testing"

func TestCholeskyDecomposition(t *testing.T) {
	a := M{{4, 12, -16}, {12, 37, -43}, {-16, -43, 98}}
	t.Log(a.CholeskyDecomposition())
}
