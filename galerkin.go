package main

import (
	//"fmt"
	"math"
	//"gonum.org/v1/gonum/mat"
	"github.com/james-bowman/sparse"
	"gonum.org/v1/gonum/mat"
)

func BinomialCoefficient(n, k int) float64 {
	if k == 0 && n == 0 {
		return 1
	}
	if 0 <= k && k < n {
		return float64(Factorial(n)) / (float64(Factorial(k)) * float64(Factorial(n-k)))
	}
	return 0
}

func Factorial(n int) int {
	if n == 0 {
		return 1
	}
	return n * Factorial(n-1)
}

func LegendreBasis(n, derivativeOrder int, x float64) float64 {
	if n < 0 || derivativeOrder < 0 || derivativeOrder > 2 {
		return 0.0
	}

	scalar := math.Pow(0.5, float64(n))
	pN := 0.0

	for k := 0; k <= n; k++ {
		scalar2 := math.Pow(BinomialCoefficient(n, k), 2)
		temp := 0.0

		switch derivativeOrder {
		case 0:
			temp = legendreBasisRegular(n, k, x)
		case 1:
			temp = legendreBasisDerivative(n, k, x)
		case 2:
			temp = legendreBasisSecondDerivative(n, k, x)
		}

		pN += temp * scalar2
	}

	return scalar * pN
}

func legendreBasisRegular(n, k int, x float64) float64 {
	return math.Pow(x-1.0, float64(n-k)) * math.Pow(x+1.0, float64(k))
}

func legendreBasisDerivative(n, k int, x float64) float64 {
	return BinomialCoefficient(n, k) * BinomialCoefficient(n, k) *
		math.Pow(x-1.0, float64(n-k-1)) * math.Pow(x+1.0, float64(k-1)) *
		(float64(n)*(x+1.0) - 2.0*float64(k))
}

func legendreBasisSecondDerivative(n, k int, x float64) float64 {
	return math.Pow(x-1, float64(n-k-2)) * math.Pow(x+1, float64(k-2)) *
		((math.Pow(float64(n), 2)-float64(n))*math.Pow(x, 2) +
			(2*math.Pow(float64(n), 2)-float64(4*k+2)*float64(n)+4*float64(k))*x +
			math.Pow(float64(n), 2)-float64(4*k+1)*float64(n)+4*math.Pow(float64(k), 2))
}

// StiffnessMatrix calculates the stiffness matrix for a given range of x values, number of nodes, and state vector.
//
// Parameters:
// - xA: the starting x value of the range
// - xB: the ending x value of the range
// - N: the number of nodes
// - u: the state vector
//
// Returns:
// - kCSR: the stiffness matrix as a sparse CSR matrix
func StiffnessMatrix(xA, xB float64, N int, u mat.VecDense) sparse.CSR {
	kDOK := sparse.NewDOK(N, N)
	for i := 0; i < N; i++ {
		for j := 0; j < N; j++ {
			kDOK.Set(i, j, Trapezoid(xA, xB, u))
		}
	}
	kCSR := kDOK.ToCSR()
	return *kCSR
}

// Trapezoid calculates the numerical approximation of an integral using the Trapezoid method.
//
// Parameters:
// - xA: the lower limit of integration.
// - xB: the upper limit of integration.
// - u: the vector of function values at equally spaced points.
//
// Returns:
// - the approximate value of the integral.
func Trapezoid(xA, xB float64, u mat.VecDense) float64 {
	//TODO the length of the state vector?
	uLen := u.Cap()
	h := (xB - xA) / float64(uLen)
	uSum := 0.0
	for i := 1; i < uLen-2; i++ {
		uSum += u.At(i, 0)
	}
	trap := h / 2.0 * (u.At(0, 0) + u.At(uLen - 1, 0) + 2.0*uSum)
	return trap
}

func CentralFluxSolver(pI, pJ float64) float64 {
	return (pI + pJ) / 2.0
}