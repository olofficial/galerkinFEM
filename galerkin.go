package main

import (
	//"fmt"
	"math"
	//"gonum.org/v1/gonum/mat"
	"github.com/james-bowman/sparse"
	"gonum.org/v1/gonum/mat"
)

// BinomialCoefficient calculates the binomial coefficient of two integers.
//
// Parameters:
// - n: an integer representing the total number of items.
// - k: an integer representing the number of items to be chosen.
//
// Return type: float64
func BinomialCoefficient(n, k int) float64 {
	if k == 0 && n == 0 {
		return 1.0
	}
	if 0 <= k && k < n {
		return float64(Factorial(n)) / (float64(Factorial(k)) * float64(Factorial(n-k)))
	}
	return 0.0
}

// Factorial calculates the factorial of a given integer.
//
// The parameter n is the integer for which the factorial is calculated.
// The function returns an integer value which is the factorial of the input.
func Factorial(n int) int {
	if n == 0 {
		return 1
	}
	return n * Factorial(n-1)
}

// LegendreBasis calculates the Legendre basis function for a given order and derivative order.
//
// Parameters:
// - n: the order of the basis function (must be non-negative)
// - derivativeOrder: the order of the derivative (must be between 0 and 2)
// - x: the input value
//
// Returns:
// - float64: the value of the Legendre basis function
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

// legendreBasisRegular calculates the Legendre basis function for regular polynomials.
//
// Parameters:
//   n - the degree of the polynomial
//   k - the index of the basis function
//   x - the input value
//
// Returns:
//   float64 - the value of the Legendre basis function
func legendreBasisRegular(n, k int, x float64) float64 {
	return math.Pow(x-1.0, float64(n-k)) * math.Pow(x+1.0, float64(k))
}

// legendreBasisDerivative calculates the derivative of the Legendre basis function.
//
// Parameters:
// - n: an integer representing the order of the Legendre basis function.
// - k: an integer representing the index of the Legendre basis function.
// - x: a float64 representing the input value.
//
// Return:
// - a float64 representing the value of the derivative of the Legendre basis function.
func legendreBasisDerivative(n, k int, x float64) float64 {
	return BinomialCoefficient(n, k) * BinomialCoefficient(n, k) *
		math.Pow(x-1.0, float64(n-k-1)) * math.Pow(x+1.0, float64(k-1)) *
		(float64(n)*(x+1.0) - 2.0*float64(k))
}

// legendreBasisSecondDerivative calculates the second derivative of the Legendre basis function.
//
// Parameters:
// - n: an integer representing the degree of the basis function
// - k: an integer representing the order of the basis function
// - x: a float64 representing the input value
//
// Return type: float64
func legendreBasisSecondDerivative(n, k int, x float64) float64 {
	floatN := float64(n)
	floatK := float64(k)

	term1 := math.Pow(x-1, floatN-floatK-2) * math.Pow(x+1, floatK-2)
	term2 := (floatN*floatN - floatN) * x * x
	term3 := (2*floatN*floatN - (4.0*floatK+2.0)*floatN + 4*floatK) * x
	term4 := floatN*floatN - (4.0*floatK+1.0)*floatN + 4*math.Pow(floatK, 2)
	return term1 * (term2 + term3 + term4)
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