package main

import (
	//"fmt"
	"math"
	//"gonum.org/v1/gonum/mat"
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

func CentralFluxSolver(pI, pJ float64) float64 {
	return (pI + pJ) / 2.0
}