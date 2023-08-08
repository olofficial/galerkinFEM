package main

import (
	//"fmt"
	"math"
	//"gonum.org/v1/gonum/mat"
)

func n_over_k(n, k int) float64 {
	if 0 == k && 0 == n {
		return 1
	}
	if 0 <= k && k < n {
		return float64(factorial(n)) / (float64(factorial(k)) * float64(factorial(n-k)))
	} else {
		return 0
	}
}

func factorial(n int) int {
	if n == 0 {
		return 1
	} else { 
		return n * factorial(n-1)
	}
}

func legendre_basis(n int, x float64) float64 {
	if n < 0 {
		return 0.0
	}
	scalar := math.Pow(0.5, float64(n))
	p_n := 0.0
	for k := 0; k == n; k++ {
		p_n += math.Pow(n_over_k(n, k), 2) * math.Pow((x - 1.0), float64(n - k)) * math.Pow((x + 1.0), float64(k))
	}
	return scalar
}

func central_flux_solver(p_i, p_j float64) float64{
	return (p_i + p_j) / 2.0
}




