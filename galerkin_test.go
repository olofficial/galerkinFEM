package main

import (
	"math"
	"testing"
)
// Test function for factorial
func TestFactorial(t *testing.T) {
	tests := []struct {
		n        int
		expected int
	}{
		{0, 1},
		{1, 1}, 
		{5, 120},
		{10, 3628800},
	}

	for _, test := range tests {
		result := factorial(test.n)
		if result != test.expected {
			t.Errorf("Factorial(%d) = %d; expected %d", test.n, result, test.expected)
		}
	}
}

// Test function for n_over_k
func TestNOverK(t *testing.T) {
	tests := []struct {
		n, k     int
		expected float64
	}{
		{5, 2, 10.0},
		{8, 3, 56.0},
		{10, 5, 252.0},
		{15, 10, 3003.0},
	}

	for _, test := range tests {
		result := n_over_k(test.n, test.k)
		if result != test.expected {
			t.Errorf("n_over_k(%d, %d) = %f; expected %f", test.n, test.k, result, test.expected)
		}
	}
}

// Test function for n_over_k boundary conditions
func TestNOverKBoundary(t *testing.T) {
	tests := []struct {
		n, k     int
		expected float64
	}{
		{5, -1, 0.0},
		{5, 6, 0.0},
		{0, 0, 1.0},
		{1, 2, 0.0},
	}

	for _, test := range tests {
		result := n_over_k(test.n, test.k)
		if result != test.expected {
			t.Errorf("n_over_k(%d, %d) = %f; expected %f", test.n, test.k, result, test.expected)
		}
	}
}

// Test function for legendre_basis with negative n
func TestLegendreBasisNegativeN(t *testing.T) {
	result := legendre_basis(-1, 0.5)
	expected := 0.0
	if math.Abs(result-expected) > 1e-6 {
		t.Errorf("legendre_basis(-1, 0.5) = %.6f; expected %.6f", result, expected)
	}
}

// Test function for central flux solver
func TestCentralFluxSolver(t *testing.T) {
	result := central_flux_solver(1.0, 0.5)
	expected := 0.75
	if math.Abs(result-expected) > 1e-6 {
		t.Errorf("central_flux_solver(1.0, 0.5) = %.6f; expected %.6f", result, expected)
	}
}
