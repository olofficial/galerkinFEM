package main

import (
	"math"
	"testing"
)

func TestBinomialCoefficient(t *testing.T) {
	tests := []struct {
		n, k     int
		expected float64
	}{
		{0, 0, 1},
		{5, 2, 10},
		{3, 5, 0},
	}

	for _, test := range tests {
		result := BinomialCoefficient(test.n, test.k)
		if result != test.expected {
			t.Errorf("BinomialCoefficient(%d, %d) = %.2f, expected %.2f",
				test.n, test.k, result, test.expected)
		}
	}
}

func TestFactorial(t *testing.T) {
	tests := []struct {
		n        int
		expected int
	}{
		{0, 1},
		{1, 1},
		{5, 120},
	}

	for _, test := range tests {
		result := Factorial(test.n)
		if result != test.expected {
			t.Errorf("Factorial(%d) = %d, expected %d",
				test.n, result, test.expected)
		}
	}
}

func TestLegendreBasis(t *testing.T) {
	tests := []struct {
		n, derivativeOrder int
		x, expected        float64
	}{
		{0, 0, 0.5, 1.0},
		{3, 0, 1.0, 0.0},
		{3, 1, 1.0, 0.0},
		{5, 2, -1.0, 0.0},
	}

	for _, test := range tests {
		result := LegendreBasis(test.n, test.derivativeOrder, test.x)
		if math.Abs(result-test.expected) > 1e-6 {
			t.Errorf("LegendreBasis(%d, %d, %.2f) = %.6f, expected %.6f",
				test.n, test.derivativeOrder, test.x, result, test.expected)
		}
	}
}

func TestCentralFluxSolver(t *testing.T) {
	pI := 2.0
	pJ := 3.0
	expected := 2.5

	result := CentralFluxSolver(pI, pJ)
	if result != expected {
		t.Errorf("CentralFluxSolver(%.2f, %.2f) = %.2f, expected %.2f",
			pI, pJ, result, expected)
	}
}