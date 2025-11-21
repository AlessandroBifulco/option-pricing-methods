# Option Pricing Methods

This repository contains two MATLAB implementations for option pricing:

## 1. American Put (Binomial Model)
- Uses a binomial tree with increasing refinement.
- Compares numerical results with the Black–Scholes closed-form formula.
- Includes a convergence plot.

## 2. American Put (Longstaff–Schwartz Method)
- Implements the Least-Squares Monte Carlo algorithm.
- Uses antithetic variates to improve simulation stability.
- Prices American options where early exercise is possible.

## Files
- `europeanPutBinomial.m` — Binomial pricing for European puts.
- `LSM.m` — Longstaff–Schwartz Monte Carlo for American puts.

## Requirements
- MATLAB R2020 or later (or equivalent compatible version).

## Author
Alessandro Bifulco
