# Fixed-Income Derivatives Pricing Engine (C++)

This project implements a production-quality pricing system in C++ for custom fixed-income instruments and derivative contracts with delivery optionality. It includes full support for analytical sensitivity calculations, rate convention modeling, and stochastic forward yield simulation.

## 🔍 Overview

The system handles two financial instruments:

- **Fixed-Income Instrument** (analogous to bonds): Calculates present value using linear, cumulative, and recursive rate conventions. Implements first and second-order analytical derivatives for price–yield sensitivity (duration and convexity).

- **Delivery-Optional Derivative** (akin to cheapest-to-deliver bond options): Selects the most economical instrument from a basket at expiry, using relative normalization and future yield modeling via a geometric Brownian motion.

## 🧠 Key Features

- **Multi-Convention Pricing**: Supports linear, compound (cumulative), and recursive compounding for price–yield transformations.
- **Sensitivity Analysis**: Computes ∂Price/∂Rate, ∂²Price/∂Rate², and their inverses fully analytically—no numerical approximations.
- **Stochastic Forward Modeling**: Uses correlated Brownian motion to simulate forward yields across instruments.
- **Quadratic Approximation**: Approximates price-to-relative-factor ratios using weighted least squares across a standard normal domain.
- **Delivery Probability Calculation**: Integrates minimum ratio quadratics to obtain exact delivery probabilities for each instrument.
- **Greeks-like Derivatives**: Implements analytical sensitivities of derivative contract price w.r.t. volatility and today’s price of each underlying.

## 📐 Mathematical Methods

- **Price–Yield Formulae** for different compounding conventions
- **Closed-form Derivatives**:
  - `dP/dY`, `d²P/dY²`
  - `dY/dP`, `d²Y/dP²`
- **Risk-adjusted Forward Yield Calculation** via quadratic root solving
- **Analytical Integration**: Weighted integration of piecewise quadratics against normal PDF for expected delivery pricing

## 🧱 Code Structure

- `ValueNote` class: Implements fixed-income instrument pricing and sensitivities.
- `DeliveryContract` class: Inherits from `Basket`; handles delivery-option pricing, stochastic modeling, and derivative sensitivities.
- `main()`: Orchestrates user input, contract setup, pricing, and outputs.

## 📊 Example Capabilities

- Price a synthetic bond for a given yield and compute sensitivities
- Model forward prices under stochastic yield evolution
- Identify and price the most cost-effective instrument to deliver at contract expiry
- Compute:
  - ∂Price/∂Volatility of underlying
  - ∂Price/∂Spot price of underlying

## 🧪 Usage

```bash
cd engine
g++ pricing_engine.cpp -o pricing_engine -std=c++17 -O2
./pricing_engine

