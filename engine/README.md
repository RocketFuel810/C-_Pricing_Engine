# ValueNote and DeliveryContract Pricing Model

This file details the working of the application.

## Overview

This C++ implementation provides a comprehensive framework for pricing ValueNotes and DeliveryContracts, including support for multiple rate conventions, volatility modeling, and optimal delivery selection.

## Features

### ValueNote Class
- **Multiple Rate Conventions**: Linear, Cumulative, and Recursive
- **Flexible Payment Frequencies**: Configurable payment schedules
- **Volatility Integration**: Brownian motion-based rate modeling
- **Price-Rate Conversions**: Newton-Raphson method for accurate conversions
- **Forward Price Calculation**: Risk-adjusted forward pricing
- **Analytical Derivatives**: First and second derivatives for sensitivity analysis

### DeliveryContract Class
- **Basket Management**: Support for multiple ValueNotes
- **Quadratic Approximation**: Weighted least squares fitting for optimal delivery
- **Intersection Analysis**: Automatic detection of optimal delivery regions
- **Probability Calculation**: Delivery probabilities for each ValueNote
- **Sensitivity Analysis**: Volatility and price sensitivity calculations

## Mathematical Framework

### Rate Conventions

1. **Linear Convention**: Simple linear discounting
   ```
   Price = Notional × (1 - EffectiveRate × Maturity / 100)
   ```

2. **Cumulative Convention**: Compound discounting with periodic payments
   ```
   Price = Σ(PaymentAmount / (1 + EffectiveRate/(100×Freq))^(Freq×Time)) + PrincipalPV
   ```

3. **Recursive Convention**: Forward-looking recursive calculation
   ```
   FV = (FV + Payment) × (1 + EffectiveRate × TimeDiff / 100)
   Price = (Notional + FV) / (1 + EffectiveRate × Maturity / 100)
   ```

### Risk-Adjusted Pricing

The model incorporates volatility through Brownian motion:
```
EffectiveRate = RiskAdjustedER × exp(σ√T × z - 0.5σ²T)
```

### Quadratic Approximation

For optimal delivery selection, the model fits quadratic functions:
```
Ratio(z) = a×z² + b×z + c
```

## Usage

### Basic Example

```cpp
#include "cla.cpp"

int main() {
    // Create ValueNotes
    auto vn1 = std::make_shared<ValueNote>(100, 5.0, 3.5, 1, 1.5);
    auto vn2 = std::make_shared<ValueNote>(100, 1.5, 2.0, 2, 2.5);
    auto vn3 = std::make_shared<ValueNote>(100, 4.5, 3.25, 1, 1.5);
    auto vn4 = std::make_shared<ValueNote>(100, 10.0, 8.0, 4, 5.0);
    
    std::vector<std::shared_ptr<ValueNote>> basket = {vn1, vn2, vn3, vn4};
    
    // Create DeliveryContract
    DeliveryContract dc(basket, 0.25, 5.0, RelativeFactorMethod::CUMULATIVE, 0.04);
    
    // Calculate price and probabilities
    double price = dc.calculatePrice();
    auto probabilities = dc.calculateDeliveryProbabilities();
    
    return 0;
}
```

### ValueNote Configuration

```cpp
// Parameters: Notional, Maturity, ValueRate, PaymentFreq, Volatility, Convention
ValueNote vn(100.0,    // Notional amount
             5.0,      // Maturity in years
             3.5,      // Value rate (%)
             1,        // Payment frequency (1 = annual)
             1.5,      // Volatility
             RateConvention::CUMULATIVE);
```

### DeliveryContract Configuration

```cpp
// Parameters: Basket, ExpiryTime, StandardizedValueRate, RFMethod, RiskFreeRate
DeliveryContract dc(basket,                    // Vector of ValueNotes
                    0.25,                      // Expiry time in years
                    5.0,                       // Standardized value rate
                    RelativeFactorMethod::CUMULATIVE,  // Relative factor method
                    0.04);                     // Risk-free rate
```

## Compilation

### Prerequisites
- C++11 or later compiler
- Standard math library

### Compile Command
```bash
g++ -std=c++11 -O2 -o cla cla.cpp -lm
```

### Run
```bash
./cla
```

## Key Classes and Methods

### ValueNote Class

#### Core Methods
- `rateToPrice(double effectiveRate)`: Convert rate to price
- `priceToRate(double price)`: Convert price to rate (Newton-Raphson)
- `forwardPrice(double currentPrice, double riskFreeRate, double timeToExpiry)`: Calculate forward price
- `calculateRiskAdjustedRate(double forwardER, double timeToExpiry)`: Risk-adjusted rate calculation

#### Sensitivity Methods
- `priceDerivative(double effectiveRate)`: First derivative
- `priceSecondDerivative(double effectiveRate)`: Second derivative

### DeliveryContract Class

#### Core Methods
- `calculatePrice()`: Calculate contract price
- `calculateDeliveryProbabilities()`: Calculate delivery probabilities
- `calculateVolatilitySensitivity()`: Volatility sensitivity analysis
- `findOptimalValueNote(double z)`: Find optimal ValueNote for given z

#### Utility Methods
- `buildQuadraticApproximations()`: Build quadratic fits
- `findIntersectionPoints()`: Find intersection points between quadratics
- `calculateRelativeFactor(const ValueNote& vn)`: Calculate relative factors

## Mathematical Utilities

The `MathUtils` namespace provides:
- `normalPDF(double x)`: Standard normal probability density function
- `normalCDF(double x)`: Standard normal cumulative distribution function
- `solveQuadratic(double a, double b, double c)`: Solve quadratic equations

## Output Format

The program outputs:
1. **DeliveryContract Price**: The calculated contract value
2. **Delivery Probabilities**: Probability of each ValueNote being delivered
3. **Relative Factors**: Relative factor for each ValueNote

Example output:
```
DeliveryContract Price: 95.123456
Delivery Probabilities:
VN1: 0.234567
VN2: 0.345678
VN3: 0.234567
VN4: 0.185188
Relative Factors:
VN1: 0.950000
VN2: 0.980000
VN3: 0.967500
VN4: 0.800000
```

## Error Handling

The implementation includes comprehensive error handling:
- Invalid quadratic equations
- Numerical convergence issues
- Invalid parameter ranges
- Memory allocation failures

## Performance Considerations

- **Quadratic Fitting**: Uses 2000 sample points for accurate approximation
- **Integration**: Analytical integration for normal-weighted quadratics
- **Optimization**: Efficient intersection point detection and sorting
- **Memory Management**: Smart pointers for automatic memory management

## Extensions and Customization

### Adding New Rate Conventions
1. Add new enum value to `RateConvention`
2. Implement pricing logic in `rateToPrice()`
3. Add derivative calculations in `priceDerivative()`

### Custom Relative Factor Methods
1. Add new enum value to `RelativeFactorMethod`
2. Implement calculation logic in `calculateRelativeFactor()`

### Volatility Models
The current implementation uses simple Brownian motion. Extensions could include:
- Stochastic volatility models
- Mean reversion
- Jump diffusion processes

## Dependencies

- **Standard C++ Library**: `<iostream>`, `<vector>`, `<cmath>`, `<algorithm>`
- **Memory Management**: `<memory>` for smart pointers
- **Error Handling**: `<stdexcept>` for exceptions
- **I/O Formatting**: `<iomanip>` for output formatting

## License

This implementation is provided for educational and research purposes.

## Contributing

To extend or improve this implementation:
1. Add comprehensive unit tests
2. Implement additional volatility models
3. Add Monte Carlo simulation capabilities
4. Optimize numerical methods for better performance
5. Add support for more complex payoff structures 
