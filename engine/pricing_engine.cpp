#include <iostream>
#include <cmath>
#include <vector>
#include <algorithm>

#define M_PI 3.14159265358979323846
#define M_SQRT1_2 0.70710678118654752440

const int NUM_SAMPLES = 2000;
const double Z_MIN = -3.0;
const double Z_MAX = 3.0;

struct QuadraticCoeffs {
    double a, b, c;  // coefficients for ax^2 + bx + c
    int vnIndex;     // index of the ValueNote this quadratic represents
};

class ValueNote {
public:
    double notional;
    double maturity; // in years
    double valueRate; // [0,100] 
    int paymentFreq; // >=1
    double effectiveRateVolatility; // [0,100]

    // Constructor
    ValueNote(double N, double M, double VR, int PF, double ERV)
        : notional(N), maturity(M), valueRate(VR), paymentFreq(PF) {
        if (N <= 0 || M <= 0 || VR <= 0 || PF <= 0) {
            throw std::invalid_argument("Invalid input parameters for ValueNote");
        }
    }

    // Q1.1

    // Linear a) Price given effective rate
    double eRateToPriceLinear(double eRate) const {
        return notional * (1 - eRate * maturity / 100.0);
    }

    // Linear b) Effective rate given price
    double priceToERateLinear(double price) const {
        return (1 - price / notional) * (100.0 / maturity);
    }

    // Cumulative a) Price given effective rate
    double eRateToPriceCumulative(double eRate) const {
        int n = maturity * paymentFreq;
        double price = 0.0;
        double interest = notional * valueRate / (100.0 * paymentFreq);
        for (int i = 1; i < n; i++) {
            price += interest / pow(1 + eRate / (paymentFreq * 100.0), i);
        }
        price += (interest + notional) / pow(1 + eRate / (paymentFreq * 100.0), n);
        return price;
    }

    // Cumulative b) Effective rate given price
    double priceToERateCumulative(double price, double tolerance = 1e-8, int maxIter = 1000) const {
        double low = 0.0;
        double high = 100.0;
        int i = 0;
        while (high - low > tolerance && i < maxIter) {
            double mid = (low + high) / 2.0;
            double diff = eRateToPriceCumulative(mid) - price;
            if (diff < 0) high = mid;
            else low = mid;
            i++;
        }
        return (low + high) / 2.0;
    }

    // Recursive a) Price given effective rate
    double eRateToPriceRecursive(double eRate) const {
        int n = maturity * paymentFreq;
        double interest = notional * valueRate / (100.0 * paymentFreq);
        double FV = 0;
        double t = 1 / static_cast<double>(paymentFreq);
        for (int i = 1; i < n; i++) {
            double t_next = static_cast<double>(std::min(i+1, n)) / static_cast<double>(paymentFreq);
            double m = t_next - t;
            FV = (FV + interest) / (1 + eRate * m / 100.0);
            t = t_next;
        }
        FV += interest;
        return (FV + notional) / (1 + eRate * maturity / 100.0);
    }

    // Recursive b) Effective rate given price
    double priceToERateRecursive(double price, double tolerance = 1e-8, int maxIter = 1000) const {
        double low = 0.0;
        double high = 100.0;
        int i = 0;
        while (high - low > tolerance && i < maxIter) {
            double mid = (low + high) / 2.0;
            double diff = eRateToPriceRecursive(mid) - price;
            if (diff < 0) high = mid;
            else low = mid;
            i++;
        }
        return (low + high) / 2.0;
    }

    // Q1.2

    // Linear a) Price sensitivity given effective rate
    double eRateToPriceSensitivityLinear(double eRate) const {
        return -notional * maturity / 100.0;
    }

    // Linear b) Effective rate sensitivity given price
    double priceToERateSensitivityLinear(double price) const {
        return -100.0 / (notional * maturity);
    }

    // Cumulative a) Price sensitivity given effective rate
    double eRateToPriceSensitivityCumulative(double eRate) const {
        int n = maturity * paymentFreq;
        double price = 0.0;
        double interest = notional * valueRate / (100.0 * paymentFreq);
        for (int i = 1; i < n; i++) {
            price += (i * interest) / pow(1 + eRate / (paymentFreq * 100.0), i + 1);
        }
        price += (n * (interest + notional)) / pow(1 + eRate / (paymentFreq * 100.0), n + 1);
        price = -price / (100.0 * paymentFreq);
        return price;
    }

    // Cumulative b) Effective rate sensitivity given price
    double priceToERateSensitivityCumulative(double price) const {
        double eRate = priceToERateCumulative(price);
        double priceSensitivity = eRateToPriceSensitivityCumulative(eRate);
        return 1.0 / priceSensitivity; 
    }

    // Recursive a) Price sensitivity given effective rate
    double eRateToPriceSensitivityRecursive(double eRate) const {
        int n = maturity * paymentFreq;
        double interest = notional * valueRate / (100.0 * paymentFreq);
        double FV = 0.0, dFV = 0.0;
        double t = 1 / static_cast<double>(paymentFreq);
        for (int i = 1; i < n; ++i) {
            double t_next = static_cast<double>(std::min(i+1, n)) / static_cast<double>(paymentFreq);
            double m = t_next - t;
            dFV = dFV / (1 + eRate * m / 100.0) + (FV + interest) * (-m / 100.0) / pow(1 + eRate * m / 100.0, 2);
            FV = (FV + interest) / (1 + eRate * m / 100.0);
            t = t_next;
        }
        dFV += -interest * maturity / 100.0 / pow(1 + eRate * maturity / 100.0, 2);
        FV += interest;
        return dFV / (1 + eRate * maturity / 100.0) - (FV + notional) * maturity / 100.0 / pow(1 + eRate * maturity / 100.0, 2);
    }

    // Recursive b) Effective rate sensitivity given price
    double priceToERateSensitivityRecursive(double price) const {
        double eRate = priceToERateRecursive(price);
        double priceSensitivity = eRateToPriceSensitivityRecursive(eRate);
        return 1.0 / priceSensitivity;
    }   

    // Q1.3

    // Linear a) Price second derivative wrt effective rate given effective rate
    double eRateToPriceSecondDerivativeLinear(double eRate) const {
        return 0.0;
    }

    // Linear b) Effective rate second derivative wrt price given price
    double priceToERateSecondDerivativeLinear(double price) const {
        return 0.0;
    }

    // Cumulative a) Price second derivative wrt effective rate given effective rate
    double eRateToPriceSecondDerivativeCumulative(double eRate) const {
        int n = maturity * paymentFreq;
        double d2P = 0.0;
        double interest = notional * valueRate / (100.0 * paymentFreq);
        for (int i = 1; i < n; ++i) {
            double denom = pow(1 + eRate / (paymentFreq * 100.0), i + 2);
            d2P += i * (i + 1) * interest / pow(100.0 * paymentFreq, 2) / denom;
        }
        double finalFlow = interest + notional;
        d2P += n * (n + 1) * finalFlow / pow(100.0 * paymentFreq, 2) / pow(1 + eRate / (paymentFreq * 100.0), n + 2);
        return d2P;
    }

    // Cumulative b) Effective rate second derivative wrt price given price
    double priceToERateSecondDerivativeCumulative(double price) const {
        double eRate = priceToERateCumulative(price);
        double f1 = eRateToPriceSensitivityCumulative(eRate);
        double f2 = eRateToPriceSecondDerivativeCumulative(eRate);
        return -f2 / (f1 * f1 * f1);
    }

    // Recursive a) Price second derivative wrt effective rate given effective rate
    double eRateToPriceSecondDerivativeRecursive(double eRate) const {
        int n = maturity * paymentFreq;
        double interest = notional * valueRate / (100.0 * paymentFreq);
        double FV = 0.0, dFV = 0.0, d2FV = 0.0;
        double t = 1 / static_cast<double>(paymentFreq);
        for (int i = 1; i < n; ++i) {
            double t_next = static_cast<double>(std::min(i+1, n)) / static_cast<double>(paymentFreq);
            double m = t_next - t;
            d2FV = d2FV / (1 + eRate * m / 100.0) +
                   2 * dFV * (-m / 100.0) / pow(1 + eRate * m / 100.0, 2) +
                   (FV + interest) * (2 * pow(m / 100.0, 2)) / pow(1 + eRate * m / 100.0, 3);
            dFV = dFV / (1 + eRate * m / 100.0) + (FV + interest) * (-m / 100.0) / pow(1 + eRate * m / 100.0, 2);
            FV = (FV + interest) / (1 + eRate * m / 100.0);
            t = t_next;
        }
        double D = 1 + eRate * maturity / 100.0;
        double dD = maturity / 100.0;
        return (d2FV / D) - 2 * dFV * dD / (D * D) + 2 * (FV + notional) * dD * dD / (D * D * D);
    }

    // Recursive b) Effective rate second derivative wrt price given price
    double priceToERateSecondDerivativeRecursive(double price) const {
        double eRate = priceToERateRecursive(price);
        double f1 = eRateToPriceSensitivityRecursive(eRate);
        double f2 = eRateToPriceSecondDerivativeRecursive(eRate);
        return -f2 / (f1 * f1 * f1);
    }

    // Q2
    double forwardPrice(double expirationTimeYears, double riskFreeRate) const {
        double VP0 = eRateToPriceCumulative(valueRate); // Price today using its own valueRate as ER
        double compounded = VP0 * (1 + (riskFreeRate / 100.0) * expirationTimeYears);

        double interest = valueRate * notional / (100.0 * paymentFreq);
        double accrued = 0.0;

        for (int i = 1; i <= maturity * paymentFreq; ++i) {
            double t_i = static_cast<double>(i) / paymentFreq;
            if (t_i <= expirationTimeYears) {
                accrued += interest / (1 + (riskFreeRate / 100.0) * t_i);
            } else {
                break;
            }
        }
        std::cout << "[DEBUG] compounded=" << compounded << " accrued=" << accrued << std::endl;
        return compounded - accrued * (1 + (riskFreeRate / 100.0) * expirationTimeYears);
    }

    double getRiskAdjustedEffectiveRateCumulative(double forwardPrice, double expirationTimeYears) const {
        double ER_T = priceToERateCumulative(forwardPrice);
        double dP = eRateToPriceSensitivityCumulative(ER_T);
        double d2P = eRateToPriceSecondDerivativeCumulative(ER_T);
        double sigma = effectiveRateVolatility / 100.0;
        double A = 0.5 * d2P * std::exp(sigma * sigma * expirationTimeYears);
        double B = dP - d2P * ER_T;
        double C = 0.5 * d2P * ER_T * ER_T - dP * ER_T;

        double discriminant = B * B - 4 * A * C;
        if (discriminant < 0) {
            std::cerr << "No real solution for risk-adjusted ER. Returning unadjusted ER.\n";
            return ER_T;
        }

        double root1 = (-B + sqrt(discriminant)) / (2 * A);
        double root2 = (-B - sqrt(discriminant)) / (2 * A);

        // Choose the root closer to ER_T
        return (fabs(root1 - ER_T) < fabs(root2 - ER_T)) ? root1 : root2;
    }

    double getRiskAdjustedEffectiveRate(double forwardPrice, double expirationTimeYears, std::string effectiveRateConvention) const {
        if (effectiveRateConvention == "Cumulative") {
            return getRiskAdjustedEffectiveRateCumulative(forwardPrice, expirationTimeYears);
        } else {
            throw std::invalid_argument("Unsupported effective rate convention.");
        }
    }

};

class Basket {
protected:
    std::vector<ValueNote> basket;

public:
    double standardizedValueRate;
    int expirationTime;
    std::string relativeFactorConvention;
    std::string effectiveRateConvention;
    double riskFreeInterestRate;
    int numberOfNotes;

    Basket(double svr, int expTime, const std::string& rfc, const std::string& erc,
           double rfr, int nNotes)
        : standardizedValueRate(svr), expirationTime(expTime),
          relativeFactorConvention(rfc), effectiveRateConvention(erc),
          riskFreeInterestRate(rfr), numberOfNotes(nNotes) {}

    void createBasket() {
        for (int i = 0; i < numberOfNotes; ++i) {
            double N, M, VR, ERV;
            int PF;
            std::cout << "Enter Notional, Maturity (years), Value Rate (%), Payment Frequency and Effective Rate Volatility for ValueNote " << (i + 1) << ": ";
            std::cin >> N >> M >> VR >> PF >> ERV;
            basket.emplace_back(N, M, VR, PF, ERV);
        }
    }

    const std::vector<ValueNote>& getBasket() const {
        return basket;
    }
};

class DeliveryContract : public Basket {
private:
    std::vector<double> relativeFactors;
    std::vector<QuadraticCoeffs> quadraticApprox;
    std::vector<double> forwardPrices; // Store forward prices for each ValueNote
    std::vector<double> adjustedEffectiveRates; // Store adjusted effective rates for each ValueNote

    // Fit quadratic using weighted least squares
    QuadraticCoeffs fitQuadratic(const std::vector<double>& z, 
                                const std::vector<double>& ratios,
                                int vnIndex) const {
        int n = z.size();
        double s0 = 0, s1 = 0, s2 = 0, s3 = 0, s4 = 0;
        double sy = 0, sxy = 0, sx2y = 0;
        
        // Calculate weighted sums
        for (int i = 0; i < n; ++i) {
            double w = std::exp(-0.5 * z[i] * z[i]); // normal weight
            double xi = z[i];
            double yi = ratios[i];
            
            s0 += w;
            s1 += w * xi;
            s2 += w * xi * xi;
            s3 += w * xi * xi * xi;
            s4 += w * xi * xi * xi * xi;
            
            sy += w * yi;
            sxy += w * xi * yi;
            sx2y += w * xi * xi * yi;
        }
        
        // Solve 3x3 system using Cramer's rule
        double det = s0 * (s2 * s4 - s3 * s3) - s1 * (s1 * s4 - s2 * s3) + s2 * (s1 * s3 - s2 * s2);
        
        if (std::abs(det) < 1e-12) {
            return {0.0, 0.0, sy / s0, vnIndex}; // Fall back to constant
        }
        
        double c = (sy * (s2 * s4 - s3 * s3) - sxy * (s1 * s4 - s2 * s3) + sx2y * (s1 * s3 - s2 * s2)) / det;
        double b = (s0 * (sxy * s4 - sx2y * s3) - sy * (s1 * s4 - s2 * s3) + sx2y * (s1 * s2 - s0 * s3)) / det;
        double a = (s0 * (s2 * sx2y - s3 * sxy) - s1 * (s1 * sx2y - s2 * sxy) + sy * (s1 * s3 - s2 * s2)) / det;
        
        return {a, b, c, vnIndex};
    }

    // Find intersection points between quadratics
    std::vector<double> findIntersections(const std::vector<QuadraticCoeffs>& quadraticApprox) const {
        std::vector<double> intersections;
        for (size_t i = 0; i < quadraticApprox.size(); ++i) {
            for (size_t j = i + 1; j < quadraticApprox.size(); ++j) {
                const auto& q1 = quadraticApprox[i];
                const auto& q2 = quadraticApprox[j];
                double da = q1.a - q2.a;
                double db = q1.b - q2.b;
                double dc = q1.c - q2.c;
                if (std::abs(da) > 1e-12) {
                    double disc = db * db - 4 * da * dc;
                    if (disc >= 0) {
                        double sqrt_disc = std::sqrt(disc);
                        double z1 = (-db + sqrt_disc) / (2 * da);
                        double z2 = (-db - sqrt_disc) / (2 * da);
                        if (z1 >= Z_MIN && z1 <= Z_MAX) intersections.push_back(z1);
                        if (z2 >= Z_MIN && z2 <= Z_MAX) intersections.push_back(z2);
                    }
                } else if (std::abs(db) > 1e-12) {
                    double z = -dc / db;
                    if (z >= Z_MIN && z <= Z_MAX) intersections.push_back(z);
                }
            }
        }
        std::sort(intersections.begin(), intersections.end());
        intersections.erase(std::unique(intersections.begin(), intersections.end()), intersections.end());
        return intersections;
    }

    // Evaluate quadratic function
    double evaluateQuadratic(const QuadraticCoeffs& q, double z) const {
        return q.a * z * z + q.b * z + q.c;
    }

    // Find optimal ValueNote for given z
    int findOptimalValueNote(double z, const std::vector<QuadraticCoeffs>& quadraticApprox) const {
        int optimalIndex = 0;
        double minRatio = quadraticApprox[0].a * z * z + quadraticApprox[0].b * z + quadraticApprox[0].c;
        for (size_t i = 1; i < quadraticApprox.size(); ++i) {
            double ratio = quadraticApprox[i].a * z * z + quadraticApprox[i].b * z + quadraticApprox[i].c;
            if (ratio < minRatio) {
                minRatio = ratio;
                optimalIndex = i;
            }
        }
        return optimalIndex;
    }

    // Integrate quadratic with normal weight
    double integrateQuadraticWithNormalWeight(double a, double b, double c, double z1, double z2) const {
        auto normalCDF = [](double x) { return 0.5 * std::erfc(-x * M_SQRT1_2); };
        auto normalPDF = [](double x) { return std::exp(-0.5 * x * x) / std::sqrt(2.0 * M_PI); };
        
        // For constant term c
        double integralC = c * (normalCDF(z2) - normalCDF(z1));
        
        // For linear term bz
        double integralB = -b * (normalPDF(z2) - normalPDF(z1));
        
        // For quadratic term az^2
        double integralA = a * ((normalCDF(z2) - normalCDF(z1)) - z2 * normalPDF(z2) + z1 * normalPDF(z1));
        
        return integralC + integralB + integralA;
    }

    // Weighted quadratic fit for a vector of y values at zSamples
    QuadraticCoeffs fitQuadraticForDeriv(const std::vector<double>& zSamples, const std::vector<double>& y, int vnIndex) const {
        int n = zSamples.size();
        double s0 = 0, s1 = 0, s2 = 0, s3 = 0, s4 = 0;
        double sy = 0, sxy = 0, sx2y = 0;
        for (int i = 0; i < n; ++i) {
            double w = std::exp(-0.5 * zSamples[i] * zSamples[i]);
            double xi = zSamples[i];
            double yi = y[i];
            s0 += w;
            s1 += w * xi;
            s2 += w * xi * xi;
            s3 += w * xi * xi * xi;
            s4 += w * xi * xi * xi * xi;
            sy += w * yi;
            sxy += w * xi * yi;
            sx2y += w * xi * xi * yi;
        }
        double det = s0 * (s2 * s4 - s3 * s3) - s1 * (s1 * s4 - s2 * s3) + s2 * (s1 * s3 - s2 * s2);
        if (std::abs(det) < 1e-12) {
            return {0.0, 0.0, sy / s0, vnIndex};
        }
        double c = (sy * (s2 * s4 - s3 * s3) - sxy * (s1 * s4 - s2 * s3) + sx2y * (s1 * s3 - s2 * s2)) / det;
        double b = (s0 * (sxy * s4 - sx2y * s3) - sy * (s1 * s4 - s2 * s3) + sx2y * (s1 * s2 - s0 * s3)) / det;
        double a = (s0 * (s2 * sx2y - s3 * sxy) - s1 * (s1 * sx2y - s2 * sxy) + sy * (s1 * s3 - s2 * s2)) / det;
        return {a, b, c, vnIndex};
    }

    // Analytical derivative of risk-adjusted effective rate with respect to VP0 for ValueNote i
    double dAdjustedER_dVP0(int i) const {
        // The quadratic equation is: d2P * ER^2 + dP * ER + (VP_T - dP * ER) = 0
        // We use the Cumulative convention only
        double expirationTimeYears = static_cast<double>(expirationTime) / 12.0;
        double VP0 = forwardPrices[i] / (1 + riskFreeInterestRate * expirationTimeYears); // invert forward price formula
        double dFwd_dVP0 = 1.0 + riskFreeInterestRate * expirationTimeYears;
        double fwdPrice = forwardPrices[i];
        double adjER = adjustedEffectiveRates[i];
        // The quadratic equation is: d2P * ER^2 + dP * ER + (fwdPrice - dP * ER) = 0
        // So Q(ER, fwdPrice) = d2P * ER^2 + dP * ER + (fwdPrice - dP * ER)
        // Partial derivatives:
        // dQ/dER = 2*d2P*ER + dP - dP = 2*d2P*ER
        // dQ/dFwd = 1
        // dFwd/dVP0 = dFwd_dVP0
        // So dER/dVP0 = - (dQ/dFwd * dFwd/dVP0) / (dQ/dER)
        double dP = basket[i].eRateToPriceSensitivityCumulative(adjER);
        double d2P = basket[i].eRateToPriceSecondDerivativeCumulative(adjER);
        double dQ_dER = 2 * d2P * adjER;
        std::cout << "[DEBUG] VN" << (i+1) << " dP=" << dP << " d2P=" << d2P << " adjER=" << adjER << " dQ_dER=" << dQ_dER << std::endl;
        double dQ_dFwd = 1.0;
        double dER_dVP0 = - (dQ_dFwd * dFwd_dVP0) / dQ_dER;
        return dER_dVP0;
    }

    // Generic integration for DeliveryContract price sensitivity wrt a parameter (e.g., sigma or VP0)
    double integrateDeliveryContractSensitivity(
        const std::vector<QuadraticCoeffs>& derivQuads,
        const std::vector<QuadraticCoeffs>& priceQuads
    ) const {
        std::vector<double> intersections = findIntersections(priceQuads);
        intersections.insert(intersections.begin(), Z_MIN);
        intersections.push_back(Z_MAX);
        double totalSensitivity = 0.0;
        for (size_t m = 0; m < intersections.size() - 1; ++m) {
            double z1 = intersections[m];
            double z2 = intersections[m + 1];
            double zMid = (z1 + z2) / 2.0;
            int optimalVN = findOptimalValueNote(zMid, priceQuads);
            const auto& dq = derivQuads[optimalVN];
            totalSensitivity += integrateQuadraticWithNormalWeight(dq.a, dq.b, dq.c, z1, z2);
        }
        return totalSensitivity;
    }

    // Volatility sensitivity for ValueNote i at given sigma_i (in percent), using piecewise integration
    double volatilitySensitivityForNote(int i, double sigma_i) const {
        double expirationTimeYears = static_cast<double>(expirationTime) / 12.0;
        // Prepare z samples
        std::vector<double> zSamples(NUM_SAMPLES);
        double step = (Z_MAX - Z_MIN) / (NUM_SAMPLES - 1);
        for (int j = 0; j < NUM_SAMPLES; ++j) {
            zSamples[j] = Z_MIN + j * step;
        }
        // Compute derivative of ratio wrt sigma_i at each z for all ValueNotes
        std::vector<std::vector<double>> derivsAll(basket.size(), std::vector<double>(NUM_SAMPLES));
        for (size_t k = 0; k < basket.size(); ++k) {
            double sigma = (k == i) ? (sigma_i / 100.0) : (basket[k].effectiveRateVolatility / 100.0);
            double ER_T = adjustedEffectiveRates[k];
            for (int j = 0; j < NUM_SAMPLES; ++j) {
                double z = zSamples[j];
                double ER_Tz = ER_T * std::exp(sigma * std::sqrt(expirationTimeYears) * z - 0.5 * sigma * sigma * expirationTimeYears);
                double dER_dsigma = (k == i) ? (ER_Tz * (std::sqrt(expirationTimeYears) * z - sigma * expirationTimeYears)) : 0.0;
                double dVP_dsigma = basket[k].eRateToPriceSensitivityCumulative(ER_Tz) * dER_dsigma;
                derivsAll[k][j] = dVP_dsigma / relativeFactors[k];
            }
        }
        // Fit quadratics to the derivatives for all ValueNotes
        std::vector<QuadraticCoeffs> derivQuads;
        for (size_t k = 0; k < basket.size(); ++k) {
            derivQuads.push_back(fitQuadraticForDeriv(zSamples, derivsAll[k], k));
        }
        // Fit price quadratics for all ValueNotes
        std::vector<std::vector<double>> ratios(basket.size(), std::vector<double>(NUM_SAMPLES));
        for (int j = 0; j < NUM_SAMPLES; ++j) {
            double z = zSamples[j];
            for (size_t k = 0; k < basket.size(); ++k) {
                double sigma_k = basket[k].effectiveRateVolatility / 100.0;
                double ER_T_k = adjustedEffectiveRates[k];
                double ER_Tz_k = ER_T_k * std::exp(sigma_k * std::sqrt(expirationTimeYears) * z - 0.5 * sigma_k * sigma_k * expirationTimeYears);
                double price = basket[k].eRateToPriceCumulative(ER_Tz_k);
                ratios[k][j] = price / relativeFactors[k];
            }
        }
        std::vector<QuadraticCoeffs> priceQuads;
        for (size_t k = 0; k < basket.size(); ++k) {
            priceQuads.push_back(fitQuadratic(zSamples, ratios[k], k));
        }
        // Integrate using the generic function
        return integrateDeliveryContractSensitivity(derivQuads, priceQuads);
    }

    // Analytical price sensitivity for ValueNote i at given VP0_i, using piecewise integration
    double priceSensitivityForNote(int i, double VP0_i) const {
        double expirationTimeYears = static_cast<double>(expirationTime) / 12.0;
        // Prepare z samples
        std::vector<double> zSamples(NUM_SAMPLES);
        double step = (Z_MAX - Z_MIN) / (NUM_SAMPLES - 1);
        for (int j = 0; j < NUM_SAMPLES; ++j) {
            zSamples[j] = Z_MIN + j * step;
        }
        // Compute derivative of ratio wrt VP0_i at each z for all ValueNotes
        std::vector<std::vector<double>> derivsAll(basket.size(), std::vector<double>(NUM_SAMPLES));
        for (size_t k = 0; k < basket.size(); ++k) {
            double adjER = adjustedEffectiveRates[k];
            double dER_dVP0 = (k == i) ? dAdjustedER_dVP0(i) : 0.0;
            double sigma = basket[k].effectiveRateVolatility / 100.0;
            for (int j = 0; j < NUM_SAMPLES; ++j) {
                double z = zSamples[j];
                double expTerm = std::exp(sigma * std::sqrt(expirationTimeYears) * z - 0.5 * sigma * sigma * expirationTimeYears);
                double ER_Tz = adjER * expTerm;
                double dER_Tz_dVP0 = dER_dVP0 * expTerm;
                double dVP_dVP0 = basket[k].eRateToPriceSensitivityCumulative(ER_Tz) * dER_Tz_dVP0;
                derivsAll[k][j] = dVP_dVP0 / relativeFactors[k];
            }
        }
        // Fit quadratics to the derivatives for all ValueNotes
        std::vector<QuadraticCoeffs> derivQuads;
        for (size_t k = 0; k < basket.size(); ++k) {
            derivQuads.push_back(fitQuadraticForDeriv(zSamples, derivsAll[k], k));
        }
        // Fit price quadratics for all ValueNotes
        std::vector<std::vector<double>> ratios(basket.size(), std::vector<double>(NUM_SAMPLES));
        for (int j = 0; j < NUM_SAMPLES; ++j) {
            double z = zSamples[j];
            for (size_t k = 0; k < basket.size(); ++k) {
                double sigma_k = basket[k].effectiveRateVolatility / 100.0;
                double ER_T_k = adjustedEffectiveRates[k];
                double ER_Tz_k = ER_T_k * std::exp(sigma_k * std::sqrt(expirationTimeYears) * z - 0.5 * sigma_k * sigma_k * expirationTimeYears);
                double price = basket[k].eRateToPriceCumulative(ER_Tz_k);
                ratios[k][j] = price / relativeFactors[k];
            }
        }
        std::vector<QuadraticCoeffs> priceQuads;
        for (size_t k = 0; k < basket.size(); ++k) {
            priceQuads.push_back(fitQuadratic(zSamples, ratios[k], k));
        }
        // Integrate using the generic function
        return integrateDeliveryContractSensitivity(derivQuads, priceQuads);
    }

public:
    DeliveryContract(double svr, int expTime, const std::string& rfc, const std::string& erc,
                     double rfr, int nNotes)
        : Basket(svr, expTime, rfc, erc, rfr, nNotes) {}

    void computeRelativeFactors() {
        relativeFactors.clear();
        if (relativeFactorConvention == "Cumulative") {
            for (const auto& note : basket) {
                double price = note.eRateToPriceCumulative(standardizedValueRate);
                relativeFactors.push_back(price / 100.0);
            }
        } else if (relativeFactorConvention == "Unity") {
            relativeFactors = std::vector<double>(basket.size(), 1.0);
        } else {
            throw std::invalid_argument("Unsupported relative factor convention.");
        }
    }

    void printRelativeFactors() const {
        std::cout << "Relative Factors:\n";
        for (size_t i = 0; i < relativeFactors.size(); ++i) {
            std::cout << "Note " << (i + 1) << ": " << relativeFactors[i] << "\n";
        }
    }

    void printForwardPrices() const {
        std::cout << "\nForward Prices at expiration (t = " << expirationTime << " months):\n";
        for (size_t i = 0; i < basket.size(); ++i) {
            double fwdPrice = basket[i].forwardPrice(static_cast<double>(expirationTime)/12.0, riskFreeInterestRate);
            std::cout << "Note " << (i + 1) << ": " << fwdPrice << "\n";
        }
    }

    void printRiskAdjustedEffectiveRates() const {
        std::cout << "\nRisk-Adjusted Effective Rates at expiration (t = " << expirationTime << " months):\n";
        for (size_t i = 0; i < basket.size(); ++i) {
            double fwdPrice = basket[i].forwardPrice(static_cast<double>(expirationTime)/12.0, riskFreeInterestRate);
            double adjER = basket[i].getRiskAdjustedEffectiveRate(fwdPrice, static_cast<double>(expirationTime)/12.0, effectiveRateConvention);
            std::cout << "Note " << (i + 1) << ": " << adjER << "\n";
        }
    }

    // Call after basket and relative factors are set up
    void computeForwardPricesAndAdjustedRates() {
        forwardPrices.clear();
        adjustedEffectiveRates.clear();
        double expirationTimeYears = static_cast<double>(expirationTime) / 12.0;
        for (size_t i = 0; i < basket.size(); ++i) {
            double fwdPrice = basket[i].forwardPrice(expirationTimeYears, riskFreeInterestRate);
            forwardPrices.push_back(fwdPrice);
            double adjER = basket[i].getRiskAdjustedEffectiveRate(fwdPrice, expirationTimeYears, effectiveRateConvention);
            adjustedEffectiveRates.push_back(adjER);
        }
    }

    // Compute and print both DeliveryContract price and probabilities together
    void computeAndPrintPriceAndProbabilities() {
        double expirationTimeYears = static_cast<double>(expirationTime) / 12.0;
        // Generate sample points and calculate ratios
        std::vector<std::vector<double>> ratios(basket.size(), std::vector<double>(NUM_SAMPLES));
        std::vector<double> zSamples(NUM_SAMPLES);
        double step = (Z_MAX - Z_MIN) / (NUM_SAMPLES - 1);
        for (int j = 0; j < NUM_SAMPLES; ++j) {
            double z = Z_MIN + j * step;
            zSamples[j] = z;
            for (size_t i = 0; i < basket.size(); ++i) {
                double sigma = basket[i].effectiveRateVolatility / 100.0;
                double ER_T = adjustedEffectiveRates[i];
                double ER_Tz = ER_T * std::exp(sigma * std::sqrt(expirationTimeYears) * z -
                                             0.5 * sigma * sigma * expirationTimeYears);
                double price = basket[i].eRateToPriceCumulative(ER_Tz);
                ratios[i][j] = price / relativeFactors[i];
            }
        }
        // Fit quadratics to the ratio data
        std::vector<QuadraticCoeffs> quadraticApproxLocal;
        for (size_t i = 0; i < basket.size(); ++i) {
            quadraticApproxLocal.push_back(fitQuadratic(zSamples, ratios[i], i));
        }
        // Find intersection points
        std::vector<double> intersections = findIntersections(quadraticApproxLocal);
        intersections.insert(intersections.begin(), Z_MIN);
        intersections.push_back(Z_MAX);
        // Calculate price and probabilities by integrating over intervals
        double price = 0.0;
        std::vector<double> probabilities(basket.size(), 0.0);
        for (size_t i = 0; i < intersections.size() - 1; ++i) {
            double z1 = intersections[i];
            double z2 = intersections[i + 1];
            double zMid = (z1 + z2) / 2.0;
            int optimalVN = findOptimalValueNote(zMid, quadraticApproxLocal);
            const auto& q = quadraticApproxLocal[optimalVN];
            double intervalPrice = integrateQuadraticWithNormalWeight(q.a, q.b, q.c, z1, z2);
            price += intervalPrice;
            // Probability for this interval
            double prob = 0.5 * (std::erfc(-z2 * M_SQRT1_2) - std::erfc(-z1 * M_SQRT1_2));
            probabilities[optimalVN] += prob;
        }
        // Print results
        std::cout << "\nDeliveryContract Price: " << price << "\n";
        std::cout << "Delivery Probabilities (decimal):\n";
        std::cout.precision(16);
        for (size_t i = 0; i < probabilities.size(); ++i) {
            std::cout << "ValueNote " << (i + 1) << ": " << probabilities[i] << "\n";
        }
    }

    // Prompt user for sigma and VP0 vectors, compute and print sensitivities using the correct chain-rule functions
    void promptAndPrintUserSensitivities() const {
        std::vector<double> sigmas(basket.size()), vp0s(basket.size());
        std::cout << "\nEnter volatility (sigma, in %) and today price (VP0) for each ValueNote:\n";
        for (size_t i = 0; i < basket.size(); ++i) {
            std::cout << "ValueNote " << (i + 1) << " (sigma, VP0): ";
            std::cin >> sigmas[i] >> vp0s[i];
        }
        std::cout << "\nUser-supplied Sensitivities (at input values):\n";
        std::cout.precision(16);
        for (size_t i = 0; i < basket.size(); ++i) {
            double volSens = volatilitySensitivityForNote(i, sigmas[i]);
            double priceSens = priceSensitivityForNote(i, vp0s[i]);
            std::cout << "ValueNote " << (i + 1) << ":\n";
            std::cout << "  sigma = " << sigmas[i] << ", VP0 = " << vp0s[i] << "\n";
            std::cout << "  Volatility Sensitivity: " << volSens << "\n";
            std::cout << "  Price Sensitivity: " << priceSens << "\n";
        }
    }

    // Print Q1.1-Q1.3 values for the first ValueNote in the basket, for all conventions
    void printQ1TableForFirstNote() const {
        if (basket.empty()) {
            std::cout << "No ValueNotes in basket.\n";
            return;
        }
        const ValueNote& note = basket[0];
        double ER0, VP0;
        std::cout << "\nEnter ER0 (effective rate, %) for Q1.1/Q1.2/Q1.3: ";
        std::cin >> ER0;
        std::cout << "Enter VP0 (price) for Q1.1/Q1.2/Q1.3: ";
        std::cin >> VP0;
        std::cout << "\n--- Q1.1 ---\n";
        std::cout << "Linear Price (ER0): " << note.eRateToPriceLinear(ER0) << "\n";
        std::cout << "Cumulative Price (ER0): " << note.eRateToPriceCumulative(ER0) << "\n";
        std::cout << "Recursive Price (ER0): " << note.eRateToPriceRecursive(ER0) << "\n";
        std::cout << "Linear Effective Rate (VP0): " << note.priceToERateLinear(VP0) << "\n";
        std::cout << "Cumulative Effective Rate (VP0): " << note.priceToERateCumulative(VP0) << "\n";
        std::cout << "Recursive Effective Rate (VP0): " << note.priceToERateRecursive(VP0) << "\n";
        std::cout << "\n--- Q1.2 ---\n";
        std::cout << "Linear Price Sensitivity (ER0): " << note.eRateToPriceSensitivityLinear(ER0) << "\n";
        std::cout << "Cumulative Price Sensitivity (ER0): " << note.eRateToPriceSensitivityCumulative(ER0) << "\n";
        std::cout << "Recursive Price Sensitivity (ER0): " << note.eRateToPriceSensitivityRecursive(ER0) << "\n";
        std::cout << "Linear Effective Rate Sensitivity (VP0): " << note.priceToERateSensitivityLinear(VP0) << "\n";
        std::cout << "Cumulative Effective Rate Sensitivity (VP0): " << note.priceToERateSensitivityCumulative(VP0) << "\n";
        std::cout << "Recursive Effective Rate Sensitivity (VP0): " << note.priceToERateSensitivityRecursive(VP0) << "\n";
        std::cout << "\n--- Q1.3 ---\n";
        std::cout << "Linear Price Second Derivative (ER0): " << note.eRateToPriceSecondDerivativeLinear(ER0) << "\n";
        std::cout << "Cumulative Price Second Derivative (ER0): " << note.eRateToPriceSecondDerivativeCumulative(ER0) << "\n";
        std::cout << "Recursive Price Second Derivative (ER0): " << note.eRateToPriceSecondDerivativeRecursive(ER0) << "\n";
        std::cout << "Linear Effective Rate Second Derivative (VP0): " << note.priceToERateSecondDerivativeLinear(VP0) << "\n";
        std::cout << "Cumulative Effective Rate Second Derivative (VP0): " << note.priceToERateSecondDerivativeCumulative(VP0) << "\n";
        std::cout << "Recursive Effective Rate Second Derivative (VP0): " << note.priceToERateSecondDerivativeRecursive(VP0) << "\n";
    }


};

int main() {
    double svr, rfr;
    int expTime, numNotes;
    std::string rfc, erc;

    std::cout << "Enter standardized value rate: ";
    std::cin >> svr;
    std::cout << "Enter expiration time(months): ";
    std::cin >> expTime;
    std::cout << "Enter relative factor convention (Cumulative or Unity): ";
    std::cin >> rfc;
    std::cout << "Enter effective rate convention: ";
    std::cin >> erc;
    std::cout << "Enter risk-free interest rate: ";
    std::cin >> rfr;
    std::cout << "Enter number of notes in the basket: ";
    std::cin >> numNotes;

    DeliveryContract dc(svr, expTime, rfc, erc, rfr, numNotes);
    try {
        dc.createBasket();
        std::cout << "[DEBUG] Finished createBasket()\n";
        dc.computeRelativeFactors();
        std::cout << "[DEBUG] Finished computeRelativeFactors()\n";
        dc.printQ1TableForFirstNote();
        std::cout << "[DEBUG] Finished printQ1TableForFirstNote()\n";
        dc.printRelativeFactors();
        std::cout << "[DEBUG] Finished printRelativeFactors()\n";
        dc.computeForwardPricesAndAdjustedRates();
        std::cout << "[DEBUG] Finished computeForwardPricesAndAdjustedRates()\n";
        dc.printForwardPrices();
        std::cout << "[DEBUG] Finished printForwardPrices()\n";
        dc.printRiskAdjustedEffectiveRates();
        std::cout << "[DEBUG] Finished printRiskAdjustedEffectiveRates()\n";
        dc.computeAndPrintPriceAndProbabilities();
        std::cout << "[DEBUG] Finished computeAndPrintPriceAndProbabilities()\n";
        dc.promptAndPrintUserSensitivities();
        std::cout << "[DEBUG] Finished promptAndPrintUserSensitivities()\n";
    } catch (const std::exception& ex) {
        std::cerr << "Exception: " << ex.what() << std::endl;
    }

    return 0;
}
