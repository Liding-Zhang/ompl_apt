#ifndef OMPL_GEOMETRIC_PLANNERS_INFORMEDTREES_APTSTAR_TOTALFORCE_
#define OMPL_GEOMETRIC_PLANNERS_INFORMEDTREES_APTSTAR_TOTALFORCE_

#include <memory>
#include <vector>
#include <mutex>

#include "ompl/base/Cost.h"
#include "ompl/base/OptimizationObjective.h"
#include "ompl/datastructures/BinaryHeap.h"

#include "ompl/geometric/planners/informedtrees/aptstar/Direction.h"
#include "ompl/geometric/planners/informedtrees/aptstar/Edge.h"
#include "ompl/geometric/planners/informedtrees/aptstar/Vertex.h"
namespace ompl
{
    namespace geometric
    {
        namespace aptstar
        {
            // Forward declare the APT* state class.
            class State;

            /** \brief The vertex class for both the forward and reverse search. */
            class TotalForce
            {
            public:
                /** \brief Constructs the vertex, which must be associated with a state. */
                TotalForce(const std::shared_ptr<State> &state, std::vector<std::shared_ptr<State>> &Allneigbors,
                           size_t dimension_, const double &radius);

                /** \brief Destructs this vertex. */
                ~TotalForce();
                
                /** \brief Set the state */
                void setState(const std::shared_ptr<State> &state);

                /** \brief Return the state */
                std::shared_ptr<State> getState() const;

                /** \brief Return the states */
                std::vector<std::shared_ptr<State>> getStates() const;

                /** \brief Return the std::vector<double> by two states */
                std::vector<double> getVector(const std::shared_ptr<State> state1,
                                              const std::shared_ptr<State> state2) const;

                /** \brief Returns distance between 2 states */
                double distance(const std::shared_ptr<State> state1, const std::shared_ptr<State> state2) const;

                /** \brief Returns all state points within the stretched ellipse */
                std::vector<std::shared_ptr<State>>RNearestEllipseticSamples(const std::shared_ptr<State> state, std::vector<double> &ForceDirection,
                                          std::vector<std::shared_ptr<State>> &Samples,std::size_t dimension_);

                /** \brief Determine whether the point is inside the stretched ellipse */
                bool isinStretchedEllipse(const std::shared_ptr<State> &state1, const std::shared_ptr<State> &state2, std::vector<double> &totalForceVec, std::size_t dimension_);

                /** \brief Normalized vector */
                std::vector<double> normalize(const std::vector<double>& v);

                /** \brief Construct an orthogonal matrix Q */
                std::vector<std::vector<double>> buildOrthogonalMatrix(const std::vector<double>& force);

                /** \brief Construct semi-axis length */
                std::vector<double> buildAxisLengths(const std::vector<double>& force, double k);

                /** \brief Returns the Coulomb force between two points */
                std::vector<double> force(const std::shared_ptr<State> &state1, const std::shared_ptr<State> &state2);

                /** \brief Returns the Coulomb force acting on the current state */
                void totalForce(const std::shared_ptr<State> &currentState, const std::vector<std::shared_ptr<State>> &Samples);

                /** \brief Returns the Coulomb force acting on the current state */
                void totalForcewithStart(const std::shared_ptr<State> &currentState, const std::shared_ptr<State> &startstate, const std::shared_ptr<State> &goalstate, bool iterateForwardSearch);

                /** \brief Return the norm of two states */
                double getNorm(const std::vector<double> &v) const;

                /** \brief Return the dot product of two Vectors */
                double dotProduct(const std::vector<double> &v1, const std::vector<double> &v2) const;

                std::vector<double> vectorProjection(const std::vector<double> &a, const std::vector<double> &b);

                /** \brief normalize the vector */
                std::vector<double> normalize(const std::vector<double> &v) const;

                /** \brief Check whether the std::vector<double> within two guided direction. */
                bool isVectorBetween(const std::vector<double> &targetVector, const std::vector<double> &goalVector, const std::shared_ptr<State> &sourceState, const std::shared_ptr<State> &neighborState) const;

                /** \brief Check whether the std::vector<double> meet the required angle. */
                bool checkAngle(const std::vector<double> &targetVector, const std::shared_ptr<State> &sourceState, const std::shared_ptr<State> &neighborState) const;

                /** \brief Sets the charge. */
                void setUseUCBCharge(bool useUCBCharge);

                /** \brief Sets the charge. */
                void setbatchsize(unsigned int numsamples);

                /** \brief filter the neighbors vector using direction info*/
                void filterNeighbors(const std::shared_ptr<State> &lastState, const std::shared_ptr<State> &state, std::vector<std::shared_ptr<State>> &states, const std::shared_ptr<State> &goalState);

                /** \brief Clear the neighbors vector of current target */
                void clearNeighbors();

                /** \brief Returns the proportion of invalid points */
                double getRatioofValidInvalidPoints(std::vector<std::shared_ptr<State>> states);

                /** \brief Set total force value */
                void settotalForcewithStartValue(std::vector<double> &totalForceVecWithStart);

                /** \brief Space dimension */
                std::size_t dimension_;

                /** \brief Total magnitude */
                mutable double totalMagnitude_;

                /** \brief Total force vector */
                mutable std::vector<double> totalForceVec_;

                /** \brief Total magnitude with start */
                mutable double totalMagnitudewithStart_;

                mutable std::vector<double> totalForceVecwithStart_;

                /** \brief small parameter */
                const double epsilon_ = 1e-9;

                std::vector<double> getValueofforceDirection();

                bool useUCBCharge_{true};

                unsigned int batchsize_;

                // Function representing an exponential decay
                double exponentialFunction(double ratio, double maxCharges, double minCharges) {
                    double decayFactor = 1 - std::exp(-5 * ratio); // Modify decay speed with constant (e.g., 5)
                    return maxCharges - (maxCharges - minCharges) * decayFactor;
                }

                // Function representing a polynomial decay
                double polynomialFunction(double ratio, double maxCharges, double minCharges, int power = 2) {
                    double decayFactor = std::pow(ratio, power); // Polynomial decay with adjustable power
                    return maxCharges - (maxCharges - minCharges) * decayFactor;
                }

                // Function representing a reciprocal decay
                double reciprocalFunction(double ratio, double maxCharges, double minCharges) {
                    if (ratio == 0) { // Avoid division by zero
                        return maxCharges;
                    }
                    double decayFactor = 1.0 / (1.0 + 10.0 * (1.0 - ratio)); // Adjust steepness with constant (e.g., 10)
                    return maxCharges - (maxCharges - minCharges) * decayFactor;
                }

                double iteratFunction(double iteration) {
                    if (iteration >= 0.1 && iteration < 0.3) {
                        return 1.9;
                    } else if (iteration >= 0.3 && iteration < 0.5) {
                        return 1.5;
                    } else if (iteration >= 0.5 && iteration < 0.8) {
                        return 0.8;
                    } else if (iteration >= 0.8 && iteration <= 1.0) {
                        return 0.1;
                    } else {
                        return 1.9;
                    }
                }

                // Function representing a logarithmic decay
                double logarithmicFunction(double ratio, double maxCharges, double minCharges) {
                    if (ratio == 0) { // Avoid log(0)
                        return minCharges;
                    }
                    double decayFactor = std::log(1.0 + 9.0 * ratio) / std::log(10.0); // Adjust base and multiplier
                    return maxCharges - (maxCharges - minCharges) * decayFactor;
                }

                // Function representing a hyperbolic tangent (tanh) decay
                // double tanhFunction(double ratio, double maxCharges, double minCharges) {
                //     double decayFactor = (std::tanh(6.0 * (ratio - 0.5)) + 1.0) / 2.0; // Shift and scale tanh
                //     return maxCharges - (maxCharges - minCharges) * decayFactor;
                // }

                double tanhFunction(double ratio, double maxCharges, double minCharges, int alpha) {
                    double decayFactor = 0.0;
                    double batchSize = ratio - 0.5; // Offset according to the formula
                    double chargeDiff = minCharges - maxCharges;
                    
                    // Taylor series expansion
                    for (int i = 1; i <= alpha; ++i) {
                        double B2i = bernoulliNumber(i);
                        double term = pow(2, 2 * i - 1) * B2i * (pow(2, 2 * i) - 1) * pow(batchSize, 2 * i - 1) / std::tgamma(2 * i + 1);
                        decayFactor += term;
                    }
                    
                    // Normalization: clamp decayFactor to [-1, 1]
                    // Scaling by tanh properties
                    double normalizedDecayFactor = std::tanh(decayFactor * 6.0);  // Adjust 6.0 to affect smoothness
                    
                    // Zoom and shift, keeping the range unchanged
                    return maxCharges - (maxCharges - minCharges) * (normalizedDecayFactor + 1) / 2.0;
                }
                double bernoulliNumber(int n) {
                    double Bn = 0.0;
                    for (int j = 0; j <= n; ++j) {
                        double innerSum = 0.0;
                        for (int k = 0; k <= j; ++k) {
                            innerSum += pow(-1, k) * (std::tgamma(j + 1) / (std::tgamma(k + 1) * std::tgamma(j - k + 1))) * pow(k, 2 * n);
                        }
                        Bn += innerSum / (j + 1);
                    }
                    return Bn;
                }
                
                // Assuming maxCharges and minCharges values for demonstration
                double maxCharges_ = 1.9;
                double minCharges_ = 0.1;

            private:
                /** \brief target state */
                std::shared_ptr<State> state_;

                /** \brief The neighbor states of target state */
                std::vector<std::shared_ptr<State>> states_;

                /** \brief The dimension of the state */

                /** \brief state lock */
                std::mutex state_mutex_;

                /** \brief The connection radius of the RGG. */
                double radius_{std::numeric_limits<double>::infinity()};
            };

        }  // namespace aptstar

    }  // namespace geometric

}  // namespace ompl
inline std::vector<double> operator-(const std::vector<double> &lhs, const std::vector<double> &rhs)
{
    assert(lhs.size() == rhs.size() && "Vectors must be of the same size for subtraction!");

    std::vector<double> result(lhs.size());
    for (size_t i = 0; i < lhs.size(); ++i)
    {
        result[i] = lhs[i] - rhs[i];
    }
    return result;
}

inline std::vector<double> operator+(const std::vector<double> &lhs, const std::vector<double> &rhs)
{
    assert(lhs.size() == rhs.size() && "Vectors must be of the same size for addition!");

    std::vector<double> result(lhs.size());
    for (size_t i = 0; i < lhs.size(); ++i)
    {
        result[i] = lhs[i] + rhs[i];
    }
    return result;
}
inline std::vector<double> operator*(const std::vector<double> &lhs, double scalar)
{
    std::vector<double> result(lhs.size());
    for (size_t i = 0; i < lhs.size(); ++i)
    {
        result[i] = lhs[i] * scalar;
    }
    return result;
}
inline std::vector<double> operator/(const std::vector<double> &lhs, double scalar)
{
    assert(scalar != 0 && "Cannot divide by zero!");

    std::vector<double> result(lhs.size());
    for (size_t i = 0; i < lhs.size(); ++i)
    {
        result[i] = lhs[i] / scalar;
    }
    return result;
}

#endif  // OMPL_GEOMETRIC_PLANNERS_INFORMEDTREES_APTSTAR_TOTALFORCE_
