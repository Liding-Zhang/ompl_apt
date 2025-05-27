#include <vector>
#include <memory>
#include <queue>
#include <functional>
#include <algorithm>
#include <limits>

#include "ompl/geometric/planners/informedtrees/aptstar/TotalForce.h"
#include "ompl/geometric/planners/informedtrees/aptstar/State.h"
#include <ompl/base/spaces/RealVectorStateSpace.h>

namespace ompl
{
    namespace geometric
    {
        namespace aptstar
        {
            // Constructor for TotalForce class
            TotalForce::TotalForce(const std::shared_ptr<State> &state, std::vector<std::shared_ptr<State>> &states,
                                   size_t dimension,
                                   const double &radius)
              : state_(state), states_(states), dimension_(dimension),radius_(radius)
            {
            }

            // Destructor for TotalForce class
            TotalForce::~TotalForce()
            {
            }

            typedef std::vector<double> Vector;

            // Set the current state
            void TotalForce::setState(const std::shared_ptr<State> &state)
            {
                if (!state)
                {
                    throw std::invalid_argument("Provided state is null");
                }

                std::lock_guard<std::mutex> lock(state_mutex_);
                state_ = state;

                return;
            }

            void TotalForce::setUseUCBCharge(bool useUCBCharge)
            {
                useUCBCharge_ = useUCBCharge;
            }


            void TotalForce::setbatchsize(unsigned int numsamples)
            {
                batchsize_ = numsamples;
            }

            // Get the current state
            std::shared_ptr<State> TotalForce::getState() const
            {
                return state_;
            }

            // Get all states
            std::vector<std::shared_ptr<State>> TotalForce::getStates() const
            {
                return states_;
            }

            // Get vector between two states
            Vector TotalForce::getVector(const std::shared_ptr<State> state1, const std::shared_ptr<State> state2) const
            {
                if (!state1)
                {
                    throw std::invalid_argument("Provided state1 is null");
                }

                if (!state2)
                {
                    throw std::invalid_argument("Provided state2 is null");
                }

                auto rstate1 = state1->raw()->as<ompl::base::RealVectorStateSpace::StateType>();
                auto rstate2 = state2->raw()->as<ompl::base::RealVectorStateSpace::StateType>();

                if (!rstate1 || !rstate2)
                {
                    throw std::runtime_error("rstate pointer is null");
                }

                std::vector<double> vector(dimension_);

                for (size_t i = 0; i < dimension_; ++i)
                {
                    vector[i] = rstate2->values[i] - rstate1->values[i];
                }

                return vector;
            }

            // Calculate the distance between two states
            double TotalForce::distance(const std::shared_ptr<State> state1, const std::shared_ptr<State> state2) const
            {
                if (!state1)
                {
                    throw std::invalid_argument("Provided state1 is null");
                }

                if (!state2)
                {
                    throw std::invalid_argument("Provided state2 is null");
                }

                auto cstate1 = state1->raw()->as<ompl::base::RealVectorStateSpace::StateType>();
                auto cstate2 = state2->raw()->as<ompl::base::RealVectorStateSpace::StateType>();

                if (!cstate1 || !cstate2)
                {
                    throw std::runtime_error("rstate pointer is null");
                }
                double theta1 = 0., theta2 = 0., dx = 0., dy = 0., dist = 0.;

                for (unsigned int i = 0; i < dimension_; ++i)
                {
                    theta1 += cstate1->values[i];
                    theta2 += cstate2->values[i];
                    dx += cos(theta1) - cos(theta2);
                    dy += sin(theta1) - sin(theta2);
                    dist += sqrt(dx * dx + dy * dy);
                }

                return dist;
            }

            // Comparator for priority queue
            struct Compare
            {
                bool operator()(const std::pair<double, std::shared_ptr<State>> &a,
                                const std::pair<double, std::shared_ptr<State>> &b)
                {
                    return a.first > b.first;
                }
            };

            // Calculate force between two states
            std::vector<double> TotalForce::force(const std::shared_ptr<State> &state1,
                                                  const std::shared_ptr<State> &state2)
            {
                // Retrieve the vector between two states
                std::vector<double> vec = getVector(state1, state2);
                // Calculate the distance
                double dist = distance(state1, state2);
                double magnitude = 1.0 / (dist * dist);
                double charge = 1.0;
                if (useUCBCharge_)
                {
                    double normalizedBatch = (static_cast<double>(batchsize_)- 1) / (199 -1);

                    // charge = exponentialFunction(normalizedBatch, maxCharges_, minCharges_);
                    // charge = polynomialFunction(normalizedBatch, maxCharges_, minCharges_);
                    // charge = logarithmicFunction(normalizedBatch, maxCharges_, minCharges_);
                    charge = tanhFunction(normalizedBatch, maxCharges_, minCharges_, 100); // Adjust alpha = 100 for different smoothness
                    // charge = iteratFunction(normalizedBatch);

                    magnitude = charge / (dist * dist);
                    std::cout << "current charge: " << charge << std::endl;
                    // std::cout << "current normalizedBatch: " << normalizedBatch << std::endl;
                    // std::cout << "current batch: " << batchsize_ << std::endl;
                }
                else
                {

                    magnitude = charge / (dist * dist);
                }
                // Coulomb's Law: F = k * (q / r^2)
                magnitude =
                    state2->isAttractive_ ? magnitude : -magnitude;  // Attractive force positive, repulsive negative

                for (size_t i = 0; i < vec.size(); ++i)
                {
                    vec[i] *= magnitude;
                }

                return vec;
            }

            // Calculate the total force acting on the current state from all samples
            void TotalForce::totalForce(const std::shared_ptr<State> &currentState,
                                        const std::vector<std::shared_ptr<State>> &Samples)
            {
                std::vector<double> totalForceVec(dimension_, 0.0);

                if (Samples.empty())
                {
                    totalForceVec_ = totalForceVec;
                }
                else
                {
                    for (const auto &sample : Samples)
                    {
                        auto f = force(currentState, sample);
                        for (size_t i = 0; i < dimension_; ++i)
                        {
                            totalForceVec[i] += f[i];
                        }
                    }
                    double totalMagnitude = 0.0;
                    for (size_t i = 0; i < dimension_; ++i)
                    {
                        totalMagnitude += totalForceVec[i] * totalForceVec[i];
                    }
                    totalMagnitude = sqrt(totalMagnitude);
                    totalMagnitude_ = totalMagnitude;

                    for (double elem : totalForceVec)
                    {
                        if (std::abs(elem) > std::numeric_limits<double>::epsilon())
                        {
                            totalForceVec_ = normalize(totalForceVec);
                        }
                        else
                            totalForceVec_ = totalForceVec;
                    }
                }
             
            }

            std::vector<std::shared_ptr<State>> TotalForce::RNearestEllipseticSamples(
                const std::shared_ptr<State> state, std::vector<double> &totalForceVec,
                std::vector<std::shared_ptr<State>> &Samples, std::size_t dimension_)
            {
                std::vector<std::shared_ptr<State>> nearest;
            
                // Iterate over all samples
                for (const auto &sample : Samples)
                {
                    // Exclude the current state from consideration
                    if (sample != state)
                    {
                        // Add sample to the queue if it's within a finite distance
                        if (isinStretchedEllipse(sample, state, totalForceVec, dimension_))
                        {
                            nearest.push_back(sample);
                        }
                    }
                }
                return nearest;
            }

            void TotalForce::totalForcewithStart(const std::shared_ptr<State> &currentState,
                                                 const std::shared_ptr<State> &startstate,
                                                 const std::shared_ptr<State> &goalstate, bool iterateForwardSearch)
            {
                // Initialize force vectors for single force and total force with start
                std::vector<double> singleForceVec(dimension_, 0.0);
                std::vector<double> totalForceVecwithStart(dimension_, 0.0);
                std::vector<double> totalForceVec = totalForceVec_;

                // Calculate the force vector if iterating forward
                if (iterateForwardSearch)
                {
                    goalstate->setAttractive();
                    if (currentState == goalstate)
                    {
                        singleForceVec = singleForceVec;
                    }
                    else
                    {
                        auto f = force(currentState, goalstate);
                        singleForceVec = normalize(f);
                        for (auto &element : singleForceVec)
                        {
                            element *= 1;
                        }
                    }
                }
                else  // Calculate the force vector if iterating backward
                {
                    startstate->setAttractive();
                    if (currentState == startstate)
                    {
                        singleForceVec = singleForceVec;
                    }
                    else
                    {
                        auto f = force(currentState, startstate);
                        singleForceVec = normalize(f);
                        for (auto &element : singleForceVec)
                        {
                            element *= 1;
                        }
                    }
                }

                // Combine and normalize the total force vector with the start force vector
                totalForceVecwithStart = normalize(singleForceVec + totalForceVec);

                // Calculate vectors from goal to start and vice versa
                std::vector<double> vectorgoaltostart = getVector(goalstate, startstate);
                vectorgoaltostart = normalize(vectorgoaltostart);
                std::vector<double> vectorstarttogoal = getVector(startstate, goalstate);
                vectorstarttogoal = normalize(vectorstarttogoal);

                // Set the final total force vector with start value
                totalForceVecwithStart_ = totalForceVecwithStart;
            }

            std::vector<double> TotalForce::getValueofforceDirection()
            {
                // Return the total force vector with start value
                return totalForceVecwithStart_;
            }

            bool TotalForce::isinStretchedEllipse(const std::shared_ptr<State> &state,
                                                  const std::shared_ptr<State> &ellipsecenter,
                                                  std::vector<double> &Vector, 
                                                  std::size_t dimension_)
            {
                std::vector<double> delta(dimension_);
                auto cstate = state->raw()->as<ompl::base::RealVectorStateSpace::StateType>();
                auto cellipsecenter = ellipsecenter->raw()->as<ompl::base::RealVectorStateSpace::StateType>();
                // Step 1: calculate (x - c)
                for (size_t i = 0; i < dimension_; ++i) {
                    delta[i] = cstate->values[i] - cellipsecenter->values[i];
                }

                // Step 2: calculate Q^T * (x - c)
                auto Q = buildOrthogonalMatrix(Vector);
                std::vector<double> rotated(dimension_, 0.0);
                for (size_t i = 0; i < dimension_; ++i) {
                    for (size_t j = 0; j < dimension_; ++j) {
                        rotated[i] += Q[j][i] * delta[j]; 
                    }
                }

                // Step 3: calculate sum(rotated[i]^2 / d[i]^2)
                double value = 0.0;
                auto D = buildAxisLengths(Vector,1);
                for (size_t i = 0; i < dimension_; ++i) {
                    value += (rotated[i] * rotated[i]) / (D[i] * D[i]);
                }

                // Step 4: Determine whether it is within the ellipsoid
                return value <= 1.0;
            }

            // Normalized vector
            std::vector<double> TotalForce::normalize(const std::vector<double>& v) 
            {
                double norm = 0.0;
                for (double val : v) {
                    norm += val * val;
                }
                norm = std::sqrt(norm);

                std::vector<double> result(v.size());
                for (size_t i = 0; i < v.size(); ++i) {
                    result[i] = v[i] / norm;
                }
                return result;
            }                      

            // Construct an orthogonal matrix Q
            std::vector<std::vector<double>> TotalForce::buildOrthogonalMatrix(const std::vector<double>& force) 
            {
                size_t dim = force.size();
                std::vector<std::vector<double>> Q(dim, std::vector<double>(dim, 0.0));

                // Step 1: Normalized force direction
                std::vector<double> u1 = normalize(force);
                for (size_t i = 0; i < dim; ++i) {
                    Q[i][0] = u1[i];
                }

                // Step 2: Randomly generate remaining orthogonal vectors
                for (size_t k = 1; k < dim; ++k) {
                    std::vector<double> v(dim, 0.0);
                    v[k] = 1.0; 

                    // Orthogonalization
                    double dot = 0.0;
                    for (size_t i = 0; i < dim; ++i) {
                        dot += v[i] * Q[i][0];
                    }
                    for (size_t i = 0; i < dim; ++i) {
                        v[i] -= dot * Q[i][0];
                    }

                    // Normalization
                    v = normalize(v);

                    for (size_t i = 0; i < dim; ++i) {
                        Q[i][k] = v[i];
                    }
                }

                return Q;
            }

            // Construct semi-axis length
            std::vector<double> TotalForce::buildAxisLengths(const std::vector<double>& force, double k) 
            {
                size_t dim = force.size();
                std::vector<double> D(dim, radius_);

                double force_norm = 0.0;
                for (double val : force) {
                    force_norm += val * val;
                }
                force_norm = std::sqrt(force_norm);

                // Along the direction of force, stretch
                D[0] = radius_ * (1.0 + k * force_norm); 
                return D;
            }
   
            // ...

            double TotalForce::getNorm(const Vector &vector) const
            {
                // Check if the provided vector is empty and throw an exception if so
                if (vector.empty())
                {
                    throw std::invalid_argument("Provided vector is null");
                }

                // Calculate the norm (magnitude) of the vector
                double sum = 0.0;
                for (auto v : vector)
                {
                    sum += v * v;  // Sum of squares
                }

                return std::sqrt(sum);  // Square root of the sum
            }

            double TotalForce::dotProduct(const Vector &v1, const Vector &v2) const
            {
                // Check for null vectors and throw an exception if found
                if (v1.empty() || v2.empty())
                {
                    throw std::invalid_argument("One of the provided vectors is null");
                }

                // Calculate the dot product of two vectors
                double dot_product = 0.0;
                for (size_t i = 0; i < dimension_; ++i)
                {
                    dot_product += v1[i] * v2[i];  // Sum of products
                }

                return dot_product;
            }

            std::vector<double> TotalForce::vectorProjection(const std::vector<double> &a, const std::vector<double> &b)
            {
                // Calculate the projection of vector a onto vector b
                double dotAB = dotProduct(a, b);  // Dot product of a and b
                double dotBB = dotProduct(b, b);  // Dot product of b with itself
                double scale = dotAB / dotBB;     // Scaling factor

                // Apply the scale to each element of vector b
                std::vector<double> projection;
                projection.reserve(b.size());
                for (auto &element : b)
                {
                    projection.push_back(element * scale);
                }
                return projection;
            }

            Vector TotalForce::normalize(const Vector &v) const
            {
                // Calculate the norm of the vector
                double norm = getNorm(v);

                // Check for zero vector and throw an exception if found
                if (std::fabs(norm) < epsilon_)
                {
                    throw std::runtime_error("Cannot normalize a zero vector");
                }

                // Normalize the vector by dividing each element by the norm
                Vector normalizedVec = v;
                for (auto &val : normalizedVec)
                {
                    val /= norm;
                }
                return normalizedVec;
            }

            bool TotalForce::isVectorBetween(const Vector &targetVector, const Vector &goalVector,
                                             const std::shared_ptr<State> &sourceState,
                                             const std::shared_ptr<State> &neighborState) const
            {
                // Safety checks for null vectors and states
                if (targetVector.empty() || !sourceState || !neighborState)
                {
                    throw std::invalid_argument("One of the provided vectors or states is null");
                }

                // Calculate the vector to check
                Vector vectorTocheck = getVector(sourceState, neighborState);

                // Normalize the vectors to focus on direction
                Vector normalizedTarget = normalize(targetVector);
                Vector normalizedGoal = normalize(goalVector);
                Vector normalizedToCheck = normalize(vectorTocheck);

                // Compute the dot products
                double dotProductWithTarget = dotProduct(normalizedToCheck, normalizedTarget);
                double dotProductWithGoal = dotProduct(normalizedToCheck, normalizedGoal);

                // Return true if the vector is between targetVector and goalVector
                return (dotProductWithTarget >= 0.0 && dotProductWithGoal >= 0.0);
            }

            bool TotalForce::checkAngle(const Vector &targetVector, const std::shared_ptr<State> &sourceState,
                                        const std::shared_ptr<State> &neighborState) const
            {
                // Safety checks for null vectors and states
                if (targetVector.empty() || !sourceState || !neighborState)
                {
                    throw std::invalid_argument("One of the provided vectors or states is null");
                }

                // Calculate the vector to check
                Vector vectorTocheck = getVector(sourceState, neighborState);

                // Calculate the dot product of target vector and vector to check
                double dotValue = dotProduct(targetVector, vectorTocheck);

                // Check for zero magnitude and throw an exception if found
                if (std::fabs(getNorm(targetVector)) < epsilon_)
                {
                    throw std::invalid_argument("Target vector has zero magnitude.");
                }

                // Calculate the cosine of the angle between vectors
                double cosAngle = dotValue / (getNorm(targetVector) * getNorm(vectorTocheck));
                double safetyBuffer = std::cos(M_PI / 3.0);  // 60 degrees threshold

                // Return true if the cosine of the angle is less than the safety buffer
                return cosAngle <= safetyBuffer;
            }

            double TotalForce::getRatioofValidInvalidPoints(std::vector<std::shared_ptr<State>> states)
            {
                // Count the number of valid and invalid points
                std::size_t l1 = 0;  // Number of valid points
                std::size_t l2 = 0;  // Number of invalid points
                for (auto state : states)
                {
                    if (state->isAttractive_)
                    {
                        l1++;  // Count valid points
                    }
                    else
                    {
                        l2++;  // Count invalid points
                    }
                }

                // Check for zero valid points and return INFINITY if found
                if (l1 == 0)
                {
                    return INFINITY;
                }
                else
                {
                    double ratio = static_cast<double>(l2) / l1;
                    return ratio;  // Return the ratio of invalid to valid points
                }
            }
        }  // namespace aptstar
    }      // namespace geometric
}  // namespace ompl
