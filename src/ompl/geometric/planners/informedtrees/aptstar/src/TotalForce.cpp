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
                                   size_t dimension)
              : state_(state), states_(states), dimension_(dimension)
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

            // Find nearest K samples to a given state
            std::vector<std::shared_ptr<State>> TotalForce::NearestKSamples(const std::shared_ptr<State> state,
                                                                            std::vector<std::shared_ptr<State>> Samples,
                                                                            int k)
            {
                std::priority_queue<std::pair<double, std::shared_ptr<State>>,
                                    std::vector<std::pair<double, std::shared_ptr<State>>>, Compare>
                    pq;

                // Iterate through all samples, calculate distance to the given state, and store in priority queue
                for (const auto &sample : Samples)
                {
                    double dist = distance(state, sample);
                    pq.push({dist, sample});
                }

                // Extract nearest k samples from the queue
                std::vector<std::shared_ptr<State>> nearest;
                for (int i = 0; i < k && !pq.empty(); ++i)
                {
                    nearest.push_back(pq.top().second);
                    pq.pop();
                }

                return nearest;
            }

            // Calculate force between two states
            std::vector<double> TotalForce::force(const std::shared_ptr<State> &state1,
                                                  const std::shared_ptr<State> &state2)
            {
                // Retrieve the vector between two states
                std::vector<double> vec = getVector(state1, state2);
                // Calculate the distance
                double dist = distance(state1, state2);
                double magnitude;
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
                // Coulomb's Law: F = k * (1 / r^2)
                magnitude = state2->isAttractive_ ? magnitude : -magnitude;  // Attractive force positive, repulsive negative
                // std::cout << "magnitude: " << magnitude << std::endl;
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
            // Find the nearest k samples in an ellipsoidal space
            std::vector<std::shared_ptr<State>> TotalForce::NearestEllipseticKSamples(
                const std::shared_ptr<State> state, std::vector<double> &totalForceVec,
                std::vector<std::shared_ptr<State>> &Samples, int k, std::size_t dimension_)
            {
                // Initialize a priority queue to store samples based on their distance
                std::priority_queue<std::pair<double, std::shared_ptr<State>>,
                                    std::vector<std::pair<double, std::shared_ptr<State>>>, Compare>
                    queue;

                // Iterate over all samples
                for (const auto &sample : Samples)
                {
                    // Exclude the current state from consideration
                    if (sample != state)
                    {
                        double dist = calculateEllipticalDistance(state, sample, totalForceVec, dimension_);
                        // Add sample to the queue if it's within a finite distance
                        if (dist <= INFINITY)
                        {
                            queue.push(std::make_pair(dist, sample));
                        }
                    }
                }

                // Extract the nearest k samples from the queue
                std::vector<std::shared_ptr<State>> nearest;
                int positiveCount = 0;  // Count of positive (attractive) samples found

                while (!queue.empty() && positiveCount < k)
                {
                    auto currentSample = queue.top().second;
                    queue.pop();
                    // Count only attractive samples
                    if (currentSample->isAttractive_)
                    {
                        positiveCount++;
                    }
                    nearest.push_back(currentSample);
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

            // Function to calculate the elliptical distance between two states
            double TotalForce::calculateEllipticalDistance(const std::shared_ptr<State> &state1,
                                                           const std::shared_ptr<State> &state2,
                                                           std::vector<double> &Vector, std::size_t dimension_)
            {
                // If the force vector is empty, initialize it as a unit circle
                if (Vector.empty())
                {
                    std::vector<double> circle(dimension_, 1.0);
                    Vector = circle;
                }

                // Check for null states and throw an exception if found
                if (!state1 || !state2)
                {
                    throw std::invalid_argument("Provided state is null");
                }

                // Convert states to real vector state space types
                auto cstate1 = state1->raw()->as<ompl::base::RealVectorStateSpace::StateType>();
                auto cstate2 = state2->raw()->as<ompl::base::RealVectorStateSpace::StateType>();

                // Calculate the elliptical distance between the two states
                double theta1 = 0., theta2 = 0., dx = 0., dy = 0., dist = 0.;
                for (unsigned int i = 0; i < dimension_; ++i)
                {
                    theta1 += cstate1->values[i];
                    theta2 += cstate2->values[i];
                    dx += cos(theta1) - cos(theta2);
                    dy += sin(theta1) - sin(theta2);
                    dist += sqrt((dx * dx + dy * dy) / (Vector[i] * Vector[i]));
                }

                return dist;
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
