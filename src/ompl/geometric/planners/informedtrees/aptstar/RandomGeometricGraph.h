/*********************************************************************
 * Software License Agreement (BSD License)
 *
 *  Copyright (c) 2023, Technical University of Munich
 *  All rights reserved.
 *
 *  Redistribution and use in source and binary forms, with or without
 *  modification, are permitted provided that the following conditions
 *  are met:
 *
 *   * Redistributions of source code must retain the above copyright
 *     notice, this list of conditions and the following disclaimer.
 *   * Redistributions in binary form must reproduce the above
 *     copyright notice, this list of conditions and the following
 *     disclaimer in the documentation and/or other materials provided
 *     with the distribution.
 *   * Neither the name of the University of Munich nor the names of its
 *     contributors may be used to endorse or promote products derived
 *     from this software without specific prior written permission.
 *
 *  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 *  "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 *  LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FDITNESS
 *  FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
 *  COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
 *  INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
 *  BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 *  LOSS OF USE, DATA, OR PROFDITS; OR BUSINESS INTERRUPTION) HOWEVER
 *  CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
 *  LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
 *  ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 *  POSSIBILITY OF SUCH DAMAGE.
 *********************************************************************/

// Authors: Marlin Strub, Liding Zhang, Xu Liang

#ifndef OMPL_GEOMETRIC_PLANNERS_INFORMEDTREES_APTSTAR_RANDOM_GEOMETRIC_GRAPH_
#define OMPL_GEOMETRIC_PLANNERS_INFORMEDTREES_APTSTAR_RANDOM_GEOMETRIC_GRAPH_

#include <limits>
#include <memory>
#include <iostream>  // This is needed for ompl's nearest neighbors struct.

#include "ompl/datastructures/NearestNeighborsGNATNoThreadSafety.h"
#include "ompl/base/samplers/InformedStateSampler.h"
#include "ompl/base/Planner.h"
#include "ompl/base/ProblemDefinition.h"
#include "ompl/base/OptimizationObjective.h"

#include "ompl/geometric/planners/informedtrees/aptstar/State.h"
#include "ompl/geometric/planners/informedtrees/aptstar/Edge.h"

namespace ompl
{
    namespace geometric
    {
        namespace aptstar
        {
            class RandomGeometricGraph
            {
            public:
                /** \brief Constructs a random geometric graph with the given space information and reference to the
                 * current solution cost. */
                RandomGeometricGraph(const std::shared_ptr<ompl::base::SpaceInformation> &spaceInfo,
                                     const ompl::base::Cost &solutionCost);

                /** \brief Destricts this random geometric graph. */
                ~RandomGeometricGraph() = default;

                /** \brief Setup the graph with the given problem definition and planner input states. */
                void setup(const std::shared_ptr<ompl::base::ProblemDefinition> &problem,
                           ompl::base::PlannerInputStates *inputStates);

                /** \brief Clears all internal planner structures but retains settings. Subsequent calls to solve() will
                 * start from scratch. */
                void clear();

                /** \brief Clears all query-specific structures, such as start and goal states. */
                void clearQuery();

                /** \brief Adds new starts and goals to the graph if available and creates a new informed sampler if
                 * necessary. */
                void updateStartAndGoalStates(const ompl::base::PlannerTerminationCondition &terminationCondition,
                                              ompl::base::PlannerInputStates *inputStates);

                /** \brief Returns the minimum possible cost for the current problem, using admissible cost estimates to
                 * calculate it. */
                ompl::base::Cost minPossibleCost() const;

                /** \brief Sets the radius factor (eta in the paper). */
                void setRadiusFactor(double factor);

                /** \brief Returns the radius factor (eta in the paper). */
                double getRadiusFactor() const;

                /** \brief Sets the effort threshold. */
                void setEffortThreshold(const unsigned int threshold);

                /** \brief Gets the effort threshold. */
                unsigned int getEffortThreshold() const;

                std::vector<std::shared_ptr<State>> getValidSamples();

                std::vector<std::shared_ptr<State>> getInValidSamples();

                std::vector<double> getTotalForceDirection();

                void setBuffer();

                void setInValidbuffer();

                /** \brief Enable pruning of the graph. */
                void enablePruning(bool prune);

                /** \brief Returns Whether pruning is enabled. */
                bool isPruningEnabled() const;

                /** \brief Enable multiquery usage of the graph. */
                void enableMultiquery(bool multiquery);

                /** \brief Returns Whether multiquery usage of the graph is enabled. */
                bool isMultiqueryEnabled() const;

                /** \brief Set whether to use a k-nearest connection model. If false, it uses an r-disc model. */
                void setUseKNearest(bool useKNearest);

                void setAllUseKNearest(bool useAllKNearest);

                void setUseRNearest(bool useRNearest);

                void nearestKEllipse(std::shared_ptr<State> state, std::numeric_limits<std::size_t> k,
                                     std::vector<std::shared_ptr<State>> neighbors);

                void nearestK(std::shared_ptr<State> state, std::numeric_limits<std::size_t> k,
                              std::vector<std::shared_ptr<State>> neighbors);

                /** \brief Returns whether the graph uses a k-nearest connection model. If false, it uses an r-disc
                 * model. */
                bool getUseKNearest() const;

                /** \brief Sets the maximum number of goals APT* will sample from sampleable goal regions. */
                void setMaxNumberOfGoals(unsigned int maxNumberOfGoals);

                /** \brief Sets the charge. */
                void setUseUCBCharge(bool useUCBCharge);

                /** \brief Returns the charge. */
                bool getUseUCBCharge();

                /** \brief Returns the maximum number of goals APT* will sample from sampleable goal regions. */
                unsigned int getMaxNumberOfGoals() const;

                /** \brief Samples random states and adds them to the graph. */
                bool addStates(std::size_t numStates,
                               const ompl::base::PlannerTerminationCondition &terminationCondition);

                bool addVAlidandInvalidStates(std::size_t numNewStates,
                                              const ompl::base::PlannerTerminationCondition &terminationCondition);

                /** \brief Prunes the graph of states that can not improve the current solution. */
                void prune();

                /** \brief Returns the neighbors of a state. */
                std::vector<std::weak_ptr<State>> getNeighbors(const std::shared_ptr<State> &state) const;

                std::vector<std::weak_ptr<State>> getAllNeighbors(const std::shared_ptr<State> &state,
                                                                  const std::shared_ptr<State> &startState,
                                                                  const std::shared_ptr<State> &goalState,
                                                                  bool iterateForwardSearch) const;

                /** \brief Returns the start states. */
                const std::vector<std::shared_ptr<State>> &getStartStates() const;

                /** \brief Returns the goal states. */
                const std::vector<std::shared_ptr<State>> &getGoalStates() const;

                /** \brief Returns the number of sampled states. */
                unsigned int getNumberOfSampledStates() const;

                /** \brief Returns the number of valid samples. */
                unsigned int getNumberOfValidSamples() const;

                /** \brief Returns the number of nearest neighbor calls. */
                unsigned int getNumberOfNearestNeighborCalls() const;

                /** \brief Sets the start state. */
                std::shared_ptr<State> registerStartState(const ompl::base::State *start);

                /** \brief Sets the goal state. */
                std::shared_ptr<State> registerGoalState(const ompl::base::State *goal);

                /** \brief Registers a whitelisted state */
                void registerWhitelistedState(const std::shared_ptr<State> &state) const;

                /** \brief Returns whether a start state is available. */
                bool hasStartState() const;

                /** \brief Returns whether a goal state is available. */
                bool hasGoalState() const;

                /** \brief Determine whether it is inside the stretched ellipse. */
                bool isInEllipticalRange(const std::shared_ptr<State> &state, const std::shared_ptr<State> &ellipsecenter, std::vector<double> &Vector, std::size_t dimension_) const;

                /** \brief Normalized vector */
                std::vector<double> normalize(const std::vector<double>& v) const;

                /** \brief Construct an orthogonal matrix Q */
                std::vector<std::vector<double>> buildOrthogonalMatrix(const std::vector<double>& force) const;

                /** \brief Construct semi-axis length */
                std::vector<double> buildAxisLengths(const std::vector<double>& force, double k) const;

                /** \brief Returns the RNN radius r */
                double getRadius() const;

                /** \brief Returns whether the given state is a start state. */
                bool isStart(const std::shared_ptr<State> &state) const;

                /** \brief Returns whether the given state is a goal state. */
                bool isGoal(const std::shared_ptr<State> &state) const;

                /** \brief Returns all sampled states (that have not been pruned). */
                std::vector<std::shared_ptr<State>> getStates() const;

                /** \brief Registers an invalid edge. */
                void registerInvalidEdge(const Edge &edge) const;

                /** \brief Returns the tag of the current RGG. */
                std::size_t getTag() const;

                /** \brief Set the seed used by the RNG and the StateSampler. The state sampler must already be
                 * allocated, as a new state sampler will not take this seed. */
                void setLocalSeed(std::uint_fast32_t localSeed)
                {
                    if (sampler_)
                    {
                        sampler_->setLocalSeed(localSeed);
                    }
                };

                /** \brief Returns the inadmissible effort to come. */
                unsigned int inadmissibleEffortToCome(const std::shared_ptr<State> &state) const;

                /** \brief Set the direction of the force. */
                void setForceDirection(std::vector<double> &totalforceDirection) const;

                /** \brief Direction of force. */
                mutable std::vector<double> totalforceDirection_;

                /** \brief All sampling points. */
                std::vector<State> allSamples_{};

               

            private:
                /** \brief Returns a sample either from the buffer or a newly generated one. */
                std::shared_ptr<State> getNewSample(const ompl::base::PlannerTerminationCondition &terminationCondition);

                /** \brief Returns all valid and invalid points. */
                std::shared_ptr<State> getValidandInValidNewSample();

                /** \brief Returns the number of states in the informed set. */
                std::size_t countSamplesInInformedSet() const;

                /** \brief Returns whether a state can be pruned because it cannot possibly be part of a solution equal
                 * to or better than the current solution. */
                bool canBePruned(const std::shared_ptr<State> &state) const;

                /** \brief Decides what start/goal vertices should become part of the graph permanently. */
                void pruneStartsAndGoals();

                /** \brief Returns the heuristic cost from the preferred start of a state. */
                ompl::base::Cost lowerBoundCostToCome(const std::shared_ptr<State> &state) const;

                /** \brief Returns the admissible effort from the preferred start of a state. */
                unsigned int lowerBoundEffortToCome(const std::shared_ptr<State> &state) const;

                /** \brief Returns the heuristic cost to the preferred goal of a state. */
                ompl::base::Cost lowerBoundCostToGo(const std::shared_ptr<State> &state) const;

                /** \brief Initializes the given state's cost and effort values. */
                void initializeState(const std::shared_ptr<State> &state);

                /** \brief Returns the number of neighbors of the k-nearest model with a given number of samples. */
                std::size_t computeNumberOfNeighbors(std::size_t numInformedSamples) const;

                /** \brief Returns the radius for the RGG. */
                double computeRadius(std::size_t numInformedSamples) const;

                /** \brief The list of whitelisted states. */
                mutable std::vector<std::shared_ptr<State>> whitelistedStates_{};

                /** \brief The tag of the current RGG. */
                std::size_t tag_{1u};

                /** \brief The buffered sampled states. */
                std::vector<std::shared_ptr<State>> buffer_;

                std::vector<std::shared_ptr<State>> InValidbuffer_;

                /** \brief The buffered start and goal states. */
                std::vector<std::shared_ptr<State>> startGoalBuffer_;

                /** \brief The current number of samples that are in use */
                std::size_t currentNumSamples_{0u};

                /** \brief The sampled states in a nearest neighbor structure. */
                NearestNeighborsGNATNoThreadSafety<std::shared_ptr<State>> samples_;

                /** \brief The new sampled states. */
                std::vector<std::shared_ptr<State>> newSamples_;

                /** \brief The state sampler. */
                std::shared_ptr<ompl::base::InformedSampler> sampler_{nullptr};

                std::shared_ptr<ompl::base::InformedSampler> sampler2_{nullptr};

                /** \brief The info about the underlying state space. */
                std::shared_ptr<ompl::base::SpaceInformation> spaceInfo_;

                /** \brief The underlying state space. */
                std::shared_ptr<ompl::base::StateSpace> space_;

                /** \brief The problem this graph is supposed to help solve. */
                std::shared_ptr<ompl::base::ProblemDefinition> problem_;

                /** \brief The optimization objective this graph is supposed to help optimize. */
                std::shared_ptr<ompl::base::OptimizationObjective> objective_;

                /** \brief The start states of the problem. */
                std::vector<std::shared_ptr<State>> startStates_;

                /** \brief The goal states of the problem. */
                std::vector<std::shared_ptr<State>> goalStates_;

                /** \brief The pruned start states of the problem. We keep these around because if a new goal is added
                 * after pruning a start, we might want to consider the pruned start again. */
                std::vector<std::shared_ptr<State>> prunedStartStates_;

                /** \brief The pruned goal states of the problem. We keep these around because if a new start is added
                 * after pruning a goal, we might want to consider the pruned goal again. */
                std::vector<std::shared_ptr<State>> prunedGoalStates_;

                /** \brief Whether pruning is enabled. */
                bool isPruningEnabled_{true};

                /** \brief Whether multiquery is enabled. */
                bool isMultiqueryEnabled_{false};

                /** \brief Whether to use a k-nearest RGG. If false, APT* uses an r-disc RGG. */
                bool useKNearest_{false};

                bool useRNearest_{false};

                bool useAllKNearest_{true};

                bool useUCBCharge_{false};

                unsigned int batchsize_;

                /** \brief The maximum number of goals APT* will sample explicitly from a sampleable goal region. */
                unsigned int maxNumGoals_{1u};

                /** \brief The number of neighbors that defines the neighborhood of a vertex if using a k-nearest graph.
                 */
                std::size_t numNeighbors_{std::numeric_limits<std::size_t>::max()};

                /** \brief A constant for the computation of the number of neighbors when using a k-nearest model. */
                std::size_t k_rgg_{std::numeric_limits<std::size_t>::max()};

                /** \brief The connection radius of the RGG. */
                 double radius_{std::numeric_limits<double>::infinity()};

                /** \brief The factor by which to scale the connection radius (eta in the paper). */
                double radiusFactor_{1.001};

                /** \brief The threshold which we use to decide if we keep a start/goal vertex. */
                unsigned int effortThreshold_{50000};

                /** \brief The dimension of the state space this graph is embedded in. */
                const double dimension_;

                /** \brief The measure of a unit ball in n dimensions. */
                const double unitNBallMeasure_;

                /** \brief The cost of the incumbent solution. */
                const ompl::base::Cost &solutionCost_;

                /** \brief The minimum possible cost for this problem. */
                ompl::base::Cost minPossibleCost_{std::numeric_limits<double>::signaling_NaN()};

                /** \brief The number of sampled states. */
                mutable unsigned int numSampledStates_{0u};

                /** \brief The number of valid samples. */
                mutable unsigned int numValidSamples_{0u};

                /** \brief The number of invalid samples. */
                mutable unsigned int numInVaildSamples_{0u};

                /** \brief The number of valid samples. */
                mutable unsigned int numNearestNeighborCalls_{0u};
            };

        }  // namespace aptstar

    }  // namespace geometric

}  // namespace ompl

#endif  // OMPL_GEOMETRIC_PLANNERS_INFORMEDTREES_APTSTAR_RANDOM_GEOMETRIC_GRAPH_
