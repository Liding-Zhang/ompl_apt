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
 *  LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
 *  FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
 *  COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
 *  INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
 *  BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 *  LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
 *  CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
 *  LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
 *  ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 *  POSSIBILITY OF SUCH DAMAGE.
 *********************************************************************/

// Authors: Marlin Strub, Liding Zhang, Sicheng Wang

#include "ompl/geometric/planners/informedtrees/APTstar.h"

#include <algorithm>
#include <memory>

#include "ompl/base/objectives/PathLengthOptimizationObjective.h"
#include "ompl/geometric/PathGeometric.h"
#include "ompl/geometric/planners/informedtrees/aptstar/AdaptiveBatchSize.h"
#include "ompl/geometric/planners/informedtrees/aptstar/TotalForce.h"

#include <fstream>
#include <sys/stat.h>

using namespace std::string_literals;
using namespace ompl::geometric::aptstar;

namespace ompl
{
    namespace geometric
    {
        APTstar::APTstar(const std::shared_ptr<ompl::base::SpaceInformation> &spaceInfo)
          : ompl::base::Planner(spaceInfo, "APT*")
          , graph_(spaceInfo, solutionCost_)
          , detectionState_(spaceInfo->allocState())
          , space_(spaceInfo->getStateSpace())
          , motionValidator_(spaceInfo->getMotionValidator())
          , solutionCost_()
        {
            // Specify APT*'s planner specs.
            specs_.recognizedGoal = base::GOAL_SAMPLEABLE_REGION;
            specs_.multithreaded = false;
            specs_.approximateSolutions = true;
            specs_.optimizingPaths = true;
            specs_.directed = true;
            specs_.provingSolutionNonExistence = false;
            specs_.canReportIntermediateSolutions = true;

            // Register the setting callbacks.
            declareParam<bool>("use_all_k_nearest", this, &APTstar::setUseKNearest, &APTstar::getUseKNearest, "0,1");
            declareParam<bool>("use_adaptive_batchsize", this, &APTstar::setUseAdaptiveBatchSize,
                               &APTstar::getUseAdaptiveBatchSize, "0,1");
            declareParam<std::size_t>("batch_size", this, &APTstar::setBatchSize, &APTstar::getBatchSize, "1:1:10000");
            declareParam<bool>("use_UCB_charge", this, &APTstar::setUseUCBCharge,
                               &APTstar::getUseUCBCharge, "0,1");
            declareParam<double>("rewire_factor", this, &APTstar::setRadiusFactor, &APTstar::getRadiusFactor,
                                 "1.0:0.01:3.0");
            declareParam<bool>("use_graph_pruning", this, &APTstar::enablePruning, &APTstar::isPruningEnabled, "0,1");
            declareParam<bool>("find_approximate_solutions", this, &APTstar::trackApproximateSolutions,
                               &APTstar::areApproximateSolutionsTracked, "0,1");
            declareParam<unsigned int>("set_max_num_goals", this, &APTstar::setMaxNumberOfGoals,
                                       &APTstar::getMaxNumberOfGoals, "1:1:1000");

            // Register the progress properties.
            addPlannerProgressProperty("iterations INTEGER", [this]() { return std::to_string(iteration_); });
            addPlannerProgressProperty("best cost DOUBLE", [this]() { return std::to_string(solutionCost_.value()); });
            addPlannerProgressProperty("state collision checks INTEGER",
                                       [this]() { return std::to_string(graph_.getNumberOfSampledStates()); });
            addPlannerProgressProperty("edge collision checks INTEGER",
                                       [this]() { return std::to_string(numCollisionCheckedEdges_); });
            addPlannerProgressProperty("nearest neighbour calls INTEGER",
                                       [this]() { return std::to_string(graph_.getNumberOfNearestNeighborCalls()); });
        }

        APTstar::~APTstar()
        {
            spaceInfo_->freeState(detectionState_);
        }

        void APTstar::setup()
        {
            // Call the base class setup.
            Planner::setup();

            // Check that the problem definition is set.
            if (static_cast<bool>(problem_))
            {
                // Default to path length optimization if no objective has been specified.
                if (!problem_->hasOptimizationObjective())
                {
                    OMPL_WARN("%s: No optimization objective has been specified. The default is optimizing path "
                              "length.",
                              name_.c_str());
                    problem_->setOptimizationObjective(
                        std::make_shared<ompl::base::PathLengthOptimizationObjective>(spaceInfo_));
                }

                if (static_cast<bool>(problem_->getGoal()))
                {
                    // If we were given a goal, make sure its of appropriate type.
                    if (!(problem_->getGoal()->hasType(ompl::base::GOAL_SAMPLEABLE_REGION)))
                    {
                        OMPL_ERROR("APT* is currently only implemented for goals that can be cast to "
                                   "ompl::base::GOAL_SAMPLEABLE_REGION.");
                        setup_ = false;
                        return;
                    }
                }

                // Pull through the optimization objective for direct access.
                objective_ = problem_->getOptimizationObjective();

                // Initialize costs to infinity.
                solutionCost_ = objective_->infiniteCost();
                reverseCost_ = objective_->infiniteCost();
                approximateSolutionCost_ = objective_->infiniteCost();
                approximateSolutionCostToGoal_ = objective_->infiniteCost();

                // Instantiate the queues.
                forwardQueue_ = std::make_unique<aptstar::ForwardQueue>(objective_, space_);
                if (isMultiqueryEnabled_)
                {
                    reverseQueue_ = std::make_unique<aptstar::ReverseQueue>(objective_, space_,
                                                                            std::isfinite(suboptimalityFactor_));
                }
                else
                {
                    reverseQueue_ = std::make_unique<aptstar::ReverseQueue>(objective_, space_, true);
                }

                // Setup the graph with the problem information.
                graph_.setup(problem_, &pis_);

                // Create the start vertices.
                for (const auto &start : graph_.getStartStates())
                {
                    startVertices_.emplace_back(start->asForwardVertex());
                }

                // Create the goal vertices.
                for (const auto &goal : graph_.getGoalStates())
                {
                    goalVertices_.emplace_back(goal->asReverseVertex());
                }

                // Populate the queues.
                expandGoalVerticesIntoReverseQueue();
                expandStartVerticesIntoForwardQueue();
            }
            else
            {
                setup_ = false;
                OMPL_WARN("%s: Unable to setup without a problem definition.", name_.c_str());
            }
        }

        ompl::base::PlannerStatus APTstar::solve(const ompl::base::PlannerTerminationCondition &terminationCondition)
        {
            // Check that the planner and state space are setup.
            auto status = ensureSetup();

            // Return early if the planner or state space are not setup.
            if (status == ompl::base::PlannerStatus::StatusType::ABORT)
            {
                return status;
            }

            // Ensure that the problem has start and goal states before solving.
            status = ensureStartAndGoalStates(terminationCondition);

            // Return early if no problem can be solved.
            if (status == ompl::base::PlannerStatus::StatusType::INVALID_START ||
                status == ompl::base::PlannerStatus::StatusType::INVALID_GOAL)
            {
                return status;
            }

            // Let the world know what we're doing.
            OMPL_INFORM("%s: Solving the given planning problem. The current best solution cost is %.4f", name_.c_str(),
                        solutionCost_.value());

            // Iterate while we should continue solving the problem.
            while (continueSolving(terminationCondition))
            {
                iterate(terminationCondition);
            }

            // for (auto p : areaList_)
            //     std::cout << "ratio: " << p.first << " times: " << p.second << std::endl;

            // Someone might call ProblemDefinition::clearSolutionPaths() between invocations of Planner::sovle(), in
            // which case previously found solutions are not registered with the problem definition anymore.
            status = updateSolution();

            // Let the caller know the status.
            informAboutPlannerStatus(status);
            return status;
        }

        void APTstar::clear()
        {
            if (forwardQueue_)
            {
                forwardQueue_->clear();
            }

            if (reverseQueue_)
            {
                reverseQueue_->clear();
            }

            startVertices_.clear();
            goalVertices_.clear();
            objective_.reset();
            graph_.clear();

            // Reset the solution costs. Cannot use infiniteCost() before resetting the objective because the objective
            // of a new problem definition objective might define that differently than the old.
            solutionCost_ = ompl::base::Cost(std::numeric_limits<double>::signaling_NaN());
            reverseCost_ = ompl::base::Cost(std::numeric_limits<double>::signaling_NaN());
            approximateSolutionCost_ = ompl::base::Cost(std::numeric_limits<double>::signaling_NaN());
            approximateSolutionCostToGoal_ = ompl::base::Cost(std::numeric_limits<double>::signaling_NaN());

            // Reset the tags.
            reverseSearchTag_ = 1u;
            startExpansionGraphTag_ = 0u;
            numSparseCollisionChecksCurrentLevel_ = initialNumSparseCollisionChecks_;
            numSparseCollisionChecksPreviousLevel_ = 0u;
            suboptimalityFactor_ = std::numeric_limits<double>::infinity();

            Planner::clear();
            setup_ = false;
        }

        void APTstar::clearQuery()
        {
            if (setup_)
            {
                forwardQueue_->clear();
                reverseQueue_->clear();
                startVertices_.clear();
                goalVertices_.clear();
                graph_.clearQuery();
                solutionCost_ = objective_->infiniteCost();
                reverseCost_ = objective_->infiniteCost();
                approximateSolutionCost_ = objective_->infiniteCost();
                approximateSolutionCostToGoal_ = objective_->infiniteCost();

                suboptimalityFactor_ = std::numeric_limits<double>::infinity();
                restartReverseSearch();

                numSparseCollisionChecksCurrentLevel_ = initialNumSparseCollisionChecks_;

                numProcessedEdges_ = 0u;

                reverseSearchTag_++;

                if (!isMultiqueryEnabled_)
                {
                    clear();
                }
                setup_ = false;
            }
        }

        ompl::base::Cost APTstar::bestCost() const
        {
            return solutionCost_;
        }

        void APTstar::setBatchSize(unsigned int numSamples)
        {
            batchSize_ = numSamples;
        }

        unsigned int APTstar::getBatchSize() const
        {
            return batchSize_;
        }

        void APTstar::setInitialNumberOfSparseCollisionChecks(std::size_t numChecks)
        {
            initialNumSparseCollisionChecks_ = numChecks;
            numSparseCollisionChecksCurrentLevel_ = numChecks;
            numSparseCollisionChecksPreviousLevel_ = 0u;
        }

        void APTstar::setRadiusFactor(double factor)
        {
            graph_.setRadiusFactor(factor);
        }

        double APTstar::getRadiusFactor() const
        {
            return graph_.getRadiusFactor();
        }

        void APTstar::adjustRadiusFactor()
        {
            if (std::isinf(solutionCost_.value()))
            {
                return graph_.setRadiusFactor(graph_.getRadiusFactor());
            }

            double a = solutionCost_.value() / 2;
            double c = graph_.minPossibleCost().value() / 2;
            double b = std::sqrt(a * a - c * c);
            double S = M_PI * a * b;
            static bool pragma = false;
            if (!pragma)
            {
                S_max_initial_ = S;
                pragma = true;
            }

            double ratio = S / S_max_initial_;

            double lambda =
                100 / spaceInfo_->getStateDimension();  // Adjust this value to make the decay faster or slower
            // Sigmoid function to smooth ratio
            double smoothedValue = 1 / (1 + exp(-10 * (ratio - 0.5)));
            double decay_factor = log(1 + lambda * smoothedValue) / log(1 + lambda);

            double rewireFactor = 0.101 + (2.001 - 0.101) * decay_factor;
            std::cout << "---------------radiusFactor of current batch is : " << rewireFactor << std::endl;
            return graph_.setRadiusFactor(rewireFactor);
        }

        void APTstar::setSuboptimalityFactor(double factor)
        {
            suboptimalityFactor_ = factor;
        }

        void APTstar::enableMultiquery(bool multiquery)
        {
            isMultiqueryEnabled_ = multiquery;
            graph_.enableMultiquery(multiquery);
        };

        bool APTstar::isMultiqueryEnabled() const
        {
            return isMultiqueryEnabled_;
        };

        void APTstar::setStartGoalPruningThreshold(unsigned int threshold)
        {
            graph_.setEffortThreshold(threshold);
        }

        unsigned int APTstar::getStartGoalPruningThreshold() const
        {
            return graph_.getEffortThreshold();
        }

        void APTstar::enablePruning(bool enable)
        {
            graph_.enablePruning(enable);
        }

        bool APTstar::isPruningEnabled() const
        {
            return graph_.isPruningEnabled();
        }

        void APTstar::trackApproximateSolutions(bool track)
        {
            trackApproximateSolutions_ = track;
        }

        bool APTstar::areApproximateSolutionsTracked() const
        {
            return trackApproximateSolutions_;
        }

        void APTstar::setUseKNearest(bool useKNearest)
        {
            graph_.setUseKNearest(useKNearest);
        }

        bool APTstar::getUseKNearest() const
        {
            return graph_.getUseKNearest();
        }

        void APTstar::setUseAdaptiveBatchSize(bool useAdaptiveBatchSize)
        {
            useAdaptiveBatchSize_ = useAdaptiveBatchSize;
        }

        bool APTstar::getUseAdaptiveBatchSize() const
        {
            return useAdaptiveBatchSize_;
        }

        void APTstar::setUseUCBCharge(bool useUCBCharge)
        {
            useUCBCharge_ = useUCBCharge;
        }

        bool APTstar::getUseUCBCharge() const
        {
            return useUCBCharge_;
        }

        void APTstar::setMaxNumberOfGoals(unsigned int numberOfGoals)
        {
            graph_.setMaxNumberOfGoals(numberOfGoals);
        }

        unsigned int APTstar::getMaxNumberOfGoals() const
        {
            return graph_.getMaxNumberOfGoals();
        }

        bool APTstar::isForwardQueueEmpty() const
        {
            assert(forwardQueue_);
            return forwardQueue_->empty();
        }

        bool APTstar::isReverseQueueEmpty() const
        {
            assert(reverseQueue_);
            return reverseQueue_->empty();
        }

        std::vector<Edge> APTstar::getForwardQueue() const
        {
            return forwardQueue_->getEdges();
        }

        std::vector<Edge> APTstar::getReverseQueue() const
        {
            return reverseQueue_->getEdges();
        }

        std::vector<Edge> APTstar::getReverseTree() const
        {
            // Prepare the return value.
            std::vector<Edge> edges;

            // Define a helper that recursively gets all reverse edges of a vertex.
            const std::function<void(const std::shared_ptr<Vertex> &)> getEdgesRecursively =
                [&edges, &getEdgesRecursively](const std::shared_ptr<Vertex> &vertex)
            {
                for (const auto &child : vertex->getChildren())
                {
                    getEdgesRecursively(child);
                }
                // Catch the root case of the recursion.
                if (auto parent = vertex->getParent().lock())
                {
                    edges.emplace_back(parent->getState(), vertex->getState());
                }
            };

            // Get the edges of all reverse roots recursively.
            for (const auto &root : goalVertices_)
            {
                getEdgesRecursively(root);
            }

            // Return all edges in the reverse tree.
            return edges;
        }

        Edge APTstar::getNextForwardEdge() const
        {
            assert(forwardQueue_);
            if (forwardQueue_->empty())
            {
                return {};
            }
            return forwardQueue_->peek(suboptimalityFactor_);
        }

        Edge APTstar::getNextReverseEdge() const
        {
            assert(reverseQueue_);
            if (reverseQueue_->empty())
            {
                return {};
            }
            return reverseQueue_->peek();
        }

        void APTstar::getPlannerData(base::PlannerData &data) const
        {
            // Get the base class data.
            Planner::getPlannerData(data);

            // Add the samples and their outgoing edges.
            for (const auto &state : graph_.getStates())
            {
                // Add the state as a vertex.
                data.addVertex(base::PlannerDataVertex(state->raw(), state->getId()));

                // If the sample is in the forward tree, add the outgoing edges.
                for (const auto &state2 : graph_.getStates())
                {
                    if (state->isWhitelisted(state2))
                    {
                        data.addEdge(base::PlannerDataVertex(state->raw(), state->getId()),
                                     base::PlannerDataVertex(state2->raw(), state2->getId()));
                    }
                }
            }
        }

        void APTstar::setLocalSeed(std::uint_fast32_t localSeed)
        {
            graph_.setLocalSeed(localSeed);
        }

        void APTstar::iterate(const ompl::base::PlannerTerminationCondition &terminationCondition)
        {
            // If we are in a multiquery setting, we do not want to search the approximation
            // only consisting of start/goals, since this completely ignores the computational effort we have already
            // invested Thus, the first thing we do in this instance is adding the first batch of samples.

            // Adjust the batch size according to the decay method
            if (useAdaptiveBatchSize_)
            {
                if (solutionCost_.value() != lastsolutionCost_.value())
                // if (true)
                {
                    lastsolutionCost_.setValue(solutionCost_.value());
                    DecayMethod decay_method_ = DecayMethod::LOG;
                    const double minPossibleCost = graph_.minPossibleCost().value();
                    AdaptiveBatchSize adaptiveBatchSize_ =
                        AdaptiveBatchSize(decay_method_, solutionCost_, minPossibleCost, batchSize_, S_max_initial_,
                                          S_min_initial_, maxSamples_, minSamples_, spaceInfo_->getStateDimension());
                    unsigned int numSamples = adaptiveBatchSize_.adjustBatchSize(decay_method_);
                    setBatchSize(numSamples);
                    // if (solutionCost_.value() != lastsolutionCost_.value())
                    // {
                    //     lastsolutionCost_.setValue(solutionCost_.value());
                    // adjustRadiusFactor();
                    // }
                }

                if (useUCBCharge_)
                {
                    graph_.setUseUCBCharge(useUCBCharge_);
                }
                
            }

            if (isMultiqueryEnabled_ &&
                graph_.getStates().size() == graph_.getStartStates().size() + graph_.getGoalStates().size())
            {
                improveApproximation(terminationCondition);
                ++iteration_;
                return;
            }

            // First check if the reverse search needs to be continued.
            if (continueReverseSearch())
            {
                // std::cout << graph_.getNumberOfSampledStates() << "hhhhhhhhhhhhhhhhhhhhhhh" << std::endl;
                // std::cout << graph_.getNumberOfValidSamples() << "ggggggggggggggggggggg" << std::endl;

                std::vector<std::shared_ptr<State>> buffer = graph_.getValidSamples();
                std::vector<std::shared_ptr<State>> inValidBuffer = graph_.getInValidSamples();

                iterateReverseSearch();
            }  // If the reverse search is suspended, check if the forward search needs to be continued.
            else if (continueForwardSearch())
            {
                iterateForwardSearch();
            }  // If neither the reverse nor the forward search needs to be continued, improve the approximation.
            else
            {
                improveApproximation(terminationCondition);
            }
            // Increment the iteration count.
            ++iteration_;
        }

        void APTstar::iterateForwardSearch()
        {
            // Ensure the forward queue is not empty.
            assert(!forwardQueue_->empty());

            // Get the top edge from the queue.
            auto edge = forwardQueue_->pop(suboptimalityFactor_);
            ++numProcessedEdges_;

            // Assert the validity of the edge.
            assert(edge.source->hasForwardVertex());
            assert(!std::isfinite(suboptimalityFactor_) || isClosed(edge.target->asReverseVertex()));

            // The edge is a freeby if its parent is already the parent of the child.
            if (isInForwardTree(edge))
            {
                // We expand the target into the queue unless the target is a goal.
                forwardQueue_->insertOrUpdate(expandUnlessGoal(edge.target));
                return;
            }

            // If the edge can not improve the forward tree, then we're done with it.
            if (!couldImproveForwardTree(edge))
            {
                return;
            }

            // The edge could possibly improve the tree, check if it is valid.
            if (isValid(edge))
            {
                // Compute the true edge cost and the target cost through this edge.
                const auto edgeCost = objective_->motionCost(edge.source->raw(), edge.target->raw());
                const auto targetCost = combine(edge.source->getCurrentCostToCome(), edgeCost);

                // Check if the edge can actually improve the forward path and tree.
                if (isBetter(targetCost, edge.target->getCurrentCostToCome()) &&
                    isBetter(combine(targetCost, edge.target->getAdmissibleCostToGo()), solutionCost_))
                {
                    // Convenience access to parent and child vertices.
                    auto source = edge.source->asForwardVertex();
                    auto target = edge.target->asForwardVertex();

                    // Update the parent of the child in the forward tree.
                    target->updateParent(source);

                    // Add the child to the parents children.
                    source->addChild(target);

                    // Set the edge cost associated with this parent.
                    target->setEdgeCost(edgeCost);

                    // Update the cost-to-come.
                    edge.target->setCurrentCostToCome(targetCost);

                    // Update the cost of the children.
                    const auto changedVertices = target->updateCurrentCostOfChildren(objective_);

                    // Reflect changes in queue and solution cost.
                    for (const auto &vertex : changedVertices)
                    {
                        // Update any edge in the queue.
                        forwardQueue_->updateIfExists({vertex->getParent().lock()->getState(), vertex->getState()});

                        // Update the solution if the vertex is a goal.
                        if (graph_.isGoal(vertex->getState()))
                        {
                            updateExactSolution(vertex->getState());
                        }
                    }

                    // Expand the outgoing edges into the queue unless this state is the goal state.
                    if (!graph_.isGoal(edge.target))
                    {
                        auto edges = expand(edge.target, true);
                        edges.erase(std::remove_if(edges.begin(), edges.end(),
                                                   [&edge](const auto &e)
                                                   { return e.target->getId() == edge.source->getId(); }),
                                    edges.end());
                        forwardQueue_->insertOrUpdate(edges);
                    }
                    else  // It is the goal state, update the solution.
                    {
                        updateExactSolution(edge.target);
                    }
                }
            }
            else
            {
                // Check if the edge is used in the reverse tree.
                const bool inReverseTree = edge.source->asReverseVertex()->isParent(edge.target->asReverseVertex()) ||
                                           edge.target->asReverseVertex()->isParent(edge.source->asReverseVertex());

                // If this edge was in the reverse search tree, it could be updated.
                if (inReverseTree)
                {
                    // Remember the current number of collision checks and increase it.
                    numSparseCollisionChecksPreviousLevel_ = numSparseCollisionChecksCurrentLevel_;
                    numSparseCollisionChecksCurrentLevel_ = (2u * numSparseCollisionChecksPreviousLevel_) + 1u;

                    // Restart the reverse search.
                    restartReverseSearch();

                    // Rebuild the forward queue.
                    forwardQueue_->rebuild();
                }
            }
        }

        void APTstar::iterateReverseSearch()
        {
            // Ensure the reverse queue is not empty.
            assert(!reverseQueue_->empty());

            // Get the top edge from the queue.
            auto edge = reverseQueue_->pop();
            auto &source = edge.source;
            auto &target = edge.target;

            // The parent vertex must have an associated vertex in the tree.
            assert(source->hasReverseVertex());

            // Register the expansion of its parent.
            source->asReverseVertex()->registerExpansionInReverseSearch(reverseSearchTag_);

            // The edge is a freeby if its parent is already the parent of the child.
            if (isInReverseTree(edge))
            {
                auto outgoingEdges = expand(target, false);
                outgoingEdges.erase(std::remove_if(outgoingEdges.begin(), outgoingEdges.end(),
                                                   [&source](const auto &e)
                                                   { return e.target->getId() == source->getId(); }),
                                    outgoingEdges.end());

                // If there are no outoing edges from the target state, we can nosider it expanded.
                if (outgoingEdges.empty())
                {
                    target->asReverseVertex()->registerExpansionInReverseSearch(reverseSearchTag_);
                }
                else  // If there are outgoing edge, add them to or update them in the reverse queue.
                {
                    reverseQueue_->insertOrUpdate(outgoingEdges);
                }

                return;
            }

            // Check whether the edge could be valid.
            if (couldBeValid(edge))
            {
                // Compute the heuristic cost.
                const auto edgeCost = objective_->motionCostHeuristic(source->raw(), target->raw());

                // std::cout << "-------------- edge cost is : " << edgeCost.value() << std::endl;

                // Incorporate the edge in the reverse tree if it provides an improvement.
                const auto effort = estimateEffortToTarget(edge);
                const bool doesDecreaseEffort = (effort < target->getEstimatedEffortToGo());

                if ((!isMultiqueryEnabled_ && doesImproveReverseTree(edge, edgeCost)) ||
                    (isMultiqueryEnabled_ &&
                     ((std::isfinite(suboptimalityFactor_) && doesImproveReverseTree(edge, edgeCost)) ||
                      (!std::isfinite(suboptimalityFactor_) && doesDecreaseEffort))))
                {
                    // Get the parent and child vertices.
                    auto parentVertex = source->asReverseVertex();
                    auto childVertex = target->asReverseVertex();

                    // The child must not be closed.
                    assert(!isClosed(childVertex));

                    // Update the parent of the child in the reverse tree.
                    childVertex->updateParent(parentVertex);

                    // Add the child to the children of the parent.
                    parentVertex->addChild(childVertex);

                    // Update the admissible cost to go.
                    target->setAdmissibleCostToGo(objective_->betterCost(
                        combine(source->getAdmissibleCostToGo(), edgeCost), edge.target->getAdmissibleCostToGo()));

                    // Update the best cost estimate of the target state if this edge can improve it.
                    target->setEstimatedCostToGo(
                        objective_->betterCost(estimateCostToTarget(edge), target->getEstimatedCostToGo()));

                    // Update the best effort estimate of the target state if this edge can improve it.
                    target->setEstimatedEffortToGo(std::min(effort, target->getEstimatedEffortToGo()));
                    //  std::cout << "-------------" << effort << std::endl;

                    // If this edge improves the reverse cost, update it.
                    if (graph_.isStart(target) && isBetter(target->getAdmissibleCostToGo(), reverseCost_))
                    {
                        reverseCost_ = target->getAdmissibleCostToGo();
                    }

                    // Update any edge in the forward queue affected by the target.
                    for (const auto &queueSource : target->getSourcesOfIncomingEdgesInForwardQueue())
                    {
                        forwardQueue_->updateIfExists({queueSource.lock(), target});
                    }

                    // Expand the target state into the reverse queue.
                    auto outgoingEdges = expand(target, false);

                    outgoingEdges.erase(
                        std::remove_if(outgoingEdges.begin(), outgoingEdges.end(),
                                       [&source, this](const auto &e)
                                       {
                                           if (e.target->getId() == source->getId())
                                           {
                                               return true;
                                           }
                                           if (objective_->isFinite(solutionCost_) &&
                                               isBetter(solutionCost_, reverseQueue_->computeAdmissibleSolutionCost(e)))
                                           {
                                               return true;
                                           }

                                           return false;
                                       }),
                        outgoingEdges.end());

                    // If there are no outoing edges from the target state, we can nosider it expanded.
                    if (outgoingEdges.empty())
                    {
                        target->asReverseVertex()->registerExpansionInReverseSearch(reverseSearchTag_);
                    }
                    else  // If there are outgoing edge, add them to or update them in the reverse queue.
                    {
                        reverseQueue_->insertOrUpdate(outgoingEdges);
                    }
                }
            }
        }

        void APTstar::improveApproximation(const ompl::base::PlannerTerminationCondition &terminationCondition)
        {
            // Add new states, also prunes states if enabled. The method returns true if all states have been added.
            if (graph_.addVAlidandInvalidStates(batchSize_, terminationCondition))
            {
                // // Reset the reverse collision detection.
                // std::cout << "__________________" << graph_.getValidSamples().size() << std::endl;

                // std::cout << "__________________" << graph_.getInValidSamples().size() << std::endl;
                numSparseCollisionChecksCurrentLevel_ = initialNumSparseCollisionChecks_;
                numSparseCollisionChecksPreviousLevel_ = 0u;

                // Restart the reverse search.
                restartReverseSearch();

                // Reinitialize the forward queue.
                forwardQueue_->clear();
                expandStartVerticesIntoForwardQueue();
            }
        }

        ompl::base::PlannerStatus::StatusType APTstar::ensureSetup()
        {
            // Call the base planners validity check. This checks if the
            // planner is setup if not then it calls setup().
            checkValidity();

            // Ensure the planner is setup.
            if (!setup_)
            {
                OMPL_ERROR("%s: Called solve without setting up the planner first.", name_.c_str());
                return ompl::base::PlannerStatus::StatusType::ABORT;
            }

            // Ensure the space is setup.
            if (!spaceInfo_->isSetup())
            {
                OMPL_ERROR("%s: Called solve without setting up the state space first.", name_.c_str());
                return ompl::base::PlannerStatus::StatusType::ABORT;
            }

            return ompl::base::PlannerStatus::StatusType::UNKNOWN;
        }

        ompl::base::PlannerStatus::StatusType
        APTstar::ensureStartAndGoalStates(const ompl::base::PlannerTerminationCondition &terminationCondition)
        {
            // If the graph currently does not have a start state, try to get one.
            if (!graph_.hasStartState())
            {
                graph_.updateStartAndGoalStates(terminationCondition, &pis_);

                // If we could not get a start state, then there's nothing to solve.
                if (!graph_.hasStartState())
                {
                    OMPL_ERROR("%s: No solution can be found as no start states are available", name_.c_str());
                    return ompl::base::PlannerStatus::StatusType::INVALID_START;
                }
            }

            // If the graph currently does not have a goal state, we wait until we get one.
            if (!graph_.hasGoalState())
            {
                graph_.updateStartAndGoalStates(terminationCondition, &pis_);

                // If the graph still doesn't have a goal after waiting, then there's nothing to solve.
                if (!graph_.hasGoalState())
                {
                    OMPL_ERROR("%s: No solution can be found as no goal states are available", name_.c_str());
                    return ompl::base::PlannerStatus::StatusType::INVALID_GOAL;
                }
            }

            // Would it be worth implementing a 'setup' or 'checked' status type?
            return ompl::base::PlannerStatus::StatusType::UNKNOWN;
        }

        void APTstar::finalPathOptimize(std::vector<std::shared_ptr<State>> &states)
        {
            size_t start = 0;
            while (start < states.size() - 1)
            {
                size_t end = states.size() - 1;
                while (end > start + 1)
                {
                    if (isValid(Edge(states[start], states[end])))
                    {
                        // Remove all states between start and end
                        states.erase(states.begin() + start + 1, states.begin() + end);
                        break;
                    }
                    end--;
                }
                start++;
            }
        }

        std::shared_ptr<ompl::geometric::PathGeometric>
        APTstar::getPathToState(const std::shared_ptr<aptstar::State> &state)
        {
            // Allocate a vector for states. The append function of the path inserts states in front of an
            // std::vector, which is not very efficient. I'll rather iterate over the vector in reverse.
            std::vector<std::shared_ptr<State>> states;
            auto current = state;

            // Collect all states in reverse order of the path (starting from the goal).
            while (!graph_.isStart(current))
            {
                assert(current->asForwardVertex()->getParent().lock());
                states.emplace_back(current);
                current = current->asForwardVertex()->getParent().lock()->getState();
            }
            states.emplace_back(current);
            // Optimize the final path
            // finalPathOptimize(states);

            // Append all states to the path in correct order (starting from the start).
            auto path = std::make_shared<ompl::geometric::PathGeometric>(spaceInfo_);
            for (auto it = states.crbegin(); it != states.crend(); ++it)
            {
                assert(*it);
                assert((*it)->raw());
                path->append((*it)->raw());
            }

            return path;
        }

        bool APTstar::continueSolving(const ompl::base::PlannerTerminationCondition &terminationCondition) const
        {
            // We stop solving the problem if:
            //   - The termination condition is satisfied; or
            //   - The current solution satisfies the objective; or
            //   - There is no better solution to be found for the current start goal pair and no new starts or
            //     goals are available
            // We continue solving the problem if we don't stop solving.
            return !(terminationCondition || objective_->isSatisfied(solutionCost_) ||
                     (!isBetter(graph_.minPossibleCost(), solutionCost_) && !pis_.haveMoreStartStates() &&
                      !pis_.haveMoreGoalStates()));
        }

        unsigned int APTstar::getForwardEffort() const
        {
            if (forwardQueue_->size() != 0u)
            {
                const auto forwardEdge = forwardQueue_->peek(suboptimalityFactor_);
                const auto forwardEffort = forwardQueue_->estimateEffort(forwardEdge);

                return forwardEffort;
            }

            return std::numeric_limits<unsigned int>::infinity();
        }

        bool APTstar::continueReverseSearch() const
        {
            // Never continue the reverse search if the reverse or forward queues are empty. There's nothing to continue
            // in the reverse search, if the reverse queue is empty. There's nothing to continue in the forward search,
            // if the forward queue is empty and continuing the reverse search cannot insert edges in the forward queue.
            if (reverseQueue_->empty() || forwardQueue_->empty())
            {
                return false;
            }

            /*

            If multiquery-planning is enabled, and we are in the process of searching for the initial solution,
            there are two conditions under which we suspend the reverse search:

               1. The best edge in the forward search has a closed target (i.e. there is a path to an edge in the
                  forward queue). Since we expand the reverse queue ordered by effort in the case that the suboptimality
                  factor is infdite, and multiquery planning is enabled, this is the lowest effort path.

               2. We already found a possible path, and there is no way that a lower effort path exists.

            If multiquery-planning is not enabled, there are three conditions under which the reverse search can be
            suspended:

               1. The best edge in the forward search has a closed target (admissible cost-to-go estimate), and the
                  reverse search cannot lead to a better solution than the potential solution of this edge.

               2. All edges in the forward queue have closed targets (admissible cost-to-go estimates).

               3. We do not care about solution cost and the least-effort edge in the forward queue is connected to the
                  reverse tree.

             */

            const auto forwardEdge = forwardQueue_->peek(suboptimalityFactor_);

            if (!std::isfinite(suboptimalityFactor_) && isMultiqueryEnabled_)
            {
                if (!forwardEdge.target->hasReverseVertex())
                {
                    return true;
                }

                const auto forwardEffort = forwardQueue_->estimateEffort(forwardEdge);
                const auto reverseEffort = reverseQueue_->peekEffort();
                unsigned int minEffortToCome = 0u;

                const auto reverseEdge = reverseQueue_->peek();
                if (!reverseEdge.target->hasForwardVertex())
                {
                    minEffortToCome = forwardQueue_->getMinEffortToCome();
                }

                return reverseEffort + minEffortToCome < forwardEffort;
            }

            const bool condition1 = isClosed(forwardEdge.target->asReverseVertex()) &&
                                    isBetter(forwardQueue_->getLowerBoundOnOptimalSolutionCost(),
                                             reverseQueue_->getLowerBoundOnOptimalSolutionCost());

            const bool condition2 = !forwardQueue_->containsOpenTargets(reverseSearchTag_);

            const bool condition3 = !std::isfinite(suboptimalityFactor_) && forwardEdge.target->hasReverseVertex();

            // The reverse search must be continued if it cannot be suspended.
            return !(condition1 || condition2 || condition3);
        }

        bool APTstar::continueForwardSearch() const
        {
            // Never continue to forward search if the forward queue is empty.
            if (forwardQueue_->empty())
            {
                return false;
            }

            // The forward search must be continued if the potential solution cost of the best edge is lower than the
            // current solution cost.
            return isBetter(forwardQueue_->getLowerBoundOnOptimalSolutionCost(), solutionCost_);
        }

        void APTstar::updateExactSolution()
        {
            for (const auto &goal : graph_.getGoalStates())
            {
                if (goal->hasForwardVertex())
                {
                    updateExactSolution(goal);
                }
            }
        }

        ompl::base::PlannerStatus::StatusType APTstar::updateSolution()
        {
            updateExactSolution();
            if (objective_->isFinite(solutionCost_))
            {
                return ompl::base::PlannerStatus::StatusType::EXACT_SOLUTION;
            }
            else if (trackApproximateSolutions_)
            {
                updateApproximateSolution();
                return ompl::base::PlannerStatus::StatusType::APPROXIMATE_SOLUTION;
            }
            else
            {
                return ompl::base::PlannerStatus::StatusType::TIMEOUT;
            }
        }

        void APTstar::restartReverseSearch()
        {
            reverseQueue_->clear();

            if (isMultiqueryEnabled_)
            {
                reverseQueue_->setCostQueueOrder(std::isfinite(suboptimalityFactor_));
            }
            else
            {
                reverseQueue_->setCostQueueOrder(true);
            }

            goalVertices_.clear();
            reverseCost_ = objective_->infiniteCost();
            for (const auto &goal : graph_.getGoalStates())
            {
                goalVertices_.emplace_back(goal->asReverseVertex());
                goal->setAdmissibleCostToGo(objective_->identityCost());
                goal->setEstimatedCostToGo(objective_->identityCost());
                goal->setEstimatedEffortToGo(0u);
            }
            expandGoalVerticesIntoReverseQueue();
            ++reverseSearchTag_;
        }

        void APTstar::updateApproximateSolution()
        {
            for (const auto &start : graph_.getStartStates())
            {
                start->asForwardVertex()->callOnBranch([this](const std::shared_ptr<aptstar::State> &state) -> void
                                                       { updateApproximateSolution(state); });
            }
        }

        void APTstar::updateCurrentCostToCome(const std::shared_ptr<aptstar::State> &state)
        {
            // There is no updating to do if the state is a start.
            if (graph_.isStart(state))
            {
                return;
            }

            // If the state is not in the forward tree, then its current cost to come is infinity.
            if (!state->hasForwardVertex())
            {
                state->setCurrentCostToCome(objective_->infiniteCost());
                return;
            }

            // If the state is in the forward tree and not the start, then update its cost.
            auto forwardVertex = state->asForwardVertex();
            state->setCurrentCostToCome(combine(forwardVertex->getParent().lock()->getState()->getCurrentCostToCome(),
                                                forwardVertex->getEdgeCost()));
        }

        void APTstar::updateExactSolution(const std::shared_ptr<aptstar::State> &goal)
        {
            // We update the current goal if
            //   1. The new goal has a better cost to come than the old goal
            //   2. Or the exact solution we found is no longer registered with the problem definition
            if (isBetter(goal->getCurrentCostToCome(), solutionCost_) || !problem_->hasExactSolution())
            {
                if (!std::isfinite(suboptimalityFactor_))
                {
                    suboptimalityFactor_ = 1.0;
                }

                // Update the best cost.
                solutionCost_ = goal->getCurrentCostToCome();

                // Register this solution with the problem definition.
                ompl::base::PlannerSolution solution(getPathToState(goal));
                solution.setPlannerName(name_);
                solution.setOptimized(objective_, solutionCost_, objective_->isSatisfied(solutionCost_));
                problem_->addSolutionPath(solution);

                if (!std::isfinite(suboptimalityFactor_))
                {
                    // If we found this solution with a suboptimality factor greater than 1, set the factor to one now.
                    // Empirically, this results in faster convergence, see associated publication for more info.
                    suboptimalityFactor_ = 1.0;

                    if (isMultiqueryEnabled_)
                    {
                        restartReverseSearch();
                        // Reinitialize the forward queue.
                        forwardQueue_->clear();
                        expandStartVerticesIntoForwardQueue();
                    }
                }

                // Let the user know about the new solution.
                informAboutNewSolution();
            }
        }

        void APTstar::updateApproximateSolution(const std::shared_ptr<aptstar::State> &state)
        {
            assert(trackApproximateSolutions_);
            if ((state->hasForwardVertex() || graph_.isStart(state)) && !graph_.isGoal(state))
            {
                const auto costToGoal = computeCostToGoToGoal(state);
                if (isBetter(costToGoal, approximateSolutionCostToGoal_) || !problem_->hasSolution())
                {
                    approximateSolutionCost_ = state->getCurrentCostToCome();
                    approximateSolutionCostToGoal_ = costToGoal;
                    ompl::base::PlannerSolution solution(getPathToState(state));
                    solution.setPlannerName(name_);

                    // Set the approximate flag.
                    solution.setApproximate(costToGoal.value());

                    // This solution is approximate and can not satisfy the objective.
                    solution.setOptimized(objective_, approximateSolutionCost_, false);

                    // Let the problem definition know that a new solution exists.
                    pdef_->addSolutionPath(solution);
                }
            }
        }

        ompl::base::Cost APTstar::computeCostToGoToGoal(const std::shared_ptr<aptstar::State> &state) const
        {
            auto bestCost = objective_->infiniteCost();
            for (const auto &goal : graph_.getGoalStates())
            {
                bestCost = objective_->betterCost(bestCost, objective_->motionCost(state->raw(), goal->raw()));
            }
            return bestCost;
        }

        void APTstar::informAboutNewSolution() const
        {
            OMPL_INFORM("%s (%u iterations): Found a new exact solution of cost %.4f. Sampled a total of %u states, %u "
                        "of which were valid samples (%.1f \%). Processed %u edges, %u of which were collision checked "
                        "(%.1f \%). The forward search tree has %u vertices. The reverse search tree has %u vertices.",
                        name_.c_str(), iteration_, solutionCost_.value(), graph_.getNumberOfSampledStates(),
                        graph_.getNumberOfValidSamples(),
                        graph_.getNumberOfSampledStates() == 0u ?
                            0.0 :
                            100.0 * (static_cast<double>(graph_.getNumberOfValidSamples()) /
                                     static_cast<double>(graph_.getNumberOfSampledStates())),
                        numProcessedEdges_, numCollisionCheckedEdges_,
                        numProcessedEdges_ == 0u ? 0.0 :
                                                   100.0 * (static_cast<float>(numCollisionCheckedEdges_) /
                                                            static_cast<float>(numProcessedEdges_)),
                        countNumVerticesInForwardTree(), countNumVerticesInReverseTree());
        }

        void APTstar::informAboutPlannerStatus(ompl::base::PlannerStatus::StatusType status) const
        {
            switch (status)
            {
                case ompl::base::PlannerStatus::StatusType::EXACT_SOLUTION:
                {
                    OMPL_INFORM("%s (%u iterations): Found an exact solution of cost %.4f.", name_.c_str(), iteration_,
                                solutionCost_.value());
                    const char* Afilename_cost_apt = "/tmp/result_cost_apt.csv";
                    struct stat sb4;
                    std::ofstream myFile_cost_apt;

                      // Check if file for COST was created

                    if (stat(Afilename_cost_apt, &sb4) == 0) {
                        OMPL_INFORM("File exists and load new COST - APT*");
                        myFile_cost_apt.open("/tmp/result_cost_apt.csv", std::ios::out | std::ios::app | std::ios::binary); 
                        myFile_cost_apt << solutionCost_.value() << ",";
                    }
                    else {
                        OMPL_INFORM("File does not exisit ->create new file - APT*");
                        myFile_cost_apt.open("/tmp/result_cost_apt.csv");
                        myFile_cost_apt << solutionCost_.value() << ",";
                    }

                    myFile_cost_apt.close();
                    break;
                }
                case ompl::base::PlannerStatus::StatusType::APPROXIMATE_SOLUTION:
                {
                    OMPL_INFORM("%s (%u iterations): Did not find an exact solution, but found an approximate solution "
                                "of cost %.4f which is %.4f away from a goal (in cost space).",
                                name_.c_str(), iteration_, approximateSolutionCost_.value(),
                                approximateSolutionCostToGoal_.value());
                    break;
                }
                case ompl::base::PlannerStatus::StatusType::TIMEOUT:
                {
                    if (trackApproximateSolutions_)
                    {
                        OMPL_INFORM("%s (%u iterations): Did not find any solution.", name_.c_str(), iteration_);
                    }
                    else
                    {
                        OMPL_INFORM("%s (%u iterations): Did not find an exact solution, and tracking approximate "
                                    "solutions is disabled.",
                                    name_.c_str(), iteration_);
                    }
                    break;
                }
                case ompl::base::PlannerStatus::StatusType::UNKNOWN:
                case ompl::base::PlannerStatus::StatusType::INVALID_START:
                case ompl::base::PlannerStatus::StatusType::INVALID_GOAL:
                case ompl::base::PlannerStatus::StatusType::UNRECOGNIZED_GOAL_TYPE:
                case ompl::base::PlannerStatus::StatusType::CRASH:
                case ompl::base::PlannerStatus::StatusType::ABORT:
                case ompl::base::PlannerStatus::StatusType::TYPE_COUNT:
                {
                    OMPL_INFORM("%s (%u iterations): Unable to solve the given planning problem.", name_.c_str(),
                                iteration_);
                }
            }
        }

        unsigned int APTstar::countNumVerticesInForwardTree() const
        {
            const auto states = graph_.getStates();
            return std::count_if(states.cbegin(), states.cend(),
                                 [](const auto &state) { return state->hasForwardVertex(); });
        }

        unsigned int APTstar::countNumVerticesInReverseTree() const
        {
            const auto states = graph_.getStates();
            return std::count_if(states.cbegin(), states.cend(),
                                 [](const auto &state) { return state->hasReverseVertex(); });
        }

        bool APTstar::couldImproveForwardPath(const Edge &edge) const
        {
            // If we currently don't have a solution, the anwer is yes.
            if (!objective_->isFinite(solutionCost_))
            {
                // Compare the costs of the full path heuristic with the current cost of the start state.
                const auto heuristicPathCost =
                    combine(edge.source->getCurrentCostToCome(),
                            objective_->motionCostHeuristic(edge.source->raw(), edge.target->raw()),
                            objective_->costToGo(edge.target->raw(), problem_->getGoal().get()));
                if (isBetter(heuristicPathCost, solutionCost_))
                {
                    return true;
                }
                else
                {
                    return false;
                }
            }
            else
            {
                return true;
            }
        }

        bool APTstar::couldImproveForwardTree(const Edge &edge) const
        {
            const auto heuristicCostToCome =
                combine(edge.source->getCurrentCostToCome(),
                        objective_->motionCostHeuristic(edge.source->raw(), edge.target->raw()));
            return isBetter(heuristicCostToCome, edge.target->getCurrentCostToCome());
        }

        bool APTstar::doesImproveForwardPath(const Edge &edge, const ompl::base::Cost &edgeCost) const
        {
            // If we don't have a solution yet, the answer is imediately true.
            if (!objective_->isFinite(solutionCost_))
            {
                return true;
            }

            // Check whether it can improve the current solution.
            return isBetter(
                combine(edge.source->getCurrentCostToCome(), edgeCost, edge.target->getLowerBoundCostToGo()),
                solutionCost_);
        }

        bool APTstar::doesImproveForwardTree(const Edge &edge, const ompl::base::Cost &edgeCost) const
        {
            return isBetter(combine(edge.source->getCurrentCostToCome(), edgeCost),
                            edge.target->getCurrentCostToCome());
        }

        ompl::base::Cost APTstar::estimateCostToTarget(const aptstar::Edge &edge) const
        {
            return combine(edge.source->getEstimatedCostToGo(),
                           objective_->motionCostBestEstimate(edge.source->raw(), edge.target->raw()));
        }

        unsigned int APTstar::estimateEffortToTarget(const aptstar::Edge &edge) const
        {
            // if we previously validated (=whitelisted) an edge, the effort to
            // check is zero
            std::size_t checksToCome = 0u;
            if (edge.source->isWhitelisted(edge.target))
            {
                checksToCome = 0u;
            }
            else
            {
                // Get the segment count for the full resolution.
                const std::size_t fullSegmentCount = space_->validSegmentCount(edge.source->raw(), edge.target->raw());
                // Get the number of checks already performed on this edge.
                const std::size_t performedChecks = edge.target->getIncomingCollisionCheckResolution(edge.source);

                checksToCome = fullSegmentCount - performedChecks;
            }

            return edge.source->getEstimatedEffortToGo() + checksToCome;
        }

        bool APTstar::isValid(const Edge &edge) const
        {
            // The number of checks required to determine whether the edge is valid is the valid segment count minus one
            // because we know that the source and target states are valid.
            const std::size_t numChecks = space_->validSegmentCount(edge.source->raw(), edge.target->raw()) - 1u;

            return isValidAtResolution(edge, numChecks);
        }

        bool APTstar::couldBeValid(const Edge &edge) const
        {
            return isValidAtResolution(edge, numSparseCollisionChecksCurrentLevel_);
        }

        bool APTstar::isValidAtResolution(const Edge &edge, std::size_t numChecks) const
        {
            // Check if the edge is whitelisted.
            if (edge.source->isWhitelisted(edge.target))
            {
                return true;
            }

            // If the edge is blacklisted.
            if (edge.source->isBlacklisted(edge.target))
            {
                return false;
            }

            // Get the segment count for the full resolution.
            const std::size_t fullSegmentCount = space_->validSegmentCount(edge.source->raw(), edge.target->raw());

            // The segment count is the number of checks on this level plus 1, capped by the full resolution segment
            // count.
            const auto segmentCount = std::min(numChecks + 1u, fullSegmentCount);

            /***
               Let's say we want to perform seven collision checks on an edge:

               position of checks: |--------x--------x--------x--------x--------x--------x--------x--------|
               indices of checks:           1        2        3        4        5        6        7
               order of testing:            4        2        5        1        6        3        7

               We create a queue that holds segments and always test the midpoint of the segments. We start
               with the outermost indices and then break the segment in half until the segment collapses to a single
               point:

                 1. indices = { (1, 7) }
                    current = (1, 7) -> test midpoint = 4 -> add (1, 3) and (5, 7) to queue

                 2. indices = { (1, 3), (5, 7) }
                    current (1, 3) -> test midpoint = 2 -> add (1, 1) and (3, 3) to queue

                 3. indices = { (5, 7), (1, 1), (3, 3) }
                    current (5, 7) -> test midpoint = 6 -> add (5, 5) and (7, 7) to queue

                 4. indices = { (1, 1), (3, 3), (5, 5), (7, 7) }
                    current (1, 1) -> test midpoint = 1 -> add nothing to the queue

                 5. indices = { (3, 3), (5, 5), (7, 7) }
                    current (3, 3) -> test midpoint = 3 -> add nothing to the queue

                 6. indices = { (5, 5) (7, 7) }
                    current (5, 5) -> test midpoint = 5 -> add nothing to the queue

                 7. indices = { (7, 7) }
                    current (7, 7) -> test midpoint = 7 -> add nothing to the queue
            ***/

            // Store the current check number.
            std::size_t currentCheck = 1u;

            // Get the number of checks already performed on this edge.
            const std::size_t performedChecks = edge.target->getIncomingCollisionCheckResolution(edge.source);

            // Initialize the queue of positions to be tested.
            std::queue<std::pair<std::size_t, std::size_t>> indices;
            indices.emplace(1u, numChecks);

            // Test states while there are states to be tested.
            while (!indices.empty())
            {
                // Get the current segment.
                const auto current = indices.front();

                // Get the midpoint of the segment.
                auto mid = (current.first + current.second) / 2;

                // Only do the detection if we haven't tested this state on a previous level.
                if (currentCheck > performedChecks)
                {
                    space_->interpolate(edge.source->raw(), edge.target->raw(),
                                        static_cast<double>(mid) / static_cast<double>(segmentCount), detectionState_);

                    if (!spaceInfo_->isValid(detectionState_))
                    {
                        // Blacklist the edge.
                        edge.source->blacklist(edge.target);
                        edge.target->blacklist(edge.source);

                        // Register it with the graph.
                        graph_.registerInvalidEdge(edge);
                        return false;
                    }
                }

                // Remove the current segment from the queue.
                indices.pop();

                // Create the first and second half of the split segment if necessary.
                if (current.first < mid)
                {
                    indices.emplace(current.first, mid - 1u);
                }
                if (current.second > mid)
                {
                    indices.emplace(mid + 1u, current.second);
                }

                // Increase the current check number.
                ++currentCheck;
            }

            // Remember at what resolution this edge was already checked. We're assuming that the number of collision
            // checks is symmetric for each edge.
            edge.source->setIncomingCollisionCheckResolution(edge.target, currentCheck - 1u);
            edge.target->setIncomingCollisionCheckResolution(edge.source, currentCheck - 1u);

            // Whitelist this edge if it was checked at full resolution.
            if (segmentCount == fullSegmentCount)
            {
                ++numCollisionCheckedEdges_;
                edge.source->whitelist(edge.target);
                edge.target->whitelist(edge.source);

                graph_.registerWhitelistedState(edge.source);
                graph_.registerWhitelistedState(edge.target);
            }

            return true;
        }

        bool APTstar::isBetter(const ompl::base::Cost &lhs, const ompl::base::Cost &rhs) const
        {
            return objective_->isCostBetterThan(lhs, rhs);
        }

        ompl::base::Cost APTstar::combine(const ompl::base::Cost &lhs, const ompl::base::Cost &rhs) const
        {
            return objective_->combineCosts(lhs, rhs);
        }

        void APTstar::expandStartVerticesIntoForwardQueue()
        {
            for (auto &vertex : startVertices_)
            {
                forwardQueue_->insertOrUpdate(expand(vertex->getState(), true));
            }
        }

        void APTstar::expandGoalVerticesIntoReverseQueue()
        {
            for (auto &vertex : goalVertices_)
            {
                reverseQueue_->insertOrUpdate(expand(vertex->getState(), false));
            }
        }

        bool APTstar::isClosed(const std::shared_ptr<Vertex> &vertex) const
        {
            return vertex->getExpandTag() == reverseSearchTag_;
        }

        bool APTstar::isInForwardTree(const Edge &edge) const
        {
            if (!edge.source->hasForwardVertex() || !edge.target->hasForwardVertex())
            {
                return false;
            }

            return edge.target->asForwardVertex()->isParent(edge.source->asForwardVertex());
        }

        bool APTstar::isInReverseTree(const Edge &edge) const
        {
            if (!edge.source->hasReverseVertex() || !edge.target->hasReverseVertex())
            {
                return false;
            }

            return edge.target->asReverseVertex()->isParent(edge.source->asReverseVertex());
        }

        bool APTstar::doesImproveReversePath(const Edge &edge) const
        {
            // If there is no reverse path the answer is ja.
            if (!objective_->isFinite(reverseCost_))
            {
                return true;
            }

            // Compare the costs of the full path heuristic with the current cost of the start state.
            const auto heuristicPathCost =
                combine(edge.source->getAdmissibleCostToGo(),
                        objective_->motionCostHeuristic(edge.source->raw(), edge.target->raw()),
                        edge.target->getLowerBoundCostToCome());

            return isBetter(heuristicPathCost, reverseCost_);
        }

        bool APTstar::doesImproveReverseTree(const Edge &edge, const ompl::base::Cost &admissibleEdgeCost) const
        {
            return isBetter(combine(edge.source->getAdmissibleCostToGo(), admissibleEdgeCost),
                            edge.target->getAdmissibleCostToGo());
        }

        std::vector<Edge> APTstar::expandUnlessGoal(const std::shared_ptr<State> &state) const
        {
            if (graph_.isGoal(state))
            {
                return {};
            }

            return expand(state, true);
        }

        std::vector<Edge> APTstar::expand(const std::shared_ptr<State> &state, bool iterateForwardSearch) const
        {
            // Only states associated with a vertex in either of the trees should be expanded.
            assert(state->hasForwardVertex() || state->hasReverseVertex());
            // Prepare the return variable.
            std::vector<Edge> outgoingEdges;

            // Get the neighbors in the current graph.
            if (iterateForwardSearch)
            {
                for (auto &neighborState : graph_.getAllNeighbors(state, graph_.getStartStates()[0],
                                                                  graph_.getGoalStates()[0], iterateForwardSearch))
                {
                    // std::vector<std::shared_ptr<State>> Allneighbors;
                    // Allneighbors.emplace_back(neighborState.lock());
                    // totalForce_.totalForce(state, Allneighbors);
                    // totalForce_.totalForcewithStart(state, graph_.getStartStates()[0], graph_.getGoalStates()[0],
                    // true);
                    double radius = graph_.getRadius();
                    //if (graph_.isInEllipticalRange(state, neighborState, graph_.totalforceDirection_, radius))
                        outgoingEdges.emplace_back(state, neighborState.lock());
                }

                // double radius = graph_.getRadius();
                // std::cout << "forward radius: "<< radius << std::endl;
                // for (int i = 0; i < graph_.totalforceDirection_.size(); i++)

                //     std::cout << "forward forceDirection  " << i << ": "<< graph_.totalforceDirection_[i] << std::endl;
            }
            else
            {
                for (auto &neighborState : graph_.getAllNeighbors(state, graph_.getStartStates()[0],
                                                                  graph_.getGoalStates()[0], iterateForwardSearch))
                {
                    double radius = graph_.getRadius();

                    // std::cout << radius << std::endl;
                    //if (graph_.isInEllipticalRange(state, neighborState, graph_.totalforceDirection_, radius))
                        outgoingEdges.emplace_back(state, neighborState.lock());
                }
                // double radius = graph_.getRadius();
                // std::cout << "reverse radius: "<< radius << std::endl;
                
                // for (int i = 0; i < graph_.totalforceDirection_.size(); i++)

                //     std::cout << "reverse graph_.totalforceDirection_:: " << graph_.totalforceDirection_[i] << std::endl;
            }

            // If the state is in the forward search tree, extra edges have to be added.
            if (state->hasForwardVertex())
            {
                // Get the vertex in the forward search tree associated with this state.
                auto forwardVertex = state->asForwardVertex();

                // Add the outgoing edge to this vertex's parent in the forward tree, if it exists.
                if (!graph_.isStart(state))
                {
                    // If this vertex is not a start state, it must have a parent.
                    assert(forwardVertex->getParent().lock());

                    // Get the state associated with the parent vertex.
                    auto forwardParentState = forwardVertex->getParent().lock()->getState();

                    // Add the edge to the forward tree parent if it has not already being added.
                    if (std::find_if(outgoingEdges.cbegin(), outgoingEdges.cend(),
                                     [&forwardParentState](const auto &edge) {
                                         return edge.target->getId() == forwardParentState->getId();
                                     }) == outgoingEdges.cend())
                    {
                        outgoingEdges.emplace_back(state, forwardParentState);
                    }
                }

                // Add the edge to the forward children.
                for (const auto &child : forwardVertex->getChildren())
                {
                    // Get the state associated with the child vertex.
                    auto forwardChildState = child->getState();

                    // Add the edge to the forward tree child if it has not already being added.
                    if (std::find_if(outgoingEdges.cbegin(), outgoingEdges.cend(),
                                     [&forwardChildState](const auto &edge) {
                                         return edge.target->getId() == forwardChildState->getId();
                                     }) == outgoingEdges.cend())
                    {
                        outgoingEdges.emplace_back(state, forwardChildState);
                    }
                }
            }

            // If the state is in the reverse search tree, extra edges have to be added.
            if (state->hasReverseVertex())
            {
                // Get the vertex in the reverse search tree associated with this state.
                auto reverseVertex = state->asReverseVertex();

                // Add the outgoing edge to this vertex's parent in the reverse tree, if it exists.
                if (!graph_.isGoal(state))
                {
                    // If this state is not a goal, it must have a parent.
                    assert(reverseVertex->getParent().lock());

                    // Get the state associated with the parent vertex.
                    auto reverseParentState = reverseVertex->getParent().lock()->getState();

                    // Add the edge to the reverse tree parent if it has not already being added.
                    if (std::find_if(outgoingEdges.cbegin(), outgoingEdges.cend(),
                                     [&reverseParentState](const auto &edge) {
                                         return edge.target->getId() == reverseParentState->getId();
                                     }) == outgoingEdges.cend())
                    {
                        outgoingEdges.emplace_back(state, reverseParentState);
                    }
                }

                // Add the edge to the reverse children.
                for (const auto &child : reverseVertex->getChildren())
                {
                    // Get the state associated with the child vertex.
                    auto reverseChildState = child->getState();

                    // Add the edge to the reverse tree parent if it has not already being added.
                    if (std::find_if(outgoingEdges.cbegin(), outgoingEdges.cend(),
                                     [&reverseChildState](const auto &edge) {
                                         return edge.target->getId() == reverseChildState->getId();
                                     }) == outgoingEdges.cend())
                    {
                        outgoingEdges.emplace_back(state, reverseChildState);
                    }
                }
            }

            return outgoingEdges;
        }
    }  // namespace geometric
}  // namespace ompl
