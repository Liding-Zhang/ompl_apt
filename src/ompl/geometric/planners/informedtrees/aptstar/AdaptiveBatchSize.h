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

// Authors: Liding Zhang, Xu Liang

#ifndef OMPL_GEOMETRIC_PLANNERS_INFORMEDTREES_APTSTAR_ADAPTIVEBATCHSIZE_
#define OMPL_GEOMETRIC_PLANNERS_INFORMEDTREES_APTSTAR_ADAPTIVEBATCHSIZE_

#include "ompl/geometric/planners/informedtrees/aptstar/RandomGeometricGraph.h"
#include "ompl/geometric/planners/informedtrees/APTstar.h"

namespace ompl
{
    namespace geometric
    {
        namespace aptstar
        {
            enum DecayMethod
            {
                ITERATION_TIMES,
                LINEAR,
                PARABOLA,
                LOG,
                BRACHISTOCHRONE,
            };

            class AdaptiveBatchSize
            {
            public:
                AdaptiveBatchSize(const DecayMethod &decay_method, ompl::base::Cost &solutionCost,
                                  const double &minPossibleCost, unsigned int &batchSize, double &S_max_initial,
                                  double &S_min_initial, const unsigned int &maxSamples,
                                  const unsigned int &minSamples, std::size_t dim);
                ~AdaptiveBatchSize() = default;

                unsigned int adjustBatchSize(DecayMethod decay_method);

            private:
                /** \brief Type of Decay method */
                unsigned int adjustBatchSizeIterationTimes();

                /** \brief Type of Decay method */
                unsigned int adjustBatchSizeLinear();

                /** \brief Type of Decay method */
                unsigned int adjustBatchSizeParabola();

                /** \brief Type of Decay method */
                unsigned int adjustBatchSizeLog();

                /** \brief Type of Decay method */
                unsigned int adjustBatchSizeBrachistochrone();

                /** \brief Type of Decay method */
                DecayMethod decay_method_;

                /** \brief Initialize the solution cost. */
                ompl::base::Cost _solutionCost_{std::numeric_limits<double>::infinity()};

                /** \brief Initialize the minial possible cost. */
                double minPossibleCost_;

                /** \brief Initialize the batch size. */
                unsigned int _batchSize_;

                /** \brief Initialize the maximal area of n-dimensional ellipsoid. */
                double _S_max_initial_;

                /** \brief Initialize the minimal area of n-dimensional ellipsoid. */
                double _S_min_initial_;

                /** \brief Initialize the maximal number per sample. */
                unsigned int _maxSamples_;

                /** \brief Initialize the minimal number per sample. */
                unsigned int _minSamples_;

                /** \brief Initialize the iteration number. */
                std::size_t _iteration_;

                /** \brief dimension of the states. */
                std::size_t dim_;
            };

        }  // namespace aptstar

    }  // namespace geometric

}  // namespace ompl

#endif  // OMPL_GEOMETRIC_PLANNERS_INFORMEDTREES_APTSTAR_ADAPTIVEBATCHSIZE_