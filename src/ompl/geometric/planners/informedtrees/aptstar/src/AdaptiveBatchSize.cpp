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

#include "ompl/geometric/planners/informedtrees/aptstar/AdaptiveBatchSize.h"

namespace ompl
{
    namespace geometric
    {
        namespace aptstar
        {
            AdaptiveBatchSize::AdaptiveBatchSize(const DecayMethod &decay_method, ompl::base::Cost &solutionCost,
                                                 const double &minPossibleCost, unsigned int &batchSize,
                                                 double &S_max_initial, double &S_min_initial,
                                                 const unsigned int &maxSamples, const unsigned int &minSamples,
                                                 std::size_t dim)
              : decay_method_(decay_method)
              , _solutionCost_(solutionCost)
              , minPossibleCost_(minPossibleCost)
              , _batchSize_(batchSize)
              , _S_max_initial_(S_max_initial)
              , _S_min_initial_(S_min_initial)
              , _maxSamples_(maxSamples)
              , _minSamples_(minSamples)
              , dim_(dim)
            {
            }
            unsigned int AdaptiveBatchSize::adjustBatchSize(DecayMethod decay_method)
            {
                switch (decay_method)
                {
                    case ITERATION_TIMES:
                        return AdaptiveBatchSize::adjustBatchSizeIterationTimes();
                    case LINEAR:
                        return AdaptiveBatchSize::adjustBatchSizeLinear();
                    case PARABOLA:
                        return AdaptiveBatchSize::adjustBatchSizeParabola();
                    case LOG:
                        return AdaptiveBatchSize::adjustBatchSizeLog();
                    case BRACHISTOCHRONE:
                        return AdaptiveBatchSize::adjustBatchSizeBrachistochrone();
                    default:
                        return _batchSize_;
                }
            }

            unsigned int AdaptiveBatchSize::adjustBatchSizeIterationTimes()
            {
                // Implementation for ITERATION_TIMES method
                if (_iteration_ >= 1000 && _iteration_ < 5000)
                {
                    return _batchSize_ = 150;
                }

                else if (_iteration_ >= 5000 && _iteration_ < 8000)
                {
                    return _batchSize_ = 100;
                }
                else if (_iteration_ >= 8000 && _iteration_ < 10000)
                {
                    return _batchSize_ = 80;
                }
                else
                {
                    return _batchSize_;
                }
            }

            unsigned int AdaptiveBatchSize::adjustBatchSizeLinear()
            {
                // Implementation for LINEAR method
                if (std::isinf(_solutionCost_.value()))
                {
                    return _batchSize_;
                }

                double a = _solutionCost_.value() / 2;
                double c = minPossibleCost_ / 2;
                double b = std::sqrt(a * a - c * c);
                double S = M_PI * a * b;
                static bool pragma = false;
                if (!pragma)
                {
                    _S_max_initial_ = S;
                    pragma = true;
                }

                double ratio = S / _S_max_initial_;

                _batchSize_ = _minSamples_ + (_maxSamples_ - _minSamples_) * ratio;

                // Clamp batchSize_ to be within range
                _batchSize_ = std::max(_minSamples_, std::min(_maxSamples_, _batchSize_));
                return _batchSize_;
            }

            unsigned int AdaptiveBatchSize::adjustBatchSizeParabola()
            {
                // Implementation for PARABOLA method
                if (std::isinf(_solutionCost_.value()))
                {
                    return _batchSize_;
                }

                double a = _solutionCost_.value() / 2;
                double c = minPossibleCost_ / 2;
                double b = std::sqrt(a * a - c * c);
                double S = M_PI * a * b;
                static bool pragma = false;
                if (!pragma)
                {
                    _S_max_initial_ = S;
                    pragma = true;
                }

                double ratio = S / _S_max_initial_;

                // parabola curve
                double parabola_factor = std::sqrt(ratio);

                _batchSize_ = _minSamples_ + (_maxSamples_ - _minSamples_) * parabola_factor;
                // Clamp batchSize_ to be within [minSamples_, maxSamples_]
                _batchSize_ = std::max(_minSamples_, std::min(_maxSamples_, _batchSize_));
                return _batchSize_;
            }

            unsigned int AdaptiveBatchSize::adjustBatchSizeLog()
            {
                // Implementation for LOG & SIGMOID method
                if (std::isinf(_solutionCost_.value()))
                {
                    return _batchSize_;
                }

                double a = _solutionCost_.value() / 2;
                double c = minPossibleCost_ / 2;
                double b = std::sqrt(a * a - c * c);
                double S = M_PI * a * b;
                static bool pragma = false;
                if (!pragma)
                {
                    _S_max_initial_ = S;
                    pragma = true;
                }

                double ratio = S / _S_max_initial_;

                // Logarithmic decay
                double lambda =
                    (_minSamples_ + _maxSamples_) / dim_;  // Adjust this value to make the decay faster or slower
                // Sigmoid function to smooth ratio
                double smoothedValue = 1 / (1 + exp(-10 * (ratio - 0.5)));
                double decay_factor = log(1 + lambda * smoothedValue) / log(1 + lambda);
                _batchSize_ = _minSamples_ + (_maxSamples_ - _minSamples_) * decay_factor;

                // Clamp batchSize_ to be within [minSamples_, maxSamples_]
                _batchSize_ = std::max(_minSamples_, std::min(_maxSamples_, _batchSize_));
                // areaList_[ratio]++;

                // std::cout << "---------------batchsize of current batch is : " << _batchSize_ << std::endl;
                return _batchSize_;
            }

            unsigned int AdaptiveBatchSize::adjustBatchSizeBrachistochrone()
            {
                // Implementation for BRACHISTOCHRONE method
                if (std::isinf(_solutionCost_.value()))
                {
                    return _batchSize_;
                }

                double a = _solutionCost_.value() / 2;
                double c = minPossibleCost_ / 2;
                double b = std::sqrt(a * a - c * c);
                double S = M_PI * a * b;
                static bool pragma = false;
                if (!pragma)
                {
                    _S_max_initial_ = S;
                    pragma = true;
                }

                double ratio = S / _S_max_initial_;

                // Sigmoid function to smooth ratio
                double smoothedValue = 1 / (1 + exp(-10 * (ratio - 0.5)));
                // Brachistochrone curve decay
                double theta = M_PI / 2 * smoothedValue;
                double brachistochrone_factor = std::sqrt(sin(theta));
                _batchSize_ = _minSamples_ + (_maxSamples_ - _minSamples_) * brachistochrone_factor;

                // Clamp batchSize_ to be within [minSamples_, maxSamples_]
                _batchSize_ = std::max(_minSamples_, std::min(_maxSamples_, _batchSize_));
                return _batchSize_;
            }

        }  // namespace aptstar

    }  // namespace geometric

}  // namespace ompl
