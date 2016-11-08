#ifndef NONSTD_STAT_H_
#define NONSTD_STAT_H_

#include <cmath>
#include <cstddef>
#include <string>

#include "nonstd/numerics.h"
#include "nonstd/runtime.h"

namespace nonstd {

  ///@name Statistics
  ///@{

  /// Compute the entropy of a container of outcomes. Each element in the
  /// must be either a frequency or a number of occurences of a given category.
  template <typename container_type>
  double
  entropy(const container_type& _c)
  {
    double sum = 0;
    for(const auto& elem : _c)
      sum += elem;
    double ent = 0;
    for(const auto& elem : _c)
      ent -= (elem == 0) ? 0 : elem * log2(static_cast<double>(elem)/sum);
    ent /= sum;
    return ent;
  }


  /// Compute the number of distinct ways to choose _k elements from a set of
  /// size _n where order is not relevant.
  inline size_t
  combinations(size_t _n, size_t _k) noexcept
  {
    assert_msg(_n >= _k, "error: called combinations with _k = " +
        std::to_string(_k) + " > _n = " + std::to_string(_n) + "!");
    return factorial(_n) / (factorial(_k) * factorial(_n - _k));
  }


  /// Compute the number of distinct ways to choose _k elements from a set of
  /// size _n where order is relevant.
  inline size_t
  permutations(size_t _n, size_t _k) noexcept
  {
    assert_msg(_n >= _k, "error: called permutations with _k = " +
        std::to_string(_k) + " > _n = " + std::to_string(_n) + "!");
    return factorial(_n) / factorial(_n - _k);
  }


  //////////////////////////////////////////////////////////////////////////////
  /// Object representation of a binomial distribution.
  //////////////////////////////////////////////////////////////////////////////
  class binomial_distribution final
  {

    ///@name Internal State
    ///@{

    size_t n; ///< The number of trials.
    double p; ///< The probability of success.

    ///@}

    public:

      ///@name Construction
      ///@{

      binomial_distribution(size_t _n, double _p);
      ~binomial_distribution() = default;

      ///@}
      ///@name Interface
      ///@{

      size_t num_trials() const noexcept;
      double probability() const noexcept;
      double probability(size_t _r) const noexcept;
      double mean() const noexcept;
      double variance() const noexcept;
      double std_dev() const noexcept;

      ///@}

  };

  ///@}
}

#endif
