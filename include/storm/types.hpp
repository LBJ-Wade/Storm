#ifndef STORM_TYPES_HPP
#define STORM_TYPES_HPP

#include <Pythia8/Pythia.h>
#include <boost/histogram.hpp>
#include <boost/histogram/axis/regular.hpp>
#include <boost/histogram/weight.hpp>

namespace storm {

using Pythia8::Vec4;

using MomentaList = std::vector<Pythia8::Vec4>;
using Spectum = std::pair<std::vector<double>, std::vector<double>>;

using SquaredMatrixElement =
    std::function<double(const std::vector<Pythia8::Vec4> &)>;

namespace bp = boost::histogram;

using LogAxis = bp::axis::regular<double, bp::axis::transform::log>;
using Histogram = bp::histogram<std::tuple<LogAxis>>;

enum Generation : size_t {
  first = 0,
  second = 1,
  third = 2,
};

} // namespace storm

#endif // STORM_TYPES_HPP
