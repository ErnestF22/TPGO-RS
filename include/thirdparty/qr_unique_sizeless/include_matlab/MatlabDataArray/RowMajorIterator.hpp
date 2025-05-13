/* Copyright 2021 The MathWorks, Inc. */

#ifndef ROW_MAJOR_ITERATOR_HPP_
#define ROW_MAJOR_ITERATOR_HPP_

#include "detail/IteratorFactory.hpp"
#include "detail/OrderedIterator.hpp"

namespace matlab {
namespace data {

template <typename T>
using RowMajorIterator = detail::OrderedIterator<T, detail::RowMajorOrder>;

using RowMajor = detail::IteratorFactory<RowMajorIterator>;

} // namespace data
} // namespace matlab

#endif
