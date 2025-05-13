/* Copyright 2021 The MathWorks, Inc. */

#ifndef COLUMN_MAJOR_ITERATOR_HPP_
#define COLUMN_MAJOR_ITERATOR_HPP_

#include "detail/IteratorFactory.hpp"
#include "detail/OrderedIterator.hpp"

namespace matlab {
namespace data {

template <typename T>
using ColumnMajorIterator =
    detail::OrderedIterator<T, detail::ColumnMajorOrder>;

using ColumnMajor = detail::IteratorFactory<ColumnMajorIterator>;

} // namespace data
} // namespace matlab

#endif
