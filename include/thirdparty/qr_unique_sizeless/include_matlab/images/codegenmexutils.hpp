/**
 * @file         codegenmexutils.hpp
 *
 * Purpose       Declaration of utility functions for code generation entry
 * point support
 *
 * Authors       Dinesh Iyer
 *
 * Copyright     2019-2020 The MathWorks, Inc.
 *
 */

#ifndef imagesdatatypescodegenmexutils_header
#define imagesdatatypescodegenmexutils_header

#include "matrix.h"

#ifdef BUILDING_IMAGESDATATYPESCODEGENMEXUTILS
#include "images/datatypes/Image.hpp"
#include "images/datatypes/codegenc/codegenc_api.hpp"
#include "images/datatypes/codegenmexutils/codegenmexutils_util.hpp"
#else
#include "Image.hpp"
#include "codegenc_api.hpp"
#include "codegenmexutils_util.hpp" /* To deal with the fact that PackNGo has no include file hierarchy during test */
#endif

// ------------------------------------------------
namespace images {
namespace datatypes {
namespace codegenmexutils {
/**
 * Conversion functions for C++ codegen
 */
IMAGESDATATYPESCODEGENMEXUTILS_API mxArray *
ImageToMxArray(const images::datatypes::Image &src);

IMAGESDATATYPESCODEGENMEXUTILS_API images::datatypes::Image
MxArrayToImage(const mxArray *im);
} // namespace codegenmexutils
} // namespace datatypes
} // namespace images

#ifdef __cplusplus
extern "C" {
#endif
EXTERN_C IMAGESDATATYPESCODEGENMEXUTILS_API mxArray *
ImagePtrToMxArray(const ImagePtr src);

EXTERN_C IMAGESDATATYPESCODEGENMEXUTILS_API ImagePtr
mxArrayToImagePtr(const mxArray *im);

#ifdef __cplusplus
}
#endif

#endif // imagesdatatypescodegenmexutils_header
