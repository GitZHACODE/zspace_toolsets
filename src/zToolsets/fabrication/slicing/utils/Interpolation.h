#pragma once
#include <cmath>
#include <limits>

namespace zSpace { namespace zTsDetail {
// Affine form preserves a constant scalar exactly. Accumulating separately
// rounded weighted float terms can exceed a constant by one float ULP.
// This routine has no fabrication bounds and never clamps the scalar value.
inline double interpolateScalar(double a,double b,double c,double u,double v,double w){
  const double sum=u+v+w;
  if(!std::isfinite(sum)||!(sum>0))return std::numeric_limits<double>::quiet_NaN();
  if(u>=v&&u>=w)return a+(v/sum)*(b-a)+(w/sum)*(c-a);
  if(v>=w)return b+(u/sum)*(a-b)+(w/sum)*(c-b);
  return c+(u/sum)*(a-c)+(v/sum)*(b-c);
}
} }
