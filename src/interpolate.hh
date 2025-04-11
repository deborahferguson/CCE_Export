#ifndef INTERPOLATE_HH
#define INTERPOLATE_HH

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"
#include "util_Table.h"

#include "vector"
#include "string"

using std::vector, std::string;

namespace CCE_export {

void Interpolate_On_Sphere_With_Derivatives(
    CCTK_ARGUMENTS, vector<CCTK_REAL> &xs, vector<CCTK_REAL> &ys,
    vector<CCTK_REAL> &zs, std::string name, vector<CCTK_REAL> &sphere_values,
    vector<CCTK_REAL> &sphere_dx, vector<CCTK_REAL> &sphere_dy,
    vector<CCTK_REAL> &sphere_dz, CCTK_INT array_size);

void Interpolate_On_Sphere(CCTK_ARGUMENTS, vector<CCTK_REAL> &xs,
                           vector<CCTK_REAL> &ys, vector<CCTK_REAL> &zs,
                           std::string name, vector<CCTK_REAL> &sphere_values,
                           CCTK_INT array_size);

} // namespace CCE_export

#endif
