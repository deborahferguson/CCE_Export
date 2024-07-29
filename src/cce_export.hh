#include <string>
#include <vector>

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"
#include "util_Table.h"

using namespace std;

extern "C" void CCE_Export(CCTK_ARGUMENTS);

void CCE_Export_Interpolate_On_Sphere(CCTK_ARGUMENTS, vector<CCTK_REAL> &xs, vector<CCTK_REAL> &ys, vector<CCTK_REAL> &zs, std::string name, vector<CCTK_REAL> &sphere_values, vector<CCTK_REAL> &sphere_dx, vector<CCTK_REAL> &sphere_dy, vector<CCTK_REAL> &sphere_dz, CCTK_INT array_size);
