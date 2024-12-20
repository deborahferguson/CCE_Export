#include <string>
#include <vector>

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"
#include "util_Table.h"

using namespace std;

extern "C" void CCE_Export(CCTK_ARGUMENTS);

void CCE_Export_Interpolate_On_Sphere_With_Derivatives(CCTK_ARGUMENTS, vector<CCTK_REAL> &xs, vector<CCTK_REAL> &ys, vector<CCTK_REAL> &zs, std::string name, vector<CCTK_REAL> &sphere_values, vector<CCTK_REAL> &sphere_dx, vector<CCTK_REAL> &sphere_dy, vector<CCTK_REAL> &sphere_dz, CCTK_INT array_size);

void CCE_Export_Interpolate_On_Sphere(CCTK_ARGUMENTS, vector<CCTK_REAL> &xs, vector<CCTK_REAL> &ys, vector<CCTK_REAL> &zs, std::string name, vector<CCTK_REAL> &sphere_values, CCTK_INT array_size);

void Decompose_Spherical_Harmonics(vector<CCTK_REAL> &th, vector<CCTK_REAL> &phi, vector<CCTK_REAL> &sphere_values, vector<CCTK_REAL> &re_data, vector<CCTK_REAL> &im_data, int array_size);

int l_m_to_index(int l, int m);

int factorial(int x);

int Binomial_Coefficient(int n, int k);

CCTK_REAL Legendre_Polynomial(int l, int m, CCTK_REAL x);

void Compute_Ylms(vector<CCTK_REAL> &th, vector<CCTK_REAL> &ph, vector<vector<CCTK_REAL>> &re_ylms, vector<vector<CCTK_REAL>> &im_ylms, int lmax, int array_size);
