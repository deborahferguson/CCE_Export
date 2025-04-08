#include <string>
#include <vector>

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"
#include "util_Table.h"

using namespace std;

extern "C" void CCE_Export(CCTK_ARGUMENTS);

void CCE_Export_Interpolate_On_Sphere_With_Derivatives(
    CCTK_ARGUMENTS, vector<CCTK_REAL> &xs, vector<CCTK_REAL> &ys,
    vector<CCTK_REAL> &zs, std::string name, vector<CCTK_REAL> &sphere_values,
    vector<CCTK_REAL> &sphere_dx, vector<CCTK_REAL> &sphere_dy,
    vector<CCTK_REAL> &sphere_dz, CCTK_INT array_size);

void CCE_Export_Interpolate_On_Sphere(CCTK_ARGUMENTS, vector<CCTK_REAL> &xs,
                                      vector<CCTK_REAL> &ys,
                                      vector<CCTK_REAL> &zs, std::string name,
                                      vector<CCTK_REAL> &sphere_values,
                                      CCTK_INT array_size);

void Decompose_Spherical_Harmonics(
    vector<CCTK_REAL> &th, vector<CCTK_REAL> &phi,
    vector<CCTK_REAL> &sphere_values, vector<CCTK_REAL> &re_data,
    vector<CCTK_REAL> &im_data, vector<vector<CCTK_REAL> > &re_ylms,
    vector<vector<CCTK_REAL> > &im_ylms, int array_size, int lmax, int ntheta,
    int nphi);

int l_m_to_index(int l, int m);

CCTK_REAL factorial(CCTK_REAL x);

CCTK_REAL Binomial_Coefficient(CCTK_REAL n, CCTK_REAL k);

CCTK_REAL Legendre_Polynomial(int l, int m, CCTK_REAL x);

void Compute_Ylms(vector<CCTK_REAL> &th, vector<CCTK_REAL> &ph,
                  vector<vector<CCTK_REAL> > &re_ylms,
                  vector<vector<CCTK_REAL> > &im_ylms, int lmax,
                  int array_size);

CCTK_REAL CCE_Export_Simpson2DIntegral(CCTK_REAL const *f, int nx, int ny,
                                       CCTK_REAL hx, CCTK_REAL hy);

void Create_Dataset(string datasetname, CCTK_REAL *data, int mode_count);

void Output_Decomposed_Metric_Data(
    CCTK_ARGUMENTS, vector<vector<vector<CCTK_REAL> > > &re_g,
    vector<vector<vector<CCTK_REAL> > > &im_g,
    vector<vector<vector<CCTK_REAL> > > &re_dr_g,
    vector<vector<vector<CCTK_REAL> > > &im_dr_g,
    vector<vector<vector<CCTK_REAL> > > &re_dt_g,
    vector<vector<vector<CCTK_REAL> > > &im_dt_g,
    vector<vector<CCTK_REAL> > &re_beta, vector<vector<CCTK_REAL> > &im_beta,
    vector<vector<CCTK_REAL> > &re_dr_beta,
    vector<vector<CCTK_REAL> > &im_dr_beta,
    vector<vector<CCTK_REAL> > &re_dt_beta,
    vector<vector<CCTK_REAL> > &im_dt_beta, vector<CCTK_REAL> &re_alpha,
    vector<CCTK_REAL> &im_alpha, vector<CCTK_REAL> &re_dr_alpha,
    vector<CCTK_REAL> &im_dr_alpha, vector<CCTK_REAL> &re_dt_alpha,
    vector<CCTK_REAL> &im_dt_alpha, float radius, int lmax);

void CCE_Export_Integrate(int array_size, int ntheta, int nphi,
                          vector<CCTK_REAL> &array1r,
                          vector<CCTK_REAL> &array1i,
                          vector<CCTK_REAL> &array2r, vector<CCTK_REAL> &th,
                          vector<CCTK_REAL> &ph, CCTK_REAL *outre,
                          CCTK_REAL *outim);
