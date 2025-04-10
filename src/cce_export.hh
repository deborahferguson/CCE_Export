#ifndef CCE_EXPORT_HH
#define CCE_EXPORT_HH

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"
#include "util_Table.h"

#include <string>
#include <vector>

using std::vector, std::string;

namespace CCE_export {

extern "C" void CCE_Export(CCTK_ARGUMENTS);

CCTK_REAL Simpson2DIntegral(CCTK_REAL const *f, int nx, int ny, CCTK_REAL hx,
                            CCTK_REAL hy);

void Integrate(int array_size, int ntheta, int nphi, vector<CCTK_REAL> &array1r,
               vector<CCTK_REAL> &array1i, vector<CCTK_REAL> &array2r,
               vector<CCTK_REAL> &th, vector<CCTK_REAL> &ph, CCTK_REAL *outre,
               CCTK_REAL *outim);

void Decompose_Spherical_Harmonics(
    vector<CCTK_REAL> &th, vector<CCTK_REAL> &phi,
    vector<CCTK_REAL> &sphere_values, vector<CCTK_REAL> &re_data,
    vector<CCTK_REAL> &im_data, vector<vector<CCTK_REAL> > &re_ylms,
    vector<vector<CCTK_REAL> > &im_ylms, int array_size, int lmax, int ntheta,
    int nphi);

void Compute_Ylms(vector<CCTK_REAL> &th, vector<CCTK_REAL> &ph,
                  vector<vector<CCTK_REAL> > &re_ylms,
                  vector<vector<CCTK_REAL> > &im_ylms, int lmax,
                  int array_size);

} // namespace CCE_export

#endif // CCE_EXPORT_HH
