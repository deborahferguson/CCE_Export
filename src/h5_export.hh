#ifndef H5_EXPORT_HH
#define H5_EXPORT_HH

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"

#ifdef HAVE_CAPABILITY_HDF5
// We currently support the HDF5 1.6 API (and when using 1.8 the
// compatibility mode introduced by H5_USE_16_API).  Several machines
// in SimFactory use HDF5 1.6, so we cannot drop support for it.  It
// seems it is hard to support both the 1.6 and 1.8 API
// simultaneously; for example H5Fopen takes a different number of
// arguments in the two versions.
#define H5_USE_16_API
#include <hdf5.h>
#endif

#include <string>
#include <vector>

using std::vector, std::string;

namespace CCE_export {

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

} // namespace CCE_export

#endif