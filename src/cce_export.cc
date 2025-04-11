#include "cce_export.hh"
#include "interpolate.hh"
#include "h5_export.hh"
#include <vector>
#include <sys/stat.h>
#include <iomanip>
#include <string.h>
#include <map>

using std::vector, std::string, std::ostringstream, std::map, std::ios,
    std::setprecision;

namespace CCE_export {

void CCE_Export(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  printf("CCE_Export: starting\n");

  static string index_to_component[] = {"x", "y", "z"};

  const int ntheta = 120;
  const int nphi = 240;
  const int array_size = (ntheta + 1) * (nphi + 1);

  // extrinsic curvature, 3d vector (3, 3, array_size)
  vector<vector<vector<CCTK_REAL> > > k(
      3, vector<vector<CCTK_REAL> >(3, vector<CCTK_REAL>(array_size)));
  vector<vector<vector<CCTK_REAL> > > dx_k(
      3, vector<vector<CCTK_REAL> >(3, vector<CCTK_REAL>(array_size)));
  vector<vector<vector<CCTK_REAL> > > dy_k(
      3, vector<vector<CCTK_REAL> >(3, vector<CCTK_REAL>(array_size)));
  vector<vector<vector<CCTK_REAL> > > dz_k(
      3, vector<vector<CCTK_REAL> >(3, vector<CCTK_REAL>(array_size)));
  // metric, 3d vector (3, 3, array_size)
  vector<vector<vector<CCTK_REAL> > > g(
      3, vector<vector<CCTK_REAL> >(3, vector<CCTK_REAL>(array_size)));
  vector<vector<vector<CCTK_REAL> > > dx_g(
      3, vector<vector<CCTK_REAL> >(3, vector<CCTK_REAL>(array_size)));
  vector<vector<vector<CCTK_REAL> > > dy_g(
      3, vector<vector<CCTK_REAL> >(3, vector<CCTK_REAL>(array_size)));
  vector<vector<vector<CCTK_REAL> > > dz_g(
      3, vector<vector<CCTK_REAL> >(3, vector<CCTK_REAL>(array_size)));
  vector<vector<vector<CCTK_REAL> > > dr_g(
      3, vector<vector<CCTK_REAL> >(3, vector<CCTK_REAL>(array_size)));
  vector<vector<vector<CCTK_REAL> > > dt_g(
      3, vector<vector<CCTK_REAL> >(3, vector<CCTK_REAL>(array_size)));
  // shift (beta), 2d vector (3, array_size)
  vector<vector<CCTK_REAL> > beta(3, vector<CCTK_REAL>(array_size));
  vector<vector<CCTK_REAL> > dx_beta(3, vector<CCTK_REAL>(array_size));
  vector<vector<CCTK_REAL> > dy_beta(3, vector<CCTK_REAL>(array_size));
  vector<vector<CCTK_REAL> > dz_beta(3, vector<CCTK_REAL>(array_size));
  vector<vector<CCTK_REAL> > dr_beta(3, vector<CCTK_REAL>(array_size));
  vector<vector<CCTK_REAL> > dt_beta(3, vector<CCTK_REAL>(array_size));
  // lapse (alpha), vector (array_size)
  vector<CCTK_REAL> alpha(array_size);
  vector<CCTK_REAL> dx_alpha(array_size);
  vector<CCTK_REAL> dy_alpha(array_size);
  vector<CCTK_REAL> dz_alpha(array_size);
  vector<CCTK_REAL> dr_alpha(array_size);
  vector<CCTK_REAL> dt_alpha(array_size);

  // x, y, z values of points on the sphere
  vector<CCTK_REAL> xs(array_size);
  vector<CCTK_REAL> ys(array_size);
  vector<CCTK_REAL> zs(array_size);
  // x, y, z unit vectors toward theta and phi points on the sphere
  vector<CCTK_REAL> xhat(array_size);
  vector<CCTK_REAL> yhat(array_size);
  vector<CCTK_REAL> zhat(array_size);
  // theta and phi points on the sphere
  vector<CCTK_REAL> th(array_size);
  vector<CCTK_REAL> ph(array_size);

  // Compute the theta and phi points as well as the corresponding x, y, z unit
  // vectors Based on the number of theta and phi points desired (ntheta, nphi)
  const CCTK_REAL PI = acos(-1.0);

  for (int theta_index = 0; theta_index <= ntheta; theta_index++) {
    for (int phi_index = 0; phi_index <= nphi; phi_index++) {
      const int array_index = theta_index + (ntheta + 1) * phi_index;

      th.at(array_index) = theta_index * PI / (ntheta);
      ph.at(array_index) = phi_index * 2 * PI / nphi;
      xhat.at(array_index) = sin(th.at(array_index)) * cos(ph.at(array_index));
      yhat.at(array_index) = sin(th.at(array_index)) * sin(ph.at(array_index));
      zhat.at(array_index) = cos(th.at(array_index));
    }
  }

  // loop through the desired radii
  for (int r = 0; r < nradii; r++) {

    // compute the values of x, y, z and the desired points on the sphere of
    // radius radius[r]
    for (int theta_index = 0; theta_index <= ntheta; theta_index++) {
      for (int phi_index = 0; phi_index <= nphi; phi_index++) {
        const int array_index = theta_index + (ntheta + 1) * phi_index;

        xs.at(array_index) = radius[r] * xhat.at(array_index);
        ys.at(array_index) = radius[r] * yhat.at(array_index);
        zs.at(array_index) = radius[r] * zhat.at(array_index);
      }
    }

    // Interpolate all the desired quantities onto the desired points on the
    // sphere
    for (int i = 0; i < 3; i++) {
      string first_component = index_to_component[i];
      for (int j = i; j < 3; j++) {
        string second_component = index_to_component[j];
        // interpolate extrinsic curvature
        Interpolate_On_Sphere_With_Derivatives(
            CCTK_PASS_CTOC, xs, ys, zs,
            "ADMBase::k" + first_component + second_component, k.at(i).at(j),
            dx_k.at(i).at(j), dy_k.at(i).at(j), dz_k.at(i).at(j), array_size);
        // interpolate metric
        Interpolate_On_Sphere_With_Derivatives(
            CCTK_PASS_CTOC, xs, ys, zs,
            "ADMBase::g" + first_component + second_component, g.at(i).at(j),
            dx_g.at(i).at(j), dy_g.at(i).at(j), dz_g.at(i).at(j), array_size);
        // compute dr_g
        for (int array_index = 0; array_index < array_size; array_index++) {
          dr_g.at(i).at(j).at(array_index) =
              (xs.at(array_index) / radius[r]) *
                  dx_g.at(i).at(j).at(array_index) +
              (ys.at(array_index) / radius[r]) *
                  dy_g.at(i).at(j).at(array_index) +
              (zs.at(array_index) / radius[r]) *
                  dz_g.at(i).at(j).at(array_index);
        }
      }
      // interpolate shift
      Interpolate_On_Sphere_With_Derivatives(
          CCTK_PASS_CTOC, xs, ys, zs, "ADMBase::beta" + first_component,
          beta.at(i), dx_beta.at(i), dy_beta.at(i), dz_beta.at(i), array_size);
      // interpolate time derivative of shift
      Interpolate_On_Sphere(CCTK_PASS_CTOC, xs, ys, zs,
                            "ADMBase::dtbeta" + first_component, dt_beta.at(i),
                            array_size);
      // compute dr_beta
      for (int array_index = 0; array_index < array_size; array_index++) {
        dr_beta.at(i).at(array_index) =
            (xs.at(array_index) / radius[r]) * dx_beta.at(i).at(array_index) +
            (ys.at(array_index) / radius[r]) * dy_beta.at(i).at(array_index) +
            (zs.at(array_index) / radius[r]) * dz_beta.at(i).at(array_index);
      }
    }
    // interpolate lapse
    Interpolate_On_Sphere_With_Derivatives(CCTK_PASS_CTOC, xs, ys, zs,
                                           "ADMBase::alp", alpha, dx_alpha,
                                           dy_alpha, dz_alpha, array_size);
    // interpolate time derivative of lapse
    Interpolate_On_Sphere(CCTK_PASS_CTOC, xs, ys, zs, "ADMBase::dtalp",
                          dt_alpha, array_size);
    // compute dr_alpha
    for (int array_index = 0; array_index < array_size; array_index++) {
      dr_alpha.at(array_index) =
          (xs.at(array_index) / radius[r]) * dx_alpha.at(array_index) +
          (ys.at(array_index) / radius[r]) * dy_alpha.at(array_index) +
          (zs.at(array_index) / radius[r]) * dz_alpha.at(array_index);
    }

    // compute time derivatives of the metric using the following:
    // d_t g_ij = -2 alpha K_ij
    //            + beta^x dx g_ij + beta^y dy g_ij + beta^z dz g_ij
    //            + g_xi dj beta^x + g_yi dj beta^y + g_zi dj beta^z
    //            + g_xj di beta^x + g_yj di beta^y + g_zj di beta^z

    // define redundant metric components
    for (int array_index = 0; array_index < array_size; array_index++) {
      g.at(1).at(0).at(array_index) = g.at(0).at(1).at(array_index);
      g.at(2).at(0).at(array_index) = g.at(0).at(2).at(array_index);
      g.at(2).at(1).at(array_index) = g.at(1).at(2).at(array_index);
    }

    // dt g_xx
    for (int array_index = 0; array_index < array_size; array_index++) {
      dt_g.at(0).at(0).at(array_index) =
          -2 * alpha.at(array_index) * k.at(0).at(0).at(array_index) +
          beta.at(0).at(array_index) * dx_g.at(0).at(0).at(array_index) +
          beta.at(1).at(array_index) * dy_g.at(0).at(0).at(array_index) +
          beta.at(2).at(array_index) * dz_g.at(0).at(0).at(array_index) +
          g.at(0).at(0).at(array_index) * dx_beta.at(0).at(array_index) +
          g.at(1).at(0).at(array_index) * dx_beta.at(1).at(array_index) +
          g.at(2).at(0).at(array_index) * dx_beta.at(2).at(array_index) +
          g.at(0).at(0).at(array_index) * dx_beta.at(0).at(array_index) +
          g.at(1).at(0).at(array_index) * dx_beta.at(1).at(array_index) +
          g.at(2).at(0).at(array_index) * dx_beta.at(2).at(array_index);

      // dt g_xy
      dt_g.at(0).at(1).at(array_index) =
          -2 * alpha.at(array_index) * k.at(0).at(1).at(array_index) +
          beta.at(0).at(array_index) * dx_g.at(0).at(1).at(array_index) +
          beta.at(1).at(array_index) * dy_g.at(0).at(1).at(array_index) +
          beta.at(2).at(array_index) * dz_g.at(0).at(1).at(array_index) +
          g.at(0).at(0).at(array_index) * dy_beta.at(0).at(array_index) +
          g.at(1).at(0).at(array_index) * dy_beta.at(1).at(array_index) +
          g.at(2).at(0).at(array_index) * dy_beta.at(2).at(array_index) +
          g.at(0).at(1).at(array_index) * dx_beta.at(0).at(array_index) +
          g.at(1).at(1).at(array_index) * dx_beta.at(1).at(array_index) +
          g.at(2).at(1).at(array_index) * dx_beta.at(2).at(array_index);

      // dt g_xz
      dt_g.at(0).at(2).at(array_index) =
          -2 * alpha.at(array_index) * k.at(0).at(2).at(array_index) +
          beta.at(0).at(array_index) * dx_g.at(0).at(2).at(array_index) +
          beta.at(1).at(array_index) * dy_g.at(0).at(2).at(array_index) +
          beta.at(2).at(array_index) * dz_g.at(0).at(2).at(array_index) +
          g.at(0).at(0).at(array_index) * dz_beta.at(0).at(array_index) +
          g.at(1).at(0).at(array_index) * dz_beta.at(1).at(array_index) +
          g.at(2).at(0).at(array_index) * dz_beta.at(2).at(array_index) +
          g.at(0).at(2).at(array_index) * dx_beta.at(0).at(array_index) +
          g.at(1).at(2).at(array_index) * dx_beta.at(1).at(array_index) +
          g.at(2).at(2).at(array_index) * dx_beta.at(2).at(array_index);

      // dt g_yy
      dt_g.at(1).at(1).at(array_index) =
          -2 * alpha.at(array_index) * k.at(1).at(1).at(array_index) +
          beta.at(0).at(array_index) * dx_g.at(1).at(1).at(array_index) +
          beta.at(1).at(array_index) * dy_g.at(1).at(1).at(array_index) +
          beta.at(2).at(array_index) * dz_g.at(1).at(1).at(array_index) +
          g.at(0).at(1).at(array_index) * dy_beta.at(0).at(array_index) +
          g.at(1).at(1).at(array_index) * dy_beta.at(1).at(array_index) +
          g.at(2).at(1).at(array_index) * dy_beta.at(2).at(array_index) +
          g.at(0).at(1).at(array_index) * dy_beta.at(0).at(array_index) +
          g.at(1).at(1).at(array_index) * dy_beta.at(1).at(array_index) +
          g.at(2).at(1).at(array_index) * dy_beta.at(2).at(array_index);

      // dt g_yz
      dt_g.at(1).at(2).at(array_index) =
          -2 * alpha.at(array_index) * k.at(1).at(2).at(array_index) +
          beta.at(0).at(array_index) * dx_g.at(1).at(2).at(array_index) +
          beta.at(1).at(array_index) * dy_g.at(1).at(2).at(array_index) +
          beta.at(2).at(array_index) * dz_g.at(1).at(2).at(array_index) +
          g.at(0).at(1).at(array_index) * dz_beta.at(0).at(array_index) +
          g.at(1).at(1).at(array_index) * dz_beta.at(1).at(array_index) +
          g.at(2).at(1).at(array_index) * dz_beta.at(2).at(array_index) +
          g.at(0).at(2).at(array_index) * dy_beta.at(0).at(array_index) +
          g.at(1).at(2).at(array_index) * dy_beta.at(1).at(array_index) +
          g.at(2).at(2).at(array_index) * dy_beta.at(2).at(array_index);

      // dt g_zz
      dt_g.at(2).at(2).at(array_index) =
          -2 * alpha.at(array_index) * k.at(2).at(2).at(array_index) +
          beta.at(0).at(array_index) * dx_g.at(2).at(2).at(array_index) +
          beta.at(1).at(array_index) * dy_g.at(2).at(2).at(array_index) +
          beta.at(2).at(array_index) * dz_g.at(2).at(2).at(array_index) +
          g.at(0).at(2).at(array_index) * dz_beta.at(0).at(array_index) +
          g.at(1).at(2).at(array_index) * dz_beta.at(1).at(array_index) +
          g.at(2).at(2).at(array_index) * dz_beta.at(2).at(array_index) +
          g.at(0).at(2).at(array_index) * dz_beta.at(0).at(array_index) +
          g.at(1).at(2).at(array_index) * dz_beta.at(1).at(array_index) +
          g.at(2).at(2).at(array_index) * dz_beta.at(2).at(array_index);
    }

    // Integrate to obtain spherical harmonic decomposition
    const int lmax = 8;
    const int mode_count = l_m_to_index(lmax, lmax) + 1;
    vector<vector<CCTK_REAL> > re_ylms(mode_count,
                                       vector<CCTK_REAL>(array_size));
    vector<vector<CCTK_REAL> > im_ylms(mode_count,
                                       vector<CCTK_REAL>(array_size));

    Compute_Ylms(th, ph, re_ylms, im_ylms, lmax, array_size);

    // Decompose g, dr_g, dt_g
    // re_g[i][j][mode], im_g[i][j][mode]
    vector<vector<vector<CCTK_REAL> > > re_g(
        3, vector<vector<CCTK_REAL> >(3, vector<CCTK_REAL>(mode_count)));
    vector<vector<vector<CCTK_REAL> > > im_g(
        3, vector<vector<CCTK_REAL> >(3, vector<CCTK_REAL>(mode_count)));
    vector<vector<vector<CCTK_REAL> > > re_dr_g(
        3, vector<vector<CCTK_REAL> >(3, vector<CCTK_REAL>(mode_count)));
    vector<vector<vector<CCTK_REAL> > > im_dr_g(
        3, vector<vector<CCTK_REAL> >(3, vector<CCTK_REAL>(mode_count)));
    vector<vector<vector<CCTK_REAL> > > re_dt_g(
        3, vector<vector<CCTK_REAL> >(3, vector<CCTK_REAL>(mode_count)));
    vector<vector<vector<CCTK_REAL> > > im_dt_g(
        3, vector<vector<CCTK_REAL> >(3, vector<CCTK_REAL>(mode_count)));
    for (int i = 0; i < 3; i++) {
      for (int j = i; j < 3; j++) {
        Decompose_Spherical_Harmonics(th, ph, g.at(i).at(j), re_g.at(i).at(j),
                                      im_g.at(i).at(j), re_ylms, im_ylms,
                                      array_size, lmax, ntheta, nphi);
        Decompose_Spherical_Harmonics(
            th, ph, dr_g.at(i).at(j), re_dr_g.at(i).at(j), im_dr_g.at(i).at(j),
            re_ylms, im_ylms, array_size, lmax, ntheta, nphi);
        Decompose_Spherical_Harmonics(
            th, ph, dt_g.at(i).at(j), re_dt_g.at(i).at(j), im_dt_g.at(i).at(j),
            re_ylms, im_ylms, array_size, lmax, ntheta, nphi);
      }
    }

    // Decompose beta, dr_beta, dt_beta
    // re_beta[i][mode]
    vector<vector<CCTK_REAL> > re_beta(3, vector<CCTK_REAL>(mode_count));
    vector<vector<CCTK_REAL> > im_beta(3, vector<CCTK_REAL>(mode_count));
    vector<vector<CCTK_REAL> > re_dr_beta(3, vector<CCTK_REAL>(mode_count));
    vector<vector<CCTK_REAL> > im_dr_beta(3, vector<CCTK_REAL>(mode_count));
    vector<vector<CCTK_REAL> > re_dt_beta(3, vector<CCTK_REAL>(mode_count));
    vector<vector<CCTK_REAL> > im_dt_beta(3, vector<CCTK_REAL>(mode_count));
    for (int i = 0; i < 3; i++) {
      Decompose_Spherical_Harmonics(th, ph, beta.at(i), re_beta.at(i),
                                    im_beta.at(i), re_ylms, im_ylms, array_size,
                                    lmax, ntheta, nphi);
      Decompose_Spherical_Harmonics(th, ph, dr_beta.at(i), re_dr_beta.at(i),
                                    im_dr_beta.at(i), re_ylms, im_ylms,
                                    array_size, lmax, ntheta, nphi);
      Decompose_Spherical_Harmonics(th, ph, dt_beta.at(i), re_dt_beta.at(i),
                                    im_dt_beta.at(i), re_ylms, im_ylms,
                                    array_size, lmax, ntheta, nphi);
    }

    // Decompose alpha, dr_alpha, dt_alpha
    vector<CCTK_REAL> re_alpha(mode_count);
    vector<CCTK_REAL> im_alpha(mode_count);
    vector<CCTK_REAL> re_dr_alpha(mode_count);
    vector<CCTK_REAL> im_dr_alpha(mode_count);
    vector<CCTK_REAL> re_dt_alpha(mode_count);
    vector<CCTK_REAL> im_dt_alpha(mode_count);
    Decompose_Spherical_Harmonics(th, ph, alpha, re_alpha, im_alpha, re_ylms,
                                  im_ylms, array_size, lmax, ntheta, nphi);
    Decompose_Spherical_Harmonics(th, ph, dr_alpha, re_dr_alpha, im_dr_alpha,
                                  re_ylms, im_ylms, array_size, lmax, ntheta,
                                  nphi);
    Decompose_Spherical_Harmonics(th, ph, dt_alpha, re_dt_alpha, im_dt_alpha,
                                  re_ylms, im_ylms, array_size, lmax, ntheta,
                                  nphi);

    // Store output in h5 file
    if (CCTK_MyProc(cctkGH) == 0) {
      Output_Decomposed_Metric_Data(
          CCTK_PASS_CTOC, re_g, im_g, re_dr_g, im_dr_g, re_dt_g, im_dt_g,
          re_beta, im_beta, re_dr_beta, im_dr_beta, re_dt_beta, im_dt_beta,
          re_alpha, im_alpha, re_dr_alpha, im_dr_alpha, re_dt_alpha,
          im_dt_alpha, radius[r], lmax);
    }
  }
}
} // namespace CCE_export
