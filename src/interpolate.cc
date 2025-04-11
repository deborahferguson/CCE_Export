#include "interpolate.hh"

namespace CCE_export {

void Interpolate_On_Sphere_With_Derivatives(
    CCTK_ARGUMENTS, vector<CCTK_REAL> &xs, vector<CCTK_REAL> &ys,
    vector<CCTK_REAL> &zs, string name, vector<CCTK_REAL> &sphere_values,
    vector<CCTK_REAL> &sphere_dx, vector<CCTK_REAL> &sphere_dy,
    vector<CCTK_REAL> &sphere_dz, CCTK_INT array_size) {

  CCTK_INT variable_index = CCTK_VarIndex(name.c_str());

  CCTK_INT num_input_arrays = 1;
  CCTK_INT num_output_arrays = 4;

  const CCTK_INT num_dims = 3;
  int ierr = -1;

  const void *interp_coords[num_dims] = {(const void *)xs.data(),
                                         (const void *)ys.data(),
                                         (const void *)zs.data()};

  const CCTK_INT input_array_indices[1] = {variable_index};

  const CCTK_INT output_array_types[4] = {
      CCTK_VARIABLE_REAL, CCTK_VARIABLE_REAL, CCTK_VARIABLE_REAL,
      CCTK_VARIABLE_REAL};

  void *output_arrays[4] = {(void *)sphere_values.data(),
                            (void *)sphere_dx.data(), (void *)sphere_dy.data(),
                            (void *)sphere_dz.data()};

  const int operator_handle =
      CCTK_InterpHandle("Hermite polynomial interpolation");

  int param_table_handle = Util_TableCreate(UTIL_TABLE_FLAGS_DEFAULT);
  ierr = Util_TableSetFromString(
      param_table_handle,
      "order=3 boundary_off_centering_tolerance={0.0 0.0 0.0 0.0 0.0 0.0} "
      "boundary_extrapolation_tolerance={0.0 0.0 0.0 0.0 0.0 0.0}");

  const CCTK_INT operand_indices[4] = {0, 0, 0, 0};
  Util_TableSetIntArray(param_table_handle, 4, operand_indices,
                        "operand_indices");

  const CCTK_INT operation_codes[4] = {0, 1, 2, 3};
  Util_TableSetIntArray(param_table_handle, 4, operation_codes,
                        "operation_codes");

  const int coord_system_handle = CCTK_CoordSystemHandle("cart3d");

  ierr = CCTK_InterpGridArrays(
      cctkGH, num_dims, operator_handle, param_table_handle,
      coord_system_handle,
      CCTK_MyProc(cctkGH) == 0 ? array_size
                               : 0, // Only the 0 processor needs the points
      CCTK_VARIABLE_REAL, interp_coords, num_input_arrays, input_array_indices,
      num_output_arrays, output_array_types, output_arrays);
}

void Interpolate_On_Sphere(CCTK_ARGUMENTS, vector<CCTK_REAL> &xs,
                           vector<CCTK_REAL> &ys, vector<CCTK_REAL> &zs,
                           string name, vector<CCTK_REAL> &sphere_values,
                           CCTK_INT array_size) {

  CCTK_INT variable_index = CCTK_VarIndex(name.c_str());

  CCTK_INT num_input_arrays = 1;
  CCTK_INT num_output_arrays = 1;

  const CCTK_INT num_dims = 3;
  int ierr = -1;

  const void *interp_coords[num_dims] = {(const void *)xs.data(),
                                         (const void *)ys.data(),
                                         (const void *)zs.data()};

  const CCTK_INT input_array_indices[1] = {variable_index};

  const CCTK_INT output_array_types[1] = {CCTK_VARIABLE_REAL};

  void *output_arrays[1] = {(void *)sphere_values.data()};

  const int operator_handle =
      CCTK_InterpHandle("Hermite polynomial interpolation");

  int param_table_handle = Util_TableCreate(UTIL_TABLE_FLAGS_DEFAULT);
  ierr = Util_TableSetFromString(
      param_table_handle,
      "order=3 boundary_off_centering_tolerance={0.0 0.0 0.0 0.0 0.0 0.0} "
      "boundary_extrapolation_tolerance={0.0 0.0 0.0 0.0 0.0 0.0}");

  const CCTK_INT operand_indices[1] = {0};
  Util_TableSetIntArray(param_table_handle, 1, operand_indices,
                        "operand_indices");

  const CCTK_INT operation_codes[1] = {0};
  Util_TableSetIntArray(param_table_handle, 1, operation_codes,
                        "operation_codes");

  const int coord_system_handle = CCTK_CoordSystemHandle("cart3d");

  ierr = CCTK_InterpGridArrays(
      cctkGH, num_dims, operator_handle, param_table_handle,
      coord_system_handle,
      CCTK_MyProc(cctkGH) == 0 ? array_size
                               : 0, // Only the 0 processor needs the points
      CCTK_VARIABLE_REAL, interp_coords, num_input_arrays, input_array_indices,
      num_output_arrays, output_array_types, output_arrays);
}

void Extract_Metric_Shift_Lapse_On_Sphere(
    CCTK_ARGUMENTS, vector<vector<vector<CCTK_REAL> > > &k,
    vector<vector<vector<CCTK_REAL> > > &dx_k,
    vector<vector<vector<CCTK_REAL> > > &dy_k,
    vector<vector<vector<CCTK_REAL> > > &dz_k,
    vector<vector<vector<CCTK_REAL> > > &g,
    vector<vector<vector<CCTK_REAL> > > &dx_g,
    vector<vector<vector<CCTK_REAL> > > &dy_g,
    vector<vector<vector<CCTK_REAL> > > &dz_g,
    vector<vector<vector<CCTK_REAL> > > &dr_g,
    vector<vector<vector<CCTK_REAL> > > &dt_g, vector<vector<CCTK_REAL> > &beta,
    vector<vector<CCTK_REAL> > &dx_beta, vector<vector<CCTK_REAL> > &dy_beta,
    vector<vector<CCTK_REAL> > &dz_beta, vector<vector<CCTK_REAL> > &dr_beta,
    vector<vector<CCTK_REAL> > &dt_beta, vector<CCTK_REAL> &alpha,
    vector<CCTK_REAL> &dx_alpha, vector<CCTK_REAL> &dy_alpha,
    vector<CCTK_REAL> &dz_alpha, vector<CCTK_REAL> &dr_alpha,
    vector<CCTK_REAL> &dt_alpha, vector<CCTK_REAL> &th, vector<CCTK_REAL> &ph,
    vector<CCTK_REAL> &xhat, vector<CCTK_REAL> &yhat, vector<CCTK_REAL> &zhat,
    vector<CCTK_REAL> &xs, vector<CCTK_REAL> &ys, vector<CCTK_REAL> &zs,
    int ntheta, int nphi, int array_size, int rad_index) {

  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  static string index_to_component[] = {"x", "y", "z"};

  // compute the values of x, y, z and the desired points on the sphere of
  // radius radius[r]
  for (int theta_index = 0; theta_index <= ntheta; theta_index++) {
    for (int phi_index = 0; phi_index <= nphi; phi_index++) {
      const int array_index = theta_index + (ntheta + 1) * phi_index;

      xs.at(array_index) = radius[rad_index] * xhat.at(array_index);
      ys.at(array_index) = radius[rad_index] * yhat.at(array_index);
      zs.at(array_index) = radius[rad_index] * zhat.at(array_index);
    }
  }

  // Interpolate all the desired quantities onto the desired points on the
  // sphere

  // interpolate metric and extrinsic curvature
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
            (xs.at(array_index) / radius[rad_index]) *
                dx_g.at(i).at(j).at(array_index) +
            (ys.at(array_index) / radius[rad_index]) *
                dy_g.at(i).at(j).at(array_index) +
            (zs.at(array_index) / radius[rad_index]) * dz_g.at(i).at(j).at(array_index);
      }
    }
  }

  // interpolate shift and its derivatives
  for (int i = 0; i < 3; i++) {
    string component = index_to_component[i];
    // interpolate shift
    Interpolate_On_Sphere_With_Derivatives(
        CCTK_PASS_CTOC, xs, ys, zs, "ADMBase::beta" + component, beta.at(i),
        dx_beta.at(i), dy_beta.at(i), dz_beta.at(i), array_size);
    // interpolate time derivative of shift
    Interpolate_On_Sphere(CCTK_PASS_CTOC, xs, ys, zs,
                          "ADMBase::dtbeta" + component, dt_beta.at(i),
                          array_size);
    // compute dr_beta
    for (int array_index = 0; array_index < array_size; array_index++) {
      dr_beta.at(i).at(array_index) =
          (xs.at(array_index) / radius[rad_index]) * dx_beta.at(i).at(array_index) +
          (ys.at(array_index) / radius[rad_index]) * dy_beta.at(i).at(array_index) +
          (zs.at(array_index) / radius[rad_index]) * dz_beta.at(i).at(array_index);
    }
  }

  // interpolate lapse
  Interpolate_On_Sphere_With_Derivatives(CCTK_PASS_CTOC, xs, ys, zs,
                                         "ADMBase::alp", alpha, dx_alpha,
                                         dy_alpha, dz_alpha, array_size);
  // interpolate time derivative of lapse
  Interpolate_On_Sphere(CCTK_PASS_CTOC, xs, ys, zs, "ADMBase::dtalp", dt_alpha,
                        array_size);
  // compute dr_alpha
  for (int array_index = 0; array_index < array_size; array_index++) {
    dr_alpha.at(array_index) =
        (xs.at(array_index) / radius[rad_index]) * dx_alpha.at(array_index) +
        (ys.at(array_index) / radius[rad_index]) * dy_alpha.at(array_index) +
        (zs.at(array_index) / radius[rad_index]) * dz_alpha.at(array_index);
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
}

} // namespace CCE_export
