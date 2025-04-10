#include "interpolate.hh"

#include <string.h>

using std::string;

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

} // namespace CCE_export
