#include "cce_export.hh"
#include <vector>

using namespace std;

void CCE_Export_Interpolate_On_Sphere(CCTK_ARGUMENTS, vector<CCTK_REAL> &xs, vector<CCTK_REAL> &ys, vector<CCTK_REAL> &zs, string name, vector<CCTK_REAL> &sphere_values, vector<CCTK_REAL> &sphere_dx, vector<CCTK_REAL> &sphere_dy, vector<CCTK_REAL> &sphere_dz, CCTK_INT array_size)
{
  printf("In interpolate\n");
  //DECLARE_CCTK_ARGUMENTS;
  //DECLARE_CCTK_PARAMETERS;

  printf("About to print name\n");
  printf("%s\n", name.c_str());

  printf("here\n");
  printf("test %d\n", 1);
  printf("%d\n", CCTK_VarIndex("ADMBase::kxx"));

  CCTK_INT variable_index = CCTK_VarIndex(name.c_str());

  printf("found variable index, %d\n", variable_index);

  CCTK_INT num_input_arrays = 1;
  CCTK_INT num_output_arrays = 4;

  const CCTK_INT num_dims = 3;
  int ierr = -1;

  const void* interp_coords[num_dims]
    = { (const void *) xs.data(),
	(const void *) ys.data(),
	(const void *) zs.data() };
  
  const CCTK_INT input_array_indices[1]
    = { variable_index };
  
  const CCTK_INT output_array_types[4]
    = { CCTK_VARIABLE_REAL,
	CCTK_VARIABLE_REAL,
	CCTK_VARIABLE_REAL,
	CCTK_VARIABLE_REAL};


  void * output_arrays[4]
    = { (void *) sphere_values.data(),
	(void *) sphere_dx.data(),
	(void *) sphere_dy.data(),
	(void *) sphere_dz.data()};

  printf("Defined variables\n");
  
  const int operator_handle = CCTK_InterpHandle("Hermite polynomial interpolation");
  
  int param_table_handle = Util_TableCreate(UTIL_TABLE_FLAGS_DEFAULT);
  ierr = Util_TableSetFromString(param_table_handle, "order=3 boundary_off_centering_tolerance={0.0 0.0 0.0 0.0 0.0 0.0} boundary_extrapolation_tolerance={0.0 0.0 0.0 0.0 0.0 0.0}");
  
  const CCTK_INT operand_indices[4] = {0, 0, 0, 0};
  Util_TableSetIntArray(param_table_handle, 4, operand_indices, "operand_indices");
  
  const CCTK_INT operation_codes[4] = {0, 1, 2, 3};
  Util_TableSetIntArray(param_table_handle, 4, operation_codes, "operation_codes");
  
  const int coord_system_handle = CCTK_CoordSystemHandle("cart3d");

  printf("about to interpolate\n");

  ierr = CCTK_InterpGridArrays(
			       cctkGH,
			       num_dims,
			       operator_handle,
			       param_table_handle,
			       coord_system_handle,
			       CCTK_MyProc(cctkGH) == 0 ? array_size : 0, // Only the 0 processor needs the points                                                                         
			       CCTK_VARIABLE_REAL,
			       interp_coords,
			       num_input_arrays,
			       input_array_indices,
			       num_output_arrays,
			       output_array_types,
			       output_arrays);
  printf("interpolated\n");

}

void CCE_Export(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  printf("CCE_Export: starting\n");

  static string index_to_component[] = {"x", "y", "z"};

  const int ntheta = 50; 
  const int nphi = 100;
  const int array_size=ntheta*nphi;

  // Interpolate onto sphere
  vector<vector<vector<CCTK_REAL>>> k(3, vector<vector<CCTK_REAL>>(3, vector<CCTK_REAL>(array_size)));
  vector<vector<vector<CCTK_REAL>>> dx_k(3, vector<vector<CCTK_REAL>>(3, vector<CCTK_REAL>(array_size)));
  vector<vector<vector<CCTK_REAL>>> dy_k(3, vector<vector<CCTK_REAL>>(3, vector<CCTK_REAL>(array_size)));
  vector<vector<vector<CCTK_REAL>>> dz_k(3, vector<vector<CCTK_REAL>>(3, vector<CCTK_REAL>(array_size)));
  vector<vector<vector<CCTK_REAL>>> g(3, vector<vector<CCTK_REAL>>(3, vector<CCTK_REAL>(array_size)));
  vector<vector<vector<CCTK_REAL>>> dx_g(3, vector<vector<CCTK_REAL>>(3, vector<CCTK_REAL>(array_size)));
  vector<vector<vector<CCTK_REAL>>> dy_g(3, vector<vector<CCTK_REAL>>(3, vector<CCTK_REAL>(array_size)));
  vector<vector<vector<CCTK_REAL>>> dz_g(3, vector<vector<CCTK_REAL>>(3, vector<CCTK_REAL>(array_size)));
  vector<vector<vector<CCTK_REAL>>> dr_g(3, vector<vector<CCTK_REAL>>(3, vector<CCTK_REAL>(array_size)));
  vector<vector<vector<CCTK_REAL>>> dt_g(3, vector<vector<CCTK_REAL>>(3, vector<CCTK_REAL>(array_size)));
  vector<vector<CCTK_REAL>> beta(3, vector<CCTK_REAL>(array_size));
  vector<vector<CCTK_REAL>> dx_beta(3, vector<CCTK_REAL>(array_size));
  vector<vector<CCTK_REAL>> dy_beta(3, vector<CCTK_REAL>(array_size));
  vector<vector<CCTK_REAL>> dz_beta(3, vector<CCTK_REAL>(array_size));
  vector<vector<CCTK_REAL>> dr_beta(3, vector<CCTK_REAL>(array_size));
  vector<vector<CCTK_REAL>> dt_beta(3, vector<CCTK_REAL>(array_size));
  vector<CCTK_REAL> alpha(array_size);
  vector<CCTK_REAL> dx_alpha(array_size);
  vector<CCTK_REAL> dy_alpha(array_size);
  vector<CCTK_REAL> dz_alpha(array_size);
  vector<CCTK_REAL> dr_alpha(array_size);
  vector<CCTK_REAL> dt_alpha(array_size);

  vector<CCTK_REAL> xs(array_size);
  vector<CCTK_REAL> ys(array_size);
  vector<CCTK_REAL> zs(array_size);
  vector<CCTK_REAL> xhat(array_size);
  vector<CCTK_REAL> yhat(array_size);
  vector<CCTK_REAL> zhat(array_size);
  vector<CCTK_REAL> th(array_size);
  vector<CCTK_REAL> ph(array_size);

  printf("CCE_Export: variables allocated\n");

  const CCTK_REAL PI = acos(-1.0);

  for (int theta_index = 0; theta_index < ntheta; theta_index++)
  {
    for (int phi_index = 0; phi_index < nphi; phi_index++)
    {
      const int array_index = theta_index + ntheta * phi_index;
      printf("theta index: %d\n", theta_index);
      printf("phi index: %d\n", phi_index);
      printf("array index: %d\n", array_index);
	
      th.at(array_index) = theta_index * PI / (ntheta);
      ph.at(array_index) = phi_index * 2 * PI / nphi;
      xhat.at(array_index) = cos(ph.at(array_index))*sin(th.at(array_index));
      yhat.at(array_index) = sin(ph.at(array_index))*sin(th.at(array_index));
      zhat.at(array_index) = cos(th.at(array_index));
    }
  }

  printf("CCE_Export: th, phi, xhat, yhat, zhat initialized\n");

  for(int r=0; r<nradii; r++)
  {
    printf("CCE_Export: in radius loop\n");

    for (int theta_index = 0; theta_index < ntheta; theta_index++)
    {
      for (int phi_index = 0; phi_index < nphi; phi_index++)
      {
	const int array_index = theta_index + ntheta * phi_index;
	
	xs.at(array_index) = radius[r] * xhat.at(array_index);
	ys.at(array_index) = radius[r] * yhat.at(array_index);
	zs.at(array_index) = radius[r] * zhat.at(array_index);
      }
    }
    printf("CCE_EXPORT: xs, ys, zs initialized\n");
    
    for(int i=0; i<3; i++){
      string first_component = index_to_component[i];
      for(int j=i; j<3; j++){
	string second_component = index_to_component[j];
	printf("About to interpolated k\n");
	CCE_Export_Interpolate_On_Sphere(CCTK_PASS_CTOC, xs, ys, zs, "ADMBase::k" + first_component + second_component , k.at(i).at(j), dx_k.at(i).at(j), dy_k.at(i).at(j), dz_k.at(i).at(j), array_size);
	printf("Inerpolated k\n");
	CCE_Export_Interpolate_On_Sphere(CCTK_PASS_CTOC, xs, ys, zs, "ADMBase::g" + first_component + second_component , g.at(i).at(j), dx_g.at(i).at(j), dy_g.at(i).at(j), dz_g.at(i).at(j), array_size);
      }
      CCE_Export_Interpolate_On_Sphere(CCTK_PASS_CTOC, xs, ys, zs, "ADMBase::beta" + first_component, beta.at(i), dx_beta.at(i), dy_beta.at(i), dz_beta.at(i), array_size);
    }
    CCE_Export_Interpolate_On_Sphere(CCTK_PASS_CTOC, xs, ys, zs, "ADMBase::alp", alpha, dx_alpha, dy_alpha, dz_alpha, array_size);

    // Integrate to obtain spherical harmonic decomposition

    // Store output in h5 file
  }

}
