#include "cce_export.hh"
#include <vector>

using namespace std;

// Copied from Multipole
#define idx(xx,yy) (assert((xx) <= nx), assert((xx) >= 0), assert((yy) <= ny), assert((yy) >= 0), ((xx) + (yy) * (nx+1)))

// Copied from Multipole
CCTK_REAL CCE_Export_Simpson2DIntegral(CCTK_REAL const *f, int nx, int ny, CCTK_REAL hx, CCTK_REAL hy)
{
  CCTK_REAL integrand_sum = 0;
  int ix = 0, iy = 0;

  assert(nx > 0); assert(ny > 0); assert (f);
  assert(nx % 2 == 0);
  assert(ny % 2 == 0);

  int px = nx / 2;
  int py = ny / 2;

  // Corners                                                                                                                                                                                          
  integrand_sum += f[idx(0,0)] + f[idx(nx,0)] + f[idx(0,ny)] + f[idx(nx,ny)];

  // Edges                                                                                                                                                                                            
  for (iy = 1; iy <= py; iy++)
    integrand_sum += 4 * f[idx(0,2*iy-1)] + 4 * f[idx(nx,2*iy-1)];

  for (iy = 1; iy <= py-1; iy++)
    integrand_sum += 2 * f[idx(0,2*iy)] + 2 * f[idx(nx,2*iy)];

  for (ix = 1; ix <= px; ix++)
    integrand_sum += 4 * f[idx(2*ix-1,0)] + 4 * f[idx(2*ix-1,ny)];

  for (ix = 1; ix <= px-1; ix++)
    integrand_sum += 2 * f[idx(2*ix,0)] + 2 * f[idx(2*ix,ny)];

  // Interior                                                                                                                                                                                         
  for (iy = 1; iy <= py; iy++)
    for (ix = 1; ix <= px; ix++)
      integrand_sum += 16 * f[idx(2*ix-1,2*iy-1)];

  for (iy = 1; iy <= py-1; iy++)
    for (ix = 1; ix <= px; ix++)
      integrand_sum += 8 * f[idx(2*ix-1,2*iy)];

  for (iy = 1; iy <= py; iy++)
    for (ix = 1; ix <= px-1; ix++)
      integrand_sum += 8 * f[idx(2*ix,2*iy-1)];

  for (iy = 1; iy <= py-1; iy++)
    for (ix = 1; ix <= px-1; ix++)
      integrand_sum += 4 * f[idx(2*ix,2*iy)];

  return ((double) 1 / (double) 9) * hx * hy * integrand_sum;
}


void CCE_Export_Interpolate_On_Sphere_With_Derivatives(CCTK_ARGUMENTS, vector<CCTK_REAL> &xs, vector<CCTK_REAL> &ys, vector<CCTK_REAL> &zs, string name, vector<CCTK_REAL> &sphere_values, vector<CCTK_REAL> &sphere_dx, vector<CCTK_REAL> &sphere_dy, vector<CCTK_REAL> &sphere_dz, CCTK_INT array_size)
{
  printf("In interpolate\n");

  CCTK_INT variable_index = CCTK_VarIndex(name.c_str());


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
  
  const int operator_handle = CCTK_InterpHandle("Hermite polynomial interpolation");
  
  int param_table_handle = Util_TableCreate(UTIL_TABLE_FLAGS_DEFAULT);
  ierr = Util_TableSetFromString(param_table_handle, "order=3 boundary_off_centering_tolerance={0.0 0.0 0.0 0.0 0.0 0.0} boundary_extrapolation_tolerance={0.0 0.0 0.0 0.0 0.0 0.0}");
  
  const CCTK_INT operand_indices[4] = {0, 0, 0, 0};
  Util_TableSetIntArray(param_table_handle, 4, operand_indices, "operand_indices");
  
  const CCTK_INT operation_codes[4] = {0, 1, 2, 3};
  Util_TableSetIntArray(param_table_handle, 4, operation_codes, "operation_codes");
  
  const int coord_system_handle = CCTK_CoordSystemHandle("cart3d");

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

}

void CCE_Export_Interpolate_On_Sphere(CCTK_ARGUMENTS, vector<CCTK_REAL> &xs, vector<CCTK_REAL> &ys, vector<CCTK_REAL> &zs, string name, vector<CCTK_REAL> &sphere_values, CCTK_INT array_size)
{
  printf("In interpolate\n");

  CCTK_INT variable_index = CCTK_VarIndex(name.c_str());


  CCTK_INT num_input_arrays = 1;
  CCTK_INT num_output_arrays = 1;

  const CCTK_INT num_dims = 3;
  int ierr = -1;

  const void* interp_coords[num_dims]
    = { (const void *) xs.data(),
	(const void *) ys.data(),
	(const void *) zs.data() };
  
  const CCTK_INT input_array_indices[1]
    = { variable_index };
  
  const CCTK_INT output_array_types[1]
    = { CCTK_VARIABLE_REAL};


  void * output_arrays[1]
    = { (void *) sphere_values.data()};
  
  const int operator_handle = CCTK_InterpHandle("Hermite polynomial interpolation");
  
  int param_table_handle = Util_TableCreate(UTIL_TABLE_FLAGS_DEFAULT);
  ierr = Util_TableSetFromString(param_table_handle, "order=3 boundary_off_centering_tolerance={0.0 0.0 0.0 0.0 0.0 0.0} boundary_extrapolation_tolerance={0.0 0.0 0.0 0.0 0.0 0.0}");
  
  const CCTK_INT operand_indices[1] = {0};
  Util_TableSetIntArray(param_table_handle, 1, operand_indices, "operand_indices");
  
  const CCTK_INT operation_codes[1] = {0};
  Util_TableSetIntArray(param_table_handle, 1, operation_codes, "operation_codes");
  
  const int coord_system_handle = CCTK_CoordSystemHandle("cart3d");

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

}

static inline int CCE_Export_Sphere_Index(int it, int ip, int ntheta)
{
  return it + (ntheta+1)*ip;
}

// Copied from Multipole
void CCE_Export_Integrate(int array_size, int ntheta, int nphi,
			 vector<CCTK_REAL> &array1r, vector<CCTK_REAL> &array1i,
			 vector<CCTK_REAL> &array2r,
			 vector<CCTK_REAL> &th, vector<CCTK_REAL> &ph,
			 CCTK_REAL *outre, CCTK_REAL *outim)
{
  DECLARE_CCTK_PARAMETERS

  int il = CCE_Export_Sphere_Index(0,0,ntheta);
  int iu = CCE_Export_Sphere_Index(1,0,ntheta);
  CCTK_REAL dth = th[iu] - th[il];
  iu = CCE_Export_Sphere_Index(0,1,ntheta);
  CCTK_REAL dph = ph[iu] - ph[il];

  static CCTK_REAL *fr = 0;
  static CCTK_REAL *fi = 0;
  static bool allocated_memory = false;

  // Construct an array for the real integrand                                                                                                                                                        
  if (!allocated_memory)
    {
      fr = new CCTK_REAL[array_size];
      fi = new CCTK_REAL[array_size];
      allocated_memory = true;
    }

  // the below calculations take the integral of conj(array1)*array2*sin(th)                                                                                                                          
  for (int i = 0; i < array_size; i++)
    {
      fr[i] = (array1r[i] * array2r[i]) * sin(th[i]);
      fi[i] = ( -1 * array1i[i] * array2r[i] ) * sin(th[i]);
    }

  if (nphi % 2 != 0 || ntheta % 2 != 0)
    {
      CCTK_WARN (CCTK_WARN_ABORT, "The Simpson integration method requires even ntheta and even nphi");
    }
  *outre = CCE_Export_Simpson2DIntegral(fr, ntheta, nphi, dth, dph);
  *outim = CCE_Export_Simpson2DIntegral(fi, ntheta, nphi, dth, dph);
} 

void Decompose_Spherical_Harmonics(vector<CCTK_REAL> &th, vector<CCTK_REAL> &phi, vector<CCTK_REAL> &sphere_values, vector<CCTK_REAL> &re_data, vector<CCTK_REAL> &im_data, vector<vector<CCTK_REAL>> &re_ylms, vector<vector<CCTK_REAL>> &im_ylms, int array_size, int lmax, int ntheta, int nphi)
{
  for(int l=0; l<lmax+1; l++){
    for(int m=-l; m<l+1; m++){
      int mode_index = l_m_to_index(l, m);
      CCE_Export_Integrate(array_size, ntheta, nphi, re_ylms.at(mode_index), im_ylms.at(mode_index), sphere_values, th, phi, &re_data.at(mode_index), &im_data.at(mode_index));
    }
  }
}

// void index_to_l_m(int index, int &l, int &m){
  
// }

int l_m_to_index(int l, int m)
{
  return l*l + l + m;
}

CCTK_REAL factorial(CCTK_REAL x)
{
  CCTK_REAL answer = 1;
  
  while(x>0){
    answer *= x;
    x-=1;
  }
  return answer;
}

static inline int imin(int a, int b)
{
  return a < b ? a : b;
}

static inline int imax(int a, int b)
{
  return a > b ? a : b;
}

static inline double combination(int n, int m)
{
  // Binomial coefficient is undefined if these conditions do not hold                                                                                                                                                                                                                                                                                                                     
  assert(n >= 0);
  assert(m >= 0);
  assert(m <= n);

  return factorial(n) / (factorial(m) * factorial(n-m));
}

void Compute_Ylms(vector<CCTK_REAL> &th, vector<CCTK_REAL> &ph, vector<vector<CCTK_REAL>> &re_ylms, vector<vector<CCTK_REAL>> &im_ylms, int lmax, int array_size)
{
  const CCTK_REAL PI = acos(-1.0);
  for(int l=0; l<lmax+1; l++){
    for(int m=-l; m<l+1; m++){
      int ylm_index = l_m_to_index(l, m);
      for(int array_index=0; array_index<array_size; array_index++){
	double all_coeff = 0, sum = 0;
	all_coeff = pow(-1.0, m);
	all_coeff *= sqrt(factorial(l+m)*factorial(l-m)*(2*l+1) / (4.*PI*factorial(l)*factorial(l)));
	sum = 0.;
	for(int i = imax(m, 0); i <= imin(l + m, l); i++){
	  double sum_coeff = combination(l, i) * combination(l, i-m);
	  sum += sum_coeff * pow(-1.0, l-i) * pow(cos(th[array_index]/2.), 2 * i - m) *
	    pow(sin(th[array_index]/2.), 2*(l-i)+m);
	}
	re_ylms.at(ylm_index).at(array_index) = all_coeff*sum*cos(m*ph[array_index]);
	im_ylms.at(ylm_index).at(array_index) = all_coeff*sum*sin(m*ph[array_index]);
      }
    }
  }
}

void CCE_Export(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  printf("CCE_Export: starting\n");

  static string index_to_component[] = {"x", "y", "z"};

  const int ntheta = 50; 
  const int nphi = 100;
  const int array_size = (ntheta+1)*(nphi+1);
  // const int array_size=ntheta*nphi;

  // extrinsic curvature, 3d vector (3, 3, array_size)
  vector<vector<vector<CCTK_REAL>>> k(3, vector<vector<CCTK_REAL>>(3, vector<CCTK_REAL>(array_size)));
  vector<vector<vector<CCTK_REAL>>> dx_k(3, vector<vector<CCTK_REAL>>(3, vector<CCTK_REAL>(array_size)));
  vector<vector<vector<CCTK_REAL>>> dy_k(3, vector<vector<CCTK_REAL>>(3, vector<CCTK_REAL>(array_size)));
  vector<vector<vector<CCTK_REAL>>> dz_k(3, vector<vector<CCTK_REAL>>(3, vector<CCTK_REAL>(array_size)));
  // metric, 3d vector (3, 3, array_size)
  vector<vector<vector<CCTK_REAL>>> g(3, vector<vector<CCTK_REAL>>(3, vector<CCTK_REAL>(array_size)));
  vector<vector<vector<CCTK_REAL>>> dx_g(3, vector<vector<CCTK_REAL>>(3, vector<CCTK_REAL>(array_size)));
  vector<vector<vector<CCTK_REAL>>> dy_g(3, vector<vector<CCTK_REAL>>(3, vector<CCTK_REAL>(array_size)));
  vector<vector<vector<CCTK_REAL>>> dz_g(3, vector<vector<CCTK_REAL>>(3, vector<CCTK_REAL>(array_size)));
  vector<vector<vector<CCTK_REAL>>> dr_g(3, vector<vector<CCTK_REAL>>(3, vector<CCTK_REAL>(array_size)));
  vector<vector<vector<CCTK_REAL>>> dt_g(3, vector<vector<CCTK_REAL>>(3, vector<CCTK_REAL>(array_size)));
  // shift (beta), 2d vector (3, array_size)
  vector<vector<CCTK_REAL>> beta(3, vector<CCTK_REAL>(array_size));
  vector<vector<CCTK_REAL>> dx_beta(3, vector<CCTK_REAL>(array_size));
  vector<vector<CCTK_REAL>> dy_beta(3, vector<CCTK_REAL>(array_size));
  vector<vector<CCTK_REAL>> dz_beta(3, vector<CCTK_REAL>(array_size));
  vector<vector<CCTK_REAL>> dr_beta(3, vector<CCTK_REAL>(array_size));
  vector<vector<CCTK_REAL>> dt_beta(3, vector<CCTK_REAL>(array_size));
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

  printf("CCE_Export: variables allocated\n");


  // Compute the theta and phi points as well as the corresponding x, y, z unit vectors 
  // Based on the number of theta and phi points desired (ntheta, nphi)
  const CCTK_REAL PI = acos(-1.0);

  for (int theta_index = 0; theta_index < ntheta; theta_index++)
  {
    for (int phi_index = 0; phi_index < nphi; phi_index++)
    {
      const int array_index = theta_index + ntheta * phi_index;
	
      th.at(array_index) = theta_index * PI / (ntheta);
      ph.at(array_index) = phi_index * 2 * PI / nphi;
      xhat.at(array_index) = sin(th.at(array_index))*cos(ph.at(array_index));
      yhat.at(array_index) = sin(th.at(array_index))*sin(ph.at(array_index));
      zhat.at(array_index) = cos(th.at(array_index));
    }
  }

  printf("CCE_Export: th, phi, xhat, yhat, zhat initialized\n");

  // loop through the desired radii
  for(int r=0; r<nradii; r++)
  {
    printf("CCE_Export: in radius loop\n");

    // compute the values of x, y, z and the desired points on the sphere of radius radius[r]
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
    

    // Interpolate all the desired quantities onto the desired points on the sphere
    for(int i=0; i<3; i++){
      string first_component = index_to_component[i];
      for(int j=i; j<3; j++){
	string second_component = index_to_component[j];
	// interpolate extrinsic curvature
	CCE_Export_Interpolate_On_Sphere_With_Derivatives(CCTK_PASS_CTOC, xs, ys, zs, "ADMBase::k" + first_component + second_component , k.at(i).at(j), dx_k.at(i).at(j), dy_k.at(i).at(j), dz_k.at(i).at(j), array_size);
	// interpolate metric
	CCE_Export_Interpolate_On_Sphere_With_Derivatives(CCTK_PASS_CTOC, xs, ys, zs, "ADMBase::g" + first_component + second_component , g.at(i).at(j), dx_g.at(i).at(j), dy_g.at(i).at(j), dz_g.at(i).at(j), array_size);
	// compute dr_g
	for(int array_index=0; array_index<array_size; array_index++){
	  dr_g.at(i).at(j).at(array_index) = (xs.at(array_index)/radius[r])*dx_g.at(i).at(j).at(array_index) + (ys.at(array_index)/radius[r])*dy_g.at(i).at(j).at(array_index) + (zs.at(array_index)/radius[r])*dz_g.at(i).at(j).at(array_index);
	}
      }
      // interpolate shift
      CCE_Export_Interpolate_On_Sphere_With_Derivatives(CCTK_PASS_CTOC, xs, ys, zs, "ADMBase::beta" + first_component, beta.at(i), dx_beta.at(i), dy_beta.at(i), dz_beta.at(i), array_size);
      // interpolate time derivative of shift
      CCE_Export_Interpolate_On_Sphere(CCTK_PASS_CTOC, xs, ys, zs, "ADMBase::dtbeta" + first_component, dt_beta.at(i), array_size);
      // compute dr_beta
      for(int array_index=0; array_index<array_size; array_index++){
	dr_beta.at(i).at(array_index) = (xs.at(array_index)/radius[r])*dx_beta.at(i).at(array_index) + (ys.at(array_index)/radius[r])*dy_beta.at(i).at(array_index) + (zs.at(array_index)/radius[r])*dz_beta.at(i).at(array_index);
      }

    }
    // interpolate lapse
    CCE_Export_Interpolate_On_Sphere_With_Derivatives(CCTK_PASS_CTOC, xs, ys, zs, "ADMBase::alp", alpha, dx_alpha, dy_alpha, dz_alpha, array_size);
    // interpolate time derivative of lapse
    CCE_Export_Interpolate_On_Sphere(CCTK_PASS_CTOC, xs, ys, zs, "ADMBase::dtalp", dt_alpha, array_size);
    // compute dr_alpha
    for(int array_index=0; array_index<array_size; array_index++){
      dr_alpha.at(array_index) = (xs.at(array_index)/radius[r])*dx_alpha.at(array_index) + (ys.at(array_index)/radius[r])*dy_alpha.at(array_index) + (zs.at(array_index)/radius[r])*dz_alpha.at(array_index);
    }

    // compute time derivatives of the metric
    // d_t g_ij = -2 alpha K_ij 
    //            + beta^x dx g_ij + beta^y dy g_ij + beta^z dz g_ij 
    //            + g_xi dj beta^x + g_yi dj beta^y + g_zi dj beta^z
    //            + g_xj di beta^x + g_yj di beta^y + g_zj di beta^z
    
    // define redundant metric components
    for(int array_index=0; array_index<array_size; array_index++){
      g.at(1).at(0).at(array_index) = g.at(0).at(1).at(array_index);
      g.at(2).at(0).at(array_index) = g.at(0).at(2).at(array_index);
      g.at(2).at(1).at(array_index) = g.at(1).at(2).at(array_index);
    }

    // dt g_xx
    for(int array_index=0; array_index<array_size; array_index++){
      dt_g.at(0).at(0).at(array_index) = -2*alpha.at(array_index)*k.at(0).at(0).at(array_index) + \
	beta.at(0).at(array_index)*dx_g.at(0).at(0).at(array_index) +	\
	beta.at(1).at(array_index)*dy_g.at(0).at(0).at(array_index) +	\
	beta.at(2).at(array_index)*dz_g.at(0).at(0).at(array_index) +	\
	g.at(0).at(0).at(array_index)*dx_beta.at(0).at(array_index) +	\
	g.at(1).at(0).at(array_index)*dx_beta.at(1).at(array_index) +	\
	g.at(2).at(0).at(array_index)*dx_beta.at(2).at(array_index) +	\
	g.at(0).at(0).at(array_index)*dx_beta.at(0).at(array_index) +	\
	g.at(1).at(0).at(array_index)*dx_beta.at(1).at(array_index) +	\
	g.at(2).at(0).at(array_index)*dx_beta.at(2).at(array_index);

      dt_g.at(0).at(1).at(array_index) = -2*alpha.at(array_index)*k.at(0).at(1).at(array_index) + \
	beta.at(0).at(array_index)*dx_g.at(0).at(1).at(array_index) +	\
	beta.at(1).at(array_index)*dy_g.at(0).at(1).at(array_index) +	\
	beta.at(2).at(array_index)*dz_g.at(0).at(1).at(array_index) +	\
	g.at(0).at(0).at(array_index)*dy_beta.at(0).at(array_index) +	\
	g.at(1).at(0).at(array_index)*dy_beta.at(1).at(array_index) +	\
	g.at(2).at(0).at(array_index)*dy_beta.at(2).at(array_index) +	\
	g.at(0).at(1).at(array_index)*dx_beta.at(0).at(array_index) +	\
	g.at(1).at(1).at(array_index)*dx_beta.at(1).at(array_index) +	\
	g.at(2).at(1).at(array_index)*dx_beta.at(2).at(array_index);

      dt_g.at(0).at(2).at(array_index) = -2*alpha.at(array_index)*k.at(0).at(2).at(array_index) + \
	beta.at(0).at(array_index)*dx_g.at(0).at(2).at(array_index) +	\
	beta.at(1).at(array_index)*dy_g.at(0).at(2).at(array_index) +	\
	beta.at(2).at(array_index)*dz_g.at(0).at(2).at(array_index) +	\
	g.at(0).at(0).at(array_index)*dz_beta.at(0).at(array_index) +	\
	g.at(1).at(0).at(array_index)*dz_beta.at(1).at(array_index) +	\
	g.at(2).at(0).at(array_index)*dz_beta.at(2).at(array_index) +	\
	g.at(0).at(2).at(array_index)*dx_beta.at(0).at(array_index) +	\
	g.at(1).at(2).at(array_index)*dx_beta.at(1).at(array_index) +	\
	g.at(2).at(2).at(array_index)*dx_beta.at(2).at(array_index);

      dt_g.at(1).at(1).at(array_index) = -2*alpha.at(array_index)*k.at(1).at(1).at(array_index) + \
	beta.at(0).at(array_index)*dx_g.at(1).at(1).at(array_index) +	\
	beta.at(1).at(array_index)*dy_g.at(1).at(1).at(array_index) +	\
	beta.at(2).at(array_index)*dz_g.at(1).at(1).at(array_index) +	\
	g.at(0).at(1).at(array_index)*dy_beta.at(0).at(array_index) +	\
	g.at(1).at(1).at(array_index)*dy_beta.at(1).at(array_index) +	\
	g.at(2).at(1).at(array_index)*dy_beta.at(2).at(array_index) +	\
	g.at(0).at(1).at(array_index)*dy_beta.at(0).at(array_index) +	\
	g.at(1).at(1).at(array_index)*dy_beta.at(1).at(array_index) +	\
	g.at(2).at(1).at(array_index)*dy_beta.at(2).at(array_index);

      dt_g.at(1).at(2).at(array_index) = -2*alpha.at(array_index)*k.at(1).at(2).at(array_index) + \
	beta.at(0).at(array_index)*dx_g.at(1).at(2).at(array_index) +	\
	beta.at(1).at(array_index)*dy_g.at(1).at(2).at(array_index) +	\
	beta.at(2).at(array_index)*dz_g.at(1).at(2).at(array_index) +	\
	g.at(0).at(1).at(array_index)*dz_beta.at(0).at(array_index) +	\
	g.at(1).at(1).at(array_index)*dz_beta.at(1).at(array_index) +	\
	g.at(2).at(1).at(array_index)*dz_beta.at(2).at(array_index) +	\
	g.at(0).at(2).at(array_index)*dy_beta.at(0).at(array_index) +	\
	g.at(1).at(2).at(array_index)*dy_beta.at(1).at(array_index) +	\
	g.at(2).at(2).at(array_index)*dy_beta.at(2).at(array_index);

      dt_g.at(2).at(2).at(array_index) = -2*alpha.at(array_index)*k.at(2).at(2).at(array_index) + \
	beta.at(0).at(array_index)*dx_g.at(2).at(2).at(array_index) +	\
	beta.at(1).at(array_index)*dy_g.at(2).at(2).at(array_index) +	\
	beta.at(2).at(array_index)*dz_g.at(2).at(2).at(array_index) +	\
	g.at(0).at(2).at(array_index)*dz_beta.at(0).at(array_index) +	\
	g.at(1).at(2).at(array_index)*dz_beta.at(1).at(array_index) +	\
	g.at(2).at(2).at(array_index)*dz_beta.at(2).at(array_index) +	\
	g.at(0).at(2).at(array_index)*dz_beta.at(0).at(array_index) +	\
	g.at(1).at(2).at(array_index)*dz_beta.at(1).at(array_index) +	\
	g.at(2).at(2).at(array_index)*dz_beta.at(2).at(array_index);

      dt_g.at(1).at(0).at(array_index) = dt_g.at(0).at(1).at(array_index);

      dt_g.at(2).at(0).at(array_index) = dt_g.at(0).at(2).at(array_index);

      dt_g.at(2).at(1).at(array_index) = dt_g.at(1).at(2).at(array_index);
    
    }
    
    // Integrate to obtain spherical harmonic decomposition
    const int lmax = 8;
    const int mode_count = l_m_to_index(lmax, lmax) + 1;
    vector<vector<CCTK_REAL>> re_ylms(mode_count, vector<CCTK_REAL>(array_size));
    vector<vector<CCTK_REAL>> im_ylms(mode_count, vector<CCTK_REAL>(array_size));
    
    printf("About to compute ylms\n");
    Compute_Ylms(th, ph, re_ylms, im_ylms, lmax, array_size);

    // // print 0 0 mode
    // printf("Ylms for 0 0 mode\n");
    // printf("th\tph\tre(ylm)\tim(ylm)\n");
    // int mode_index = l_m_to_index(0, 0);
    // for(int array_index=0; array_index<array_size; array_index++){
    //   printf("%f\t%f\t%f\t%f\n", th[array_index], ph[array_index], re_ylms.at(mode_index).at(array_index), im_ylms.at(mode_index).at(array_index));
    // }

    // // print 1 -1 mode
    // printf("Ylms for 1 -1 mode\n");
    // printf("th\tph\tre(ylm)\tim(ylm)\n");
    // mode_index = l_m_to_index(1, -1);
    // for(int array_index=0; array_index<array_size; array_index++){
    //   printf("%f\t%f\t%f\t%f\n", th[array_index], ph[array_index], re_ylms.at(mode_index).at(array_index), im_ylms.at(mode_index).at(array_index));
    // }

    // // print 1 0 mode
    // printf("Ylms for 1 0 mode\n");
    // printf("th\tph\tre(ylm)\tim(ylm)\n");
    // mode_index = l_m_to_index(1, 0);
    // for(int array_index=0; array_index<array_size; array_index++){
    //   printf("%f\t%f\t%f\t%f\n", th[array_index], ph[array_index], re_ylms.at(mode_index).at(array_index), im_ylms.at(mode_index).at(array_index));
    // }

    // // print 1 1 mode
    // printf("Ylms for 1 1 mode\n");
    // printf("th\tph\tre(ylm)\tim(ylm)\n");
    // mode_index = l_m_to_index(1, 1);
    // for(int array_index=0; array_index<array_size; array_index++){
    //   printf("%f\t%f\t%f\t%f\n", th[array_index], ph[array_index], re_ylms.at(mode_index).at(array_index), im_ylms.at(mode_index).at(array_index));
    // }

    // // print 2 2 mode
    // printf("Ylms for 2 2 mode\n");
    // printf("th\tph\tre(ylm)\tim(ylm)\n");
    // mode_index = l_m_to_index(2, 2);
    // for(int array_index=0; array_index<array_size; array_index++){
    //   printf("%f\t%f\t%f\t%f\n", th[array_index], ph[array_index], re_ylms.at(mode_index).at(array_index), im_ylms.at(mode_index).at(array_index));
    // }

    // Decompose g, dr_g, dt_g
    // re_g[i][j][mode], im_g[i][j][mode]
    vector<vector<vector<CCTK_REAL>>> re_g(3, vector<vector<CCTK_REAL>>(3, vector<CCTK_REAL>(mode_count)));
    vector<vector<vector<CCTK_REAL>>> im_g(3, vector<vector<CCTK_REAL>>(3, vector<CCTK_REAL>(mode_count)));
    vector<vector<vector<CCTK_REAL>>> re_dr_g(3, vector<vector<CCTK_REAL>>(3, vector<CCTK_REAL>(mode_count)));
    vector<vector<vector<CCTK_REAL>>> im_dr_g(3, vector<vector<CCTK_REAL>>(3, vector<CCTK_REAL>(mode_count)));
    vector<vector<vector<CCTK_REAL>>> re_dt_g(3, vector<vector<CCTK_REAL>>(3, vector<CCTK_REAL>(mode_count)));
    vector<vector<vector<CCTK_REAL>>> im_dt_g(3, vector<vector<CCTK_REAL>>(3, vector<CCTK_REAL>(mode_count)));
    for(int i=0; i<3; i++){
      for(int j=i; j<3; j++){
	Decompose_Spherical_Harmonics(th, ph, g.at(i).at(j), re_g.at(i).at(j), im_g.at(i).at(j), re_ylms, im_ylms, array_size, lmax, ntheta, nphi);
	Decompose_Spherical_Harmonics(th, ph, dr_g.at(i).at(j), re_dr_g.at(i).at(j), im_dr_g.at(i).at(j), re_ylms, im_ylms, array_size, lmax, ntheta, nphi);
	Decompose_Spherical_Harmonics(th, ph, dt_g.at(i).at(j), re_dt_g.at(i).at(j), im_dt_g.at(i).at(j), re_ylms, im_ylms, array_size, lmax, ntheta, nphi);
      }
    }

    printf("Decomposed metric data\n");

    // Decompose beta, dr_beta, dt_beta
    // re_beta[i][mode]
    vector<vector<CCTK_REAL>> re_beta(3, vector<CCTK_REAL>(mode_count));
    vector<vector<CCTK_REAL>> im_beta(3,  vector<CCTK_REAL>(mode_count));
    vector<vector<CCTK_REAL>> re_dr_beta(3, vector<CCTK_REAL>(mode_count));
    vector<vector<CCTK_REAL>> im_dr_beta(3, vector<CCTK_REAL>(mode_count));
    vector<vector<CCTK_REAL>> re_dt_beta(3, vector<CCTK_REAL>(mode_count));
    vector<vector<CCTK_REAL>> im_dt_beta(3, vector<CCTK_REAL>(mode_count));
    for(int i=0; i<3; i++){
      Decompose_Spherical_Harmonics(th, ph, beta.at(i), re_beta.at(i), im_beta.at(i), re_ylms, im_ylms, array_size, lmax, ntheta, nphi);
      Decompose_Spherical_Harmonics(th, ph, dr_beta.at(i), re_dr_beta.at(i), im_dr_beta.at(i), re_ylms, im_ylms, array_size, lmax, ntheta, nphi);
      Decompose_Spherical_Harmonics(th, ph, dt_beta.at(i), re_dt_beta.at(i), im_dt_beta.at(i), re_ylms, im_ylms, array_size, lmax, ntheta, nphi);
    }

    printf("Decomposed shift data\n");

    // Decompose alpha, dr_alpha, dt_alpha
    vector<CCTK_REAL> re_alpha(mode_count);
    vector<CCTK_REAL> im_alpha(mode_count);
    vector<CCTK_REAL> re_dr_alpha(mode_count);
    vector<CCTK_REAL> im_dr_alpha(mode_count);
    vector<CCTK_REAL> re_dt_alpha(mode_count);
    vector<CCTK_REAL> im_dt_alpha(mode_count);
    Decompose_Spherical_Harmonics(th, ph, alpha, re_alpha, im_alpha, re_ylms, im_ylms, array_size, lmax, ntheta, nphi);
    Decompose_Spherical_Harmonics(th, ph, dr_alpha, re_dr_alpha, im_dr_alpha, re_ylms, im_ylms, array_size, lmax, ntheta, nphi);
    Decompose_Spherical_Harmonics(th, ph, dt_alpha, re_dt_alpha, im_dt_alpha, re_ylms, im_ylms, array_size, lmax, ntheta, nphi);

    printf("Decomposed lapse data\n");

    // Store output in h5 file
  }

}
