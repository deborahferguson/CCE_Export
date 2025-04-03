#include "cce_export.hh"
#include <vector>
//#define H5_USE_16_API
//#include <hdf5.h>
#include <sys/stat.h>
#include <iomanip>
#include <string.h>
#include <map>

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


#define HDF5_ERROR(fn_call)						\
  do {                                                                  \
  hid_t _error_code = fn_call;						\
  									\
  									\
  if (_error_code < 0)							\
    {                                                                   \
      CCTK_VWarn (0, __LINE__, __FILE__, CCTK_THORNSTRING,              \
		  "HDF5 call '%s' returned error code %d",              \
		  #fn_call, (int)_error_code);                          \
    }                                                                   \
  } while (0)


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

#ifdef HAVE_CAPABILITY_HDF5

static bool file_exists(const string &name)
{
  struct stat sts;
  return !(stat(name.c_str(), &sts) == -1 && errno == ENOENT);
}

static bool dataset_exists(hid_t file, const string &dataset_name)
{
  // To test whether a dataset exists, the recommended way in API 1.6
  // is to use H5Gget_objinfo, but this prints an error to stderr if
  // the dataset does not exist.  We explicitly avoid this by wrapping
  // the call in H5E_BEGIN_TRY/H5E_END_TRY statements.  In 1.8,
  // H5Gget_objinfo is deprecated, and H5Lexists does the job.  See
  // http://www.mail-archive.com/hdf-forum@hdfgroup.org/msg00125.html

  #if 1
  bool exists;
  H5E_BEGIN_TRY
    {
      exists = H5Gget_objinfo(file, dataset_name.c_str(), 1, NULL) >= 0;
    }
  H5E_END_TRY;
  return exists;
  #else
  return H5Lexists(file, dataset_name.c_str(), H5P_DEFAULT);
  #endif
}

void Create_Dataset(CCTK_ARGUMENTS, hid_t file, string datasetname, CCTK_REAL *data, int lmax)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  int mode_count = l_m_to_index(lmax, lmax)+1;

  hid_t dataset = -1;
  printf(datasetname.c_str());
  printf("\n");
  
  if (dataset_exists(file, datasetname))
    {
    dataset = H5Dopen(file, datasetname.c_str());
  }
  else
  {
    hsize_t dims[2] = {0, (hsize_t)(2*mode_count + 1)};
    hsize_t maxdims[2] = {H5S_UNLIMITED, (hsize_t)(2*mode_count + 1)};
    hid_t dataspace = H5Screate_simple(2, dims, maxdims);
    
    hid_t cparms = -1;
    hsize_t chunk_dims[2] = {hsize_t(hdf5_chunk_size), (hsize_t)(2*mode_count + 1)};
    cparms = H5Pcreate(H5P_DATASET_CREATE);
    printf("Property list handle: %lld\n", (long long)cparms);
    HDF5_ERROR(H5Pset_chunk(cparms, 2, chunk_dims));
    
    dataset = H5Dcreate(file, datasetname.c_str(), H5T_NATIVE_DOUBLE, dataspace, cparms);
    H5Pclose(cparms);

    // Create the legend 

    // Create variable-length string datatype
    hid_t str_type = H5Tcopy(H5T_C_S1);
    printf("string type: %lld\n", (long long)str_type);
    H5Tset_size(str_type, H5T_VARIABLE);

    // Create dataspace for array of strings
    hsize_t legend_dims[] = {(hsize_t)(2*mode_count + 1)}; 
    hid_t space = H5Screate_simple(1, legend_dims, NULL);
    printf("space: %lld\n", (long long)space);

    // Create attribute
    hid_t attr = H5Acreate(dataset, "Legend", str_type, space, H5P_DEFAULT);
    printf("attr: %lld\n", (long long)attr);

    // Write array of strings
    char* legend[2*mode_count + 1];
    legend[0] = "time";
    for(int l=0; l<lmax+1; l++){
      for(int m=-l; m<l+1; m++){
    	int mode_index = l_m_to_index(l, m);
    	ostringstream re_label;
    	re_label << "Re(" << l << "," << m <<")";
	printf("%s\n", re_label.str().c_str());
    	legend[2*mode_index+1] = strdup(re_label.str().c_str());
    	ostringstream im_label;
    	im_label << "Im(" << l << "," << m <<")";
	printf("%s\n", im_label.str().c_str());
    	legend[2*mode_index+2] = strdup(im_label.str().c_str());
      }
    }
    printf("Made the legend array\n");
    // for(int i = 0; i < 2*mode_count+1; i++) {
    //   printf("%s\n", legend[i]);
    // }

    H5Awrite(attr, str_type, legend);
  }
  
  printf("About to try to write to dataset\n");
  
  hid_t filespace = H5Dget_space (dataset);
  
  printf("%ld\n", (long)filespace);  
  
  hsize_t filedims[2];
  hsize_t maxdims[2];
  HDF5_ERROR(H5Sget_simple_extent_dims(filespace, filedims, maxdims));
  
  filedims[0] += 1;
  hsize_t size[2] = {filedims[0], filedims[1]};
  HDF5_ERROR(H5Dextend (dataset, size));
  HDF5_ERROR(H5Sclose(filespace));
  
  /* Select a hyperslab  */
  hsize_t offset[2] = {filedims[0]-1, 0};
  // hsize_t dims2[2] = {1, 3};
  hsize_t dims2[2] = {1, (hsize_t)(2*mode_count + 1)};
  filespace = H5Dget_space (dataset);
  HDF5_ERROR(H5Sselect_hyperslab (filespace, H5S_SELECT_SET, offset, NULL,
				  dims2, NULL));
  
  hid_t memdataspace = H5Screate_simple(2, dims2, NULL);
  
  /* Write the data to the hyperslab  */
  HDF5_ERROR(H5Dwrite (dataset, H5T_NATIVE_DOUBLE, memdataspace, filespace,
		       H5P_DEFAULT, data));
  
  HDF5_ERROR(H5Dclose(dataset));
  HDF5_ERROR(H5Sclose(filespace));
  HDF5_ERROR(H5Sclose(memdataspace));
}
	      
void Output_Decomposed_Metric_Data(CCTK_ARGUMENTS, vector<vector<vector<CCTK_REAL>>> &re_g, vector<vector<vector<CCTK_REAL>>> &im_g, 
				   vector<vector<vector<CCTK_REAL>>> &re_dr_g, vector<vector<vector<CCTK_REAL>>> &im_dr_g, 
				   vector<vector<vector<CCTK_REAL>>> &re_dt_g, vector<vector<vector<CCTK_REAL>>> &im_dt_g, 
				   vector<vector<CCTK_REAL>> &re_beta, vector<vector<CCTK_REAL>> &im_beta, 
				   vector<vector<CCTK_REAL>> &re_dr_beta, vector<vector<CCTK_REAL>> &im_dr_beta, 
				   vector<vector<CCTK_REAL>> &re_dt_beta, vector<vector<CCTK_REAL>> &im_dt_beta, 
				   vector<CCTK_REAL> &re_alpha, vector<CCTK_REAL> &im_alpha, 
				   vector<CCTK_REAL> &re_dr_alpha, vector<CCTK_REAL> &im_dr_alpha, 
				   vector<CCTK_REAL> &re_dt_alpha, vector<CCTK_REAL> &im_dt_alpha, 
				   float rad, int lmax)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  int mode_count = l_m_to_index(lmax, lmax) + 1;

  const char *my_out_dir = strcmp(out_dir, "") ? out_dir : io_out_dir;
  if (CCTK_CreateDirectory(0755, my_out_dir) < 0)
    CCTK_VWarn(CCTK_WARN_ABORT, __LINE__, __FILE__, CCTK_THORNSTRING,
               "Multipole output directory %s could not be created",
               my_out_dir);

  static map<string,bool> checked; 
  // static bool checked; // Has the given file been checked
                                   // for truncation? bool
                                   // defaults to false  

  ostringstream basename;
  basename << "CCE_Export_R" << setiosflags(ios::fixed) << setprecision(2) << rad << ".h5";
  string output_name = my_out_dir + string("/") + basename.str();

  printf(output_name.c_str());
  printf("\n");

  hid_t file;

  if (!file_exists(output_name) || (!checked[output_name] && IO_TruncateOutputFiles(cctkGH)))
    {
      printf("File does not exist\n");
      file = H5Fcreate(output_name.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    }
  else
    {
      printf("File exists\n");
      file = H5Fopen(output_name.c_str(), H5F_ACC_RDWR, H5P_DEFAULT);
    }

  checked[output_name] = true;

  // store metric data
  for(int i = 0; i<3; i++){
    for(int j=i; j<3; j++){
      string components;
      if(i==0 && j==0){components = "xx";}
      if(i==0 && j==1){components = "xy";}
      if(i==0 && j==2){components = "xz";}
      if(i==1 && j==1){components = "yy";}
      if(i==1 && j==2){components = "yz";}
      if(i==2 && j==2){components = "zz";}

      string datasetname = "g" + components + ".dat";
      string dt_datasetname = "Dtg" + components + ".dat";
      string dr_datasetname = "Drg" + components + ".dat";

      CCTK_REAL data[2*mode_count + 1];
      CCTK_REAL dt_data[2*mode_count + 1];
      CCTK_REAL dr_data[2*mode_count + 1];
      data[0] = cctk_time;
      dt_data[0] = cctk_time;
      dr_data[0] = cctk_time;
      for(int l=0; l<=lmax; l++){
      	for(int m=-l; m<l+1; m++){
      	  int mode_index = l_m_to_index(l, m);
      	  // printf("%f, %f, %f\n", l, m, mode_index);
	  data[2*mode_index+1] = re_g[i][j][mode_index];
	  data[2*mode_index+2] = im_g[i][j][mode_index];
	  dt_data[2*mode_index+1] = re_dt_g[i][j][mode_index];
	  dt_data[2*mode_index+2] = im_dt_g[i][j][mode_index];
	  dr_data[2*mode_index+1] = re_dr_g[i][j][mode_index];
	  dr_data[2*mode_index+2] = im_dr_g[i][j][mode_index];
	}
      }

      Create_Dataset(CCTK_PASS_CTOC, file, datasetname, data, lmax);
      Create_Dataset(CCTK_PASS_CTOC, file, dt_datasetname, dt_data, lmax);
      Create_Dataset(CCTK_PASS_CTOC, file, dr_datasetname, dr_data, lmax);

    }
  }

  // store shift data
  for(int i=0; i<3; i++){
    string component;
    if(i==0){component = "x";}
    if(i==1){component = "y";}
    if(i==2){component = "z";}
    
    string datasetname = "Shift" + component + ".dat";
    string dt_datasetname = "DtShift" + component + ".dat";
    string dr_datasetname = "DrShift" + component + ".dat";

    CCTK_REAL data[2*mode_count + 1];
    CCTK_REAL dt_data[2*mode_count + 1];
    CCTK_REAL dr_data[2*mode_count + 1];
    data[0] = cctk_time;
    dt_data[0] = cctk_time;
    dr_data[0] = cctk_time;
    //for(int mode_index = 0; mode_index<mode_count; mode_index++){
    for(int l=0; l<=lmax; l++){
      for(int m=-l; m<l+1; m++){
    	int mode_index = l_m_to_index(l, m);
	data[2*mode_index+1] = re_beta[i][mode_index];
	data[2*mode_index+2] = im_beta[i][mode_index];
	dt_data[2*mode_index+1] = re_dt_beta[i][mode_index];
	dt_data[2*mode_index+2] = im_dt_beta[i][mode_index];
	dr_data[2*mode_index+1] = re_dr_beta[i][mode_index];
	dr_data[2*mode_index+2] = im_dr_beta[i][mode_index];
      }
    } 
	  
    Create_Dataset(CCTK_PASS_CTOC, file, datasetname, data, lmax);
    Create_Dataset(CCTK_PASS_CTOC, file, dt_datasetname, dt_data, lmax);
    Create_Dataset(CCTK_PASS_CTOC, file, dr_datasetname, dr_data, lmax);
  }

  // store lapse data
  string datasetname = "Lapse.dat";
  string dt_datasetname = "DtLapse.dat";
  string dr_datasetname = "DrLapse.dat";

  CCTK_REAL data[2*mode_count + 1];
  CCTK_REAL dt_data[2*mode_count + 1];
  CCTK_REAL dr_data[2*mode_count + 1];
  data[0] = cctk_time;
  dt_data[0] = cctk_time;
  dr_data[0] = cctk_time;
  for(int l=0; l<=lmax; l++){
    for(int m=-l; m<l+1; m++){
      int mode_index = l_m_to_index(l, m);
      data[2*mode_index+1] = re_alpha[mode_index];
      data[2*mode_index+2] = im_alpha[mode_index];
      dt_data[2*mode_index+1] = re_dt_alpha[mode_index];
      dt_data[2*mode_index+2] = im_dt_alpha[mode_index];
      dr_data[2*mode_index+1] = re_dr_alpha[mode_index];
      dr_data[2*mode_index+2] = im_dr_alpha[mode_index];
    }
  }

  Create_Dataset(CCTK_PASS_CTOC, file, datasetname, data, lmax);
  Create_Dataset(CCTK_PASS_CTOC, file, dt_datasetname, dt_data, lmax);
  Create_Dataset(CCTK_PASS_CTOC, file, dr_datasetname, dr_data, lmax);  

}

#else

void  Output_Decomposed_Metric_Data(CCTK_ARGUMENTS, vector<vector<vector<CCTK_REAL>>> &re_g, vector<vector<vector<CCTK_REAL>>> &im_g,
				    vector<vector<vector<CCTK_REAL>>> &re_dr_g, vector<vector<vector<CCTK_REAL>>> &im_dr_g,
				    vector<vector<vector<CCTK_REAL>>> &re_dt_g, vector<vector<vector<CCTK_REAL>>> &im_dt_g,
				    vector<vector<CCTK_REAL>> &re_beta, vector<vector<CCTK_REAL>> &im_beta,
				    vector<vector<CCTK_REAL>> &re_dr_beta, vector<vector<CCTK_REAL>> &im_dr_beta,
				    vector<vector<CCTK_REAL>> &re_dt_beta, vector<vector<CCTK_REAL>> &im_dt_beta,
				    vector<CCTK_REAL> &re_alpha, vector<CCTK_REAL> &im_alpha,
				    vector<CCTK_REAL> &re_dr_alpha, vector<CCTK_REAL> &im_dr_alpha,
				    vector<CCTK_REAL> &re_dt_alpha, vector<CCTK_REAL> &im_dt_alpha,
				    float rad, int lmax)
{
  CCTK_WARN(0,"HDF5 output has been requested but Cactus has been compiled without HDF5 support");
}

#endif

void CCE_Export(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  printf("CCE_Export: starting\n");

  static string index_to_component[] = {"x", "y", "z"};

  //const int ntheta = 50; 
  //const int nphi = 100;
  const int ntheta = 120; 
  const int nphi = 240;
  const int array_size = (ntheta+1)*(nphi+1);
  //const int array_size=ntheta*nphi;

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

  for (int theta_index = 0; theta_index <= ntheta; theta_index++)
  {
    for (int phi_index = 0; phi_index <= nphi; phi_index++)
    {
      const int array_index = theta_index + (ntheta+1) * phi_index;
	
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
    for (int theta_index = 0; theta_index <= ntheta; theta_index++)
    {
      for (int phi_index = 0; phi_index <= nphi; phi_index++)
      {
	const int array_index = theta_index + (ntheta+1) * phi_index;
	
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

      // dt_g.at(1).at(0).at(array_index) = dt_g.at(0).at(1).at(array_index);

      // dt_g.at(2).at(0).at(array_index) = dt_g.at(0).at(2).at(array_index);

      // dt_g.at(2).at(1).at(array_index) = dt_g.at(1).at(2).at(array_index);
    
    }
    
    // // print the values of gxx on the sphere
    // printf("CCE_Export: theta\tphi\tgxx\n");
    // for(int array_index=0; array_index<array_size; array_index++){
    //   printf("%f\t%f\t%f\n", th[array_index], ph[array_index], g.at(0).at(0).at(array_index));
    // }

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

    // printf("CCE_Export: l\tm\tre_gxx_lm\tim_gxx_lm\n");
    // for(int l=0; l<lmax+1; l++){
    //   for(int m=-l; m<l+1; m++){
    // 	int mode_index = l_m_to_index(l, m);
    // 	printf("%d\t%d\t%f\t%f\n", l, m, re_g.at(0).at(0).at(mode_index), im_g.at(0).at(0).at(mode_index));
    //   }
    // }

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
    if (CCTK_MyProc(cctkGH) == 0)
    {
      Output_Decomposed_Metric_Data(CCTK_PASS_CTOC, re_g, im_g, re_dr_g, im_dr_g, re_dt_g, im_dt_g, 
				    re_beta, im_beta, re_dr_beta, im_dr_beta, re_dt_beta, im_dt_beta, 
				    re_alpha, im_alpha, re_dr_alpha, im_dr_alpha, re_dt_alpha, im_dt_alpha, radius[r], lmax);
    }
  }

}
