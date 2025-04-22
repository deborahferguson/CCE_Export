#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "hdf5.h"
#ifndef TRUE
#define TRUE 1
#define FALSE 0
#endif
#define CHECK_ERROR(err) \
  if (err < 0) { \
  fprintf(stderr, "Error at line %d\n", __LINE__); \
  return 1; \
  }

int main(int argc, char *argv[]) {
  if (argc != 2) {
    fprintf(stderr, "Usage: %s <input_h5_file>\n", argv[0]);
    return 1;
  }

  // Initialize HDF5 library
  CHECK_ERROR(H5dont_atexit());
  hid_t file_id = -1, group_id = -1, dataset_id = -1;
  herr_t status;

  // Open file
  file_id = H5Fopen(argv[1], H5F_ACC_RDONLY, H5P_DEFAULT);
  CHECK_ERROR(file_id);

  // Open root group
  group_id = H5Gopen2(file_id, "/", H5P_DEFAULT);
  CHECK_ERROR(group_id);

  hsize_t num_objects;
  status = H5Gget_num_objs(group_id, &num_objects);
  CHECK_ERROR(status);

  for (hsize_t i = 0; i < num_objects; i++) {
    char obj_name[256];
    H5O_info_t obj_info;

    // Get object name
    status = H5Lget_name_by_idx(group_id, ".", H5_INDEX_NAME, H5_ITER_INC, i,
				obj_name, sizeof(obj_name), H5P_DEFAULT);
    CHECK_ERROR(status);

    // Print dataset name
    printf("Dataset: %s\n", obj_name);

    // Try to open as dataset
    dataset_id = H5Dopen2(group_id, obj_name, H5P_DEFAULT);
    if (dataset_id >= 0) {
      hid_t dataspace_id = -1, datatype_id = -1;
      hsize_t dims[2];
      double *data = NULL;

      // Get dataspace
      dataspace_id = H5Dget_space(dataset_id);
      CHECK_ERROR(dataspace_id);

      // Get dimensions
      status = H5Sget_simple_extent_dims(dataspace_id, dims, NULL);
      CHECK_ERROR(status);

      // Get datatype
      datatype_id = H5Dget_type(dataset_id);
      CHECK_ERROR(datatype_id);

      // Check if it's double precision
      H5T_class_t class_id = H5Tget_class(datatype_id);
      if (class_id == H5T_FLOAT) {
	// Allocate memory for double precision data
	size_t total_size = sizeof(double) * dims[0] * dims[1];
	data = (double*)malloc(total_size);
	if (!data) {
	  fprintf(stderr, "Memory allocation failed\n");
	  goto cleanup_dataset;
	}

	// Read data
	status = H5Dread(dataset_id, datatype_id, H5S_ALL, H5S_ALL,
			 H5P_DEFAULT, data);
	CHECK_ERROR(status);

	// Print data to stdout
	for (hsize_t row = 0; row < dims[0]; row++) {
	  for (hsize_t col = 0; col < dims[1]; col++) {
	    printf("%.6e ", data[row * dims[1] + col]);
	  }
	  printf("\n");
	}
	free(data);
      }
    cleanup_dataset:
      if (datatype_id >= 0) H5Tclose(datatype_id);
      if (dataspace_id >= 0) H5Sclose(dataspace_id);
      if (dataset_id >= 0) H5Dclose(dataset_id);
    }
  }

  // Close remaining resources
  if (group_id >= 0) H5Gclose(group_id);
  if (file_id >= 0) H5Fclose(file_id);
  return 0;
}
