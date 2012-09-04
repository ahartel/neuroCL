#include "errorMessage.h"
#include <iostream>
#include <fstream>
#include <cassert>

using namespace std;

/*
 * Constants
 */
const int N = 1000;
const float v_thresh = 30/*mV*/;

void get_device_info(cl_device_id* device)
{
	/*
	 * get device info
	 */
	cl_int error = 0;   // Used to handle error codes

	size_t param_value_size_ret;
	// query device name
	char device_name[100];
	error = clGetDeviceInfo(*device,CL_DEVICE_NAME,100,&device_name,&param_value_size_ret);
	check_error("getting device name",error);
	assert (error == CL_SUCCESS);
	std::cout << "Device name: " << device_name << std::endl;
	// openCL version
	char opencl_version[100];
	error = clGetDeviceInfo(*device,CL_DRIVER_VERSION,100,&opencl_version,&param_value_size_ret);
	check_error("getting openCL version",error);
	assert (error == CL_SUCCESS);
	std::cout << "openCL Version: " << opencl_version << std::endl;
	// max. compute units
	cl_uint max_comp_units;
	error = clGetDeviceInfo(*device,CL_DEVICE_MAX_COMPUTE_UNITS,sizeof(cl_uint),&max_comp_units,&param_value_size_ret);
	check_error("getting max. compute units",error);
	assert (error == CL_SUCCESS);
	std::cout << "max. compute units: " << max_comp_units << std::endl;
	// query work_item dimensions
	cl_uint max_work_item_dimensions;
	error = clGetDeviceInfo(*device,CL_DEVICE_MAX_WORK_ITEM_DIMENSIONS,sizeof(cl_uint),&max_work_item_dimensions,NULL);
	check_error("getting device information",error);
	assert (error == CL_SUCCESS);
	std::cout << "max. work_item dimensions: " << max_work_item_dimensions << std::endl;
	// query max. work_group sizes
	size_t max_work_group_size;
	error = clGetDeviceInfo(*device,CL_DEVICE_MAX_WORK_GROUP_SIZE,sizeof(size_t),&max_work_group_size,NULL);
	check_error("getting max. work group size",error);
	assert (error == CL_SUCCESS);
	std::cout << "max. work_group size: " << max_work_group_size << std::endl;
	// query work_item sizes
	size_t max_work_item_sizes[max_work_item_dimensions];
	clGetDeviceInfo(*device,CL_DEVICE_MAX_WORK_ITEM_SIZES,max_work_item_dimensions*sizeof(size_t),&max_work_item_sizes,NULL);
	check_error("getting device information",error);
	assert (error == CL_SUCCESS);
	for (int i = 0; i < max_work_item_dimensions; i++)
		std::cout << "max. work_item size for dim. " << i << ": " << max_work_item_sizes[i] << std::endl;
	// max memory allocation size
	unsigned long long max_mem_alloc_size;
	clGetDeviceInfo(*device,CL_DEVICE_MAX_MEM_ALLOC_SIZE,sizeof(unsigned long long),&max_mem_alloc_size,NULL);
	check_error("getting max. mem allocation size",error);
	assert (error == CL_SUCCESS);
	cout << "Max. memory allocation size: " << max_mem_alloc_size/1024/1024 << "M" << endl;
}

void init_neurons(float* membranes)
{
	for (int i=0; i<N; i++) {
		membranes[i] = (float)rand()/float(RAND_MAX)*50.0;
	}
}

int main() {
	cl_int error = 0;   // Used to handle error codes
	cl_platform_id platform;
	cl_context context;
	cl_command_queue queue;
	cl_device_id device;

	/*
	 * set up system
	 */
	// Platform
	error = clGetPlatformIDs( 1, &platform, NULL );
	if (error != CL_SUCCESS) {
	   cout << "Error getting platform id: " << errorMessage(error) << endl;
	   exit(error);
	}
	// Device
	error = clGetDeviceIDs(platform, CL_DEVICE_TYPE_ALL, 1, &device, NULL);
	if (error != CL_SUCCESS) {
	   cout << "Error getting device ids: " << errorMessage(error) << endl;
	   exit(error);
	}
	// Context
	context = clCreateContext(0, 1, &device, NULL, NULL, &error);
	if (error != CL_SUCCESS) {
	   cout << "Error creating context: " << errorMessage(error) << endl;
	   exit(error);
	}
	// Command-queue
	queue = clCreateCommandQueue(context, device, 0, &error);
	if (error != CL_SUCCESS) {
	   cout << "Error creating command queue: " << errorMessage(error) << endl;
	   exit(error);
	}

	// print some information about the device to stdout
	get_device_info(&device);

	/*
	 * load and compile kernel
	 */
	std::ifstream file("source/kernels/neuron_fired.c");
	std::string source;
	while(!file.eof()){
		char line[256];
		file.getline(line,255);
		source += line;
	}
	const char* str = source.c_str();
	cl_program program = clCreateProgramWithSource( context,
                         1,
                         &str,
                         NULL, NULL );
	cl_int result = clBuildProgram( program, 1, &device, NULL, NULL, NULL );
	if ( result )
	{
		std::cout << "Error during compilation! (" << result << ")" << std::endl;

		// Shows the log
		char* build_log;
		size_t log_size;
		// First call to know the proper size
		clGetProgramBuildInfo(program, device, CL_PROGRAM_BUILD_LOG, 0, NULL, &log_size);
		build_log = new char[log_size+1];
		// Second call to get the log
		clGetProgramBuildInfo(program, device, CL_PROGRAM_BUILD_LOG,
			log_size, build_log, NULL);

		build_log[log_size] = '\0';
		cout << build_log << endl;
		delete[] build_log;
	}

	// Extracting the kernel
	cl_kernel neuron_fired = clCreateKernel(program, "neuron_fired", &error);
	assert(error == CL_SUCCESS);

	/*
	 * What we need here:
	 * Data structures
	 * - list of neurons of size N
	 * - list of spikes of size N
	 * Kernel
	 * - iterate over all neurons, check v against threshold
	 * - add event to spike buffer
	 */
	size_t local_ws,global_ws,num_work_groups;
	local_ws = 64;
	global_ws = 1024;
	num_work_groups = 16;

	float* membranes = new float[N];
	unsigned int* spikes = new unsigned int[N];
	unsigned int* k = new unsigned int[1];
	k[0] = 0;

	init_neurons(membranes);

	cl_mem cl_membranes = clCreateBuffer(context, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR, sizeof(float)*N, membranes, &error);
	check_error("setting kernel args",error);
	assert(error == CL_SUCCESS);
	cl_mem cl_spikes = clCreateBuffer(context, CL_MEM_WRITE_ONLY, sizeof(unsigned int)*N, NULL, &error);
	check_error("setting kernel args",error);
	assert(error == CL_SUCCESS);
	cl_mem cl_k = clCreateBuffer(context, CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR, sizeof(unsigned int), k, &error);
	check_error("setting kernel args",error);
	assert(error == CL_SUCCESS);

	// Enqueuing parameters
	// Note that we inform the size of the cl_mem object, not the size of the memory pointed by it
	error = clSetKernelArg(neuron_fired, 0, sizeof(cl_mem), (void*)&cl_membranes);
	error |= clSetKernelArg(neuron_fired, 1, sizeof(cl_mem), (void*)&cl_spikes);
	error |= clSetKernelArg(neuron_fired, 2, sizeof(cl_mem), (void*)&cl_k);
	error |= clSetKernelArg(neuron_fired, 3, sizeof(cl_float), (void*)&v_thresh);
	error |= clSetKernelArg(neuron_fired, 4, sizeof(cl_uint), (void*)&N);
	check_error("setting kernel args",error);
	assert(error == CL_SUCCESS);

	error = clEnqueueNDRangeKernel(queue, neuron_fired, 1, NULL, &global_ws, &local_ws, 0, NULL, NULL);
	check_error("enqueueing NDRangeKernel",error);
	assert(error == CL_SUCCESS);

	// Reading back
	unsigned int* check_k;
	error = clEnqueueReadBuffer(queue, cl_spikes, CL_TRUE, 0, sizeof(unsigned int)*N, spikes, 0, NULL, NULL);
	error = clEnqueueReadBuffer(queue, cl_k, CL_TRUE, 0, sizeof(unsigned int), check_k, 0, NULL, NULL);
	check_error("enqueueing read buffer",error);
	assert(error == CL_SUCCESS);

	cout << "Results:" << endl;
	unsigned int num_fired=0;
	for (int i=0; i<N; i++) {
		if (membranes[i] > v_thresh) {
			num_fired++;
		}
	}
	cout << "expected: " << num_fired << endl;
	cout << "found: " << check_k[0] << endl;

	delete[] k;
	delete[] spikes;
	delete[] membranes;
	clReleaseKernel(neuron_fired);
	clReleaseCommandQueue(queue);
	clReleaseContext(context);
	clReleaseMemObject(cl_k);
	clReleaseMemObject(cl_spikes);
	clReleaseMemObject(cl_membranes);
}
