#include "errorMessage.h"
#include <iostream>
#include <fstream>
#include <cassert>

using namespace std;

#define MAX_NUM_PLATFORMS 10
#define MAX_NUM_DEVICES 10

#include "helpers.h"

void get_platform_info(cl_platform_id* platform)
{
	/*
	 * get device info
	 */
	cl_int error = 0;   // Used to handle error codes

	size_t param_value_size_ret;
	// query device name
	char platform_profile[100];
	error = clGetPlatformInfo(*platform,CL_PLATFORM_PROFILE,100,&platform_profile,&param_value_size_ret);
	check_error("getting platform profile",error);
	assert (error == CL_SUCCESS);
	std::cout << "Platform profile: " << platform_profile << std::endl;
	// query platform extensions
	/*
	char platform_extensions[512];
	error = clGetPlatformInfo(*platform,CL_PLATFORM_EXTENSIONS,512,&platform_extensions,&param_value_size_ret);
	check_error("getting platform extensions",error);
	assert (error == CL_SUCCESS);
	std::cout << "Platform extenions: " << platform_extensions << std::endl;
	*/
}

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
	// Does the device support cl_khr_fp64?
	cl_uint double_vector_size;
	error = clGetDeviceInfo(*device,CL_DEVICE_PREFERRED_VECTOR_WIDTH_DOUBLE,sizeof(cl_uint),&double_vector_size,NULL);
	check_error("getting device information",error);
	assert (error == CL_SUCCESS);
	std::cout << "preferred vector width for double: " << double_vector_size << std::endl;
	// Device extensions
	/*
	char extensions[512];
	error = clGetDeviceInfo(*device,CL_DEVICE_EXTENSIONS,512,&extensions,&param_value_size_ret);
	check_error("getting device information",error);
	assert (error == CL_SUCCESS);
	std::cout << "Device extensions: " << extensions << std::endl;
	*/
}

void setup(cl_platform_id* platforms,
	cl_context* context,
	cl_command_queue* queue,
	cl_device_id* devices)
{
	cl_int error = 0;   // Used to handle error codes
	/*
	 * set up system
	 */
	// Platform
	cl_uint num_platforms;
	error = clGetPlatformIDs(MAX_NUM_PLATFORMS, platforms, &num_platforms );
	if (error != CL_SUCCESS) {
	   cout << "Error getting platform id: " << errorMessage(error) << endl;
	   exit(error);
	}
	std::cout << "Number of platforms found: " << num_platforms << std::endl;
	std::cout << "Using platform 0" << std::endl;
	// print some information about the platform to stdout
	get_platform_info(&(platforms[0]));
	cout << "-----------------------------------" << endl;
	// Device
	cl_uint num_devices;
	error = clGetDeviceIDs(platforms[0], CL_DEVICE_TYPE_ALL, MAX_NUM_DEVICES, devices, &num_devices);
	if (error != CL_SUCCESS) {
	   cout << "Error getting device ids: " << errorMessage(error) << endl;
	   exit(error);
	}
	std::cout << "Number of devices found: " << num_devices << std::endl;
	std::cout << "Using device 0" << std::endl;
	// Context
	*context = clCreateContext(0, 1, &(devices[0]), NULL, NULL, &error);
	if (error != CL_SUCCESS) {
	   cout << "Error creating context: " << errorMessage(error) << endl;
	   exit(error);
	}
	// Command-queue
	*queue = clCreateCommandQueue(*context, devices[0], 0, &error);
	if (error != CL_SUCCESS) {
	   cout << "Error creating command queue: " << errorMessage(error) << endl;
	   exit(error);
	}

	// print some information about the device to stdout
	get_device_info(&(devices[0]));
	cout << "-----------------------------------" << endl;

}

void load_and_compile_kernel(cl_context* context,
	cl_device_id* device,
	const char* filename,
	const char* name,
	cl_kernel* kernel)
{
	cl_int error = 0;   // Used to handle error codes

	/*
	 * load and compile kernel
	 * First step: load source code from <filename>
	 */
	std::ifstream file(filename);
	assert (file.good());
	// read file contents
	std::string source;
	while(file.good()){
		char line[256];
		file.getline(line,255);
		source += line;
		source += string("\n");
	}
	file.close();
	const char* str = source.c_str();

	cl_program program = clCreateProgramWithSource( *context,
                         1,
                         &str,
                         NULL, &error );
	check_error("Compiling kernel "+string(name),error);
	assert(error == CL_SUCCESS);

	//cl_int result = clBuildProgram( program, 1, device, NULL, NULL, NULL );
	cl_int result = clBuildProgram( program, 0, NULL, NULL, NULL, NULL );
	if ( result )
	{
		std::cout << "Error during compilation! (" << result << ")" << std::endl;
	}
	else
	{
		std::cout << "Program in " << filename << " successfully compiled" << endl;
	}

	// Shows the log
	char* build_log;
	size_t log_size;
	// First call to know the proper size
	clGetProgramBuildInfo(program, *device, CL_PROGRAM_BUILD_LOG, 0, NULL, &log_size);
	build_log = new char[log_size+1];
	// Second call to get the log
	clGetProgramBuildInfo(program, *device, CL_PROGRAM_BUILD_LOG,
		log_size, build_log, NULL);

	build_log[log_size] = '\0';
	cout << build_log << endl;
	delete[] build_log;

	// Extracting the kernel
	*kernel = clCreateKernel(program, name, &error);
	check_error("Creating kernel "+string(name)+" from source file "+string(filename),error);
	assert(error == CL_SUCCESS);
}

void load_compiled_kernel(cl_context* context,
	cl_device_id* device,
	const char* filename,
	const char* name,
	cl_kernel* kernel)
{
	cl_int error = 0;   // Used to handle error codes
	cl_int binary_status[1];

	std::ifstream file(filename);
	assert (file.good());
	// read file contents
	unsigned char** binaries = new unsigned char*[1];

	char read_bytes[4096];
	file.read(read_bytes,4096);
	binaries[0] = (unsigned char*)read_bytes;

	size_t lengths[1];
	lengths[0] = file.gcount();

	const unsigned char** test = (const unsigned char**)binaries;

	for (unsigned int i=0; i<lengths[0]; i++)
		std::cout << test[0][i];
    std::cout << std::endl;

	file.close();

	cl_program program;

	program = clCreateProgramWithBinary( *context,
                         1,
                         device,
						 lengths,
						 test,
                         binary_status, &error );

	//delete[] binaries[0];
	delete[] binaries;
	cout << binary_status << endl;
	check_error("clCreateProgramWithBinary: binary_status",binary_status[0]);
	assert(binary_status == CL_SUCCESS);
	check_error("clCreateProgramWithBinary",error);
	assert(error == CL_SUCCESS);

	// Extracting the kernel
	*kernel = clCreateKernel(program, name, &error);
	check_error("clCreateKernel",error);
	assert(error == CL_SUCCESS);
}

int main()
{
	INIT_TIMER(complete)
	START_TIMER(complete)

	cl_int error = 0;   // Used to handle error codes
	cl_platform_id platforms[MAX_NUM_PLATFORMS];
	cl_context context;
	cl_command_queue queue;
	cl_device_id devices[MAX_NUM_DEVICES];

	setup(platforms,&context,&queue,devices);

	cl_kernel neuron_fired;
	load_and_compile_kernel(&context,&(devices[0]),"source/kernels/neuron_fired.c","neuron_fired",&neuron_fired);
	//load_compiled_kernel(&context,&(devices[0]),"source/kernels/neuron_fired.ir","neuron_fired",&neuron_fired);
	cl_kernel evolve_neuron;
	load_and_compile_kernel(&context,&(devices[0]),"source/kernels/evolve_neuron.c","evolve_neuron",&evolve_neuron);
	//load_compiled_kernel(&context,&device,"source/kernels/evolve_neuron.ir","evolve_neuron",&evolve_neuron);

	/*
	 Data structures
	  - list of neurons of size N
	  - list of spikes of size N
	*/
	size_t local_ws,global_ws;
	local_ws = 512;
	global_ws = 512*(int(N/512)+1);

	unsigned int memory_used = 0;
	float membranes[N];
		memory_used += N*sizeof(float);
	float u[N];
		memory_used += N*sizeof(float);
	float d[N];
		memory_used += N*sizeof(float);
	float a[N];
		memory_used += N*sizeof(float);
	float I[N][D];
		memory_used += D*N*sizeof(float);

	// network parameters
	unsigned int delay_start[N][D];
		memory_used += N*D*sizeof(unsigned int);
	unsigned int delay_count[N][D];
		memory_used += N*D*sizeof(unsigned int);
	vector<float> weights[N];
		memory_used += N*D*sizeof(float);
	vector<unsigned int> post_neurons[N];
		memory_used += N*D*sizeof(unsigned int);
	unsigned int num_post[N];
		memory_used += N*sizeof(unsigned int);

	size_t spike_array_length = int(N*T*h/10);
	unsigned int spikes[spike_array_length];
		memory_used += spike_array_length*sizeof(unsigned int);
	unsigned int k[T];
		memory_used += T*sizeof(unsigned int);
	k[0] = 0;

	cout << "Max. memory usage: " << memory_used/1024 << "kB" << endl;

	init_neurons_sparse(membranes,u,d,a,I,weights,delay_start,delay_count,post_neurons,num_post);

	return 0;

	cl_mem cl_membranes = clCreateBuffer(context, CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR, sizeof(float)*N, membranes, &error);
	check_error("setting kernel args",error);
	assert(error == CL_SUCCESS);
	cl_mem cl_u = clCreateBuffer(context, CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR, sizeof(float)*N, u, &error);
	check_error("setting kernel args",error);
	assert(error == CL_SUCCESS);
	cl_mem cl_d = clCreateBuffer(context, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR, sizeof(float)*N, d, &error);
	check_error("setting kernel args",error);
	assert(error == CL_SUCCESS);
	cl_mem cl_a = clCreateBuffer(context, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR, sizeof(float)*N, a, &error);
	check_error("setting kernel args",error);
	assert(error == CL_SUCCESS);
	cl_mem cl_I = clCreateBuffer(context, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR, sizeof(float)*N, I, &error);
	check_error("setting kernel args",error);
	assert(error == CL_SUCCESS);
	cl_mem cl_spikes = clCreateBuffer(context, CL_MEM_WRITE_ONLY, sizeof(unsigned int)*spike_array_length, NULL, &error);
	check_error("setting kernel args",error);
	assert(error == CL_SUCCESS);
	cl_mem cl_k = clCreateBuffer(context, CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR, sizeof(unsigned int), &(k[0]), &error);
	check_error("setting kernel args",error);
	assert(error == CL_SUCCESS);
	cl_mem cl_num_post = clCreateBuffer(context, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR, sizeof(unsigned int)*N, num_post, &error);
	check_error("setting kernel args",error);
	assert(error == CL_SUCCESS);
	/*
	cl_mem cl_pre_fired = clCreateBuffer(context, CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR, sizeof(bool)*T*M*N, pre_fired, &error);
	check_error("setting kernel args",error);
	assert(error == CL_SUCCESS);
	cl_mem cl_exc_input = clCreateBuffer(context, CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR, sizeof(float)*N, exc_input, &error);
	check_error("setting kernel args",error);
	assert(error == CL_SUCCESS);
	*/

	// Enqueuing parameters
	// Note that we inform the size of the cl_mem object, not the size of the memory pointed by it
	error = clSetKernelArg(neuron_fired, 0, sizeof(cl_mem), (void*)&cl_membranes);
	error |= clSetKernelArg(neuron_fired, 1, sizeof(cl_mem), (void*)&cl_u);
	error |= clSetKernelArg(neuron_fired, 2, sizeof(cl_mem), (void*)&cl_d);
	error |= clSetKernelArg(neuron_fired, 3, sizeof(cl_mem), (void*)&cl_a);
	error |= clSetKernelArg(neuron_fired, 4, sizeof(cl_mem), (void*)&cl_spikes);
	error |= clSetKernelArg(neuron_fired, 5, sizeof(cl_mem), (void*)&cl_k);
	error |= clSetKernelArg(neuron_fired, 6, sizeof(cl_float), (void*)&v_thresh);
	error |= clSetKernelArg(neuron_fired, 7, sizeof(cl_float), (void*)&v_reset);
	error |= clSetKernelArg(neuron_fired, 8, sizeof(cl_uint), (void*)&N);
	check_error("setting kernel args for neuron_fired",error);
	assert(error == CL_SUCCESS);

	// Enqueuing parameters
	// Note that we inform the size of the cl_mem object, not the size of the memory pointed by it
	error = clSetKernelArg(evolve_neuron, 0, sizeof(cl_mem), (void*)&cl_membranes);
	error |= clSetKernelArg(evolve_neuron, 1, sizeof(cl_mem), (void*)&cl_u);
	error |= clSetKernelArg(evolve_neuron, 2, sizeof(cl_mem), (void*)&cl_a);
	error |= clSetKernelArg(evolve_neuron, 3, sizeof(cl_mem), (void*)&cl_I);
	error |= clSetKernelArg(evolve_neuron, 4, sizeof(cl_uint), (void*)&N);
	error |= clSetKernelArg(evolve_neuron, 5, sizeof(cl_float), (void*)&h);
	check_error("setting kernel args for evolve_neuron",error);
	assert(error == CL_SUCCESS);

	//unsigned int check_k[1];

#ifdef WATCH_NEURONS
	float watched_membrane[num_watched_neurons][T];
#ifdef WATCH_ADAPTATION
	float watched_us[num_watched_neurons][T];
#endif
#endif

	INIT_TIMER(kernels)
	START_TIMER(kernels)
	for (unsigned int t=0; t<T; t++)
	{
		/* Step 1
		// random thalamic input for 1 in 1000 neurons
		*/
		error = clEnqueueNDRangeKernel(queue, neuron_fired, 1, NULL, &global_ws, &local_ws, 0, NULL, NULL);
		check_error("enqueueing NDRangeKernel",error);
		assert(error == CL_SUCCESS);
		error = clFinish(queue);
		check_error("waiting for NDRangeKernel",error);
		assert(error == CL_SUCCESS);
		//cout << "finished kernel neuron_fired" << endl;

		/* Step 1.5
		// read number recently detected spikes
		*/
		error = clEnqueueReadBuffer(queue, cl_k, CL_TRUE, 0, sizeof(unsigned int), &(k[t]), 0, NULL, NULL);
		check_error("enqueueing read buffer",error);
		assert(error == CL_SUCCESS);
		//cout << "Num spikes after " << t << " ms: " << k[t] << endl;

		/* Step 2
		// transmit loop
		*/

		/* Step 3
		// evolve neurons after spike transmission
		*/
		error = clEnqueueNDRangeKernel(queue, evolve_neuron, 1, NULL, &global_ws, &local_ws, 0, NULL, NULL);
		check_error("enqueueing NDRangeKernel",error);
		assert(error == CL_SUCCESS);
		error = clFinish(queue);
		check_error("waiting for NDRangeKernel",error);
		assert(error == CL_SUCCESS);
		//cout << "finished kernel evolve_neuron" << endl;

		/* Step 3.5
		// read watched membranes
		*/
#ifdef WATCH_NEURONS
		error = clEnqueueReadBuffer(queue, cl_membranes, CL_TRUE, 0, sizeof(float)*N, membranes, 0, NULL, NULL);
		check_error("enqueueing read buffer",error);
		assert(error == CL_SUCCESS);
#ifdef WATCH_ADAPTATION
		error = clEnqueueReadBuffer(queue, cl_u, CL_TRUE, 0, sizeof(float)*N, u, 0, NULL, NULL);
		check_error("enqueueing read buffer",error);
		assert(error == CL_SUCCESS);
#endif
		for (unsigned int i=0; i < num_watched_neurons; i++)
		{
			watched_membrane[i][t] = membranes[neurons_tobe_watched[i]];
#ifdef WATCH_ADAPTATION
			watched_us[i][t] = u[neurons_tobe_watched[i]];
#endif
		}
#endif

	}
	STOP_TIMER("kernels",kernels)

	/*
	* Read stuff back from openCL memories
	*/
	error = clEnqueueReadBuffer(queue, cl_spikes, CL_TRUE, 0, sizeof(unsigned int)*spike_array_length, spikes, 0, NULL, NULL);
	//error = clEnqueueReadBuffer(queue, cl_membranes, CL_TRUE, 0, sizeof(float)*N, membranes, 0, NULL, NULL);
	//error = clEnqueueReadBuffer(queue, cl_u, CL_TRUE, 0, sizeof(float)*N, u, 0, NULL, NULL);
	//error = clEnqueueReadBuffer(queue, cl_k, CL_TRUE, 0, sizeof(unsigned int), check_k, 0, NULL, NULL);
	check_error("enqueueing read buffer",error);
	assert(error == CL_SUCCESS);

	/*
	* Output results
	*/
	cout << "Results:" << endl;
	cout << "found: " << k[T-1] << endl;

	ofstream myfile;
#ifdef WATCH_NEURONS
	myfile.open("results/membrane.txt");
	for (int t=0; t<T; t++) {
		myfile << t;
		for (unsigned int n=0; n < num_watched_neurons; n++)
		{
			myfile << "," << watched_membrane[n][t];
#ifdef WATCH_ADAPTATION
			myfile << "," << watched_us[n][t];
#endif
		}
		myfile << endl;
	}
	myfile.close();
#endif

	myfile.open("results/spikes.txt");
	for (int t=1; t<T; t++)
	{
		for (unsigned int n=k[t-1]; n < k[t]; n++)
		{
			myfile << t;
			myfile << "," << spikes[n];
			myfile << endl;
		}
	}
	myfile.close();


	clReleaseKernel(neuron_fired);
	clReleaseKernel(evolve_neuron);
	clReleaseCommandQueue(queue);
	clReleaseContext(context);
	clReleaseMemObject(cl_spikes);

	STOP_TIMER("complete",complete)
}
