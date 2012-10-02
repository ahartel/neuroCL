__kernel void evolve_neuron (
			__global float* membranes,
			__global float* u,
			__global const float* a,
			__global const float* I,
			const float v_thresh,
			const float v_reset,
			const int num)
{
	/* get_global_id(0) returns the ID of the thread in execution.
	As many threads are launched at the same time, executing the same kernel,
	each one will receive a different ID, and consequently perform a different computation.*/
	const int idx = get_global_id(0);

	/* Now each work-item asks itself: "is my ID inside the vector's range?"
	If the answer is YES, the work-item performs the corresponding computation*/
	if (idx < num)
	{
		membranes[idx] += 0.5 * ( ( 0.04 * membranes[idx] + 5.0 ) * membranes[idx] + 140.0 - u[idx] + I[idx] );
		membranes[idx] += 0.5 * ( ( 0.04 * membranes[idx] + 5.0 ) * membranes[idx] + 140.0 - u[idx] + I[idx] );
		u[idx]+=a[idx]*(0.2*membranes[idx]-u[idx]);
	}
}
