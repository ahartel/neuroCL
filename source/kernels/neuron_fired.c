
__kernel void neuron_fired (
			__global const float* membranes,
			__global unsigned int* spikes,
			__global unsigned int* k,
			const float threshold,
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
		if (membranes[idx] > threshold)
		{
			unsigned int p = atomic_add(k,1);
			spikes[p] = idx;
		}
	}
}
