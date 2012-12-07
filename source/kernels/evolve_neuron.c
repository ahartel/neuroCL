
__kernel void evolve_neuron (
			__global float* membranes,
			__global float* u,
			__global const float* a,
			__global const float* I,
			const int num,
			const float h
		)
{
	/* get_global_id(0) returns the ID of the thread in execution.
	As many threads are launched at the same time, executing the same kernel,
	each one will receive a different ID, and consequently perform a different computation.*/
	const int idx = get_global_id(0);

	/* Now each work-item asks itself: "is my ID inside the vector's range?"
	If the answer is YES, the work-item performs the corresponding computation*/
	if (idx < num)
	{
		membranes[idx] += h * 0.5 * (( 0.04 * membranes[idx] * membranes[idx]) + (5.0 * membranes[idx] + (140.0 - u[idx] + I[idx] )));
		membranes[idx] += h * 0.5 * (( 0.04 * membranes[idx] * membranes[idx]) + (5.0 * membranes[idx] + (140.0 - u[idx] + I[idx] )));
/*
		float y0_,ya,ya_,yb,yb_,yc,yc_;
		y0_ = ( 0.04 * membranes[idx] + 5.0 ) * membranes[idx] + 140.0 - u[idx] + I[idx];
		ya = membranes[idx] + h * 0.5 * y0_;
		ya_ = ( 0.04 * ya + 5.0 ) * ya + 140.0 - u[idx] + I[idx];
		yb = membranes[idx] + h * 0.5 * ya_;
		yb_ = ( 0.04 * yb + 5.0 ) * yb + 140.0 - u[idx] + I[idx];
		yc = membranes[idx] + h * yb_;
		yc_ = ( 0.04 * yc + 5.0 ) * yc + 140.0 - u[idx] + I[idx];
		membranes[idx] = membranes[idx] + h/6.0 * (y0_ + 2 * (ya_+yb_) + yc_);
*/
		u[idx] += h * a[idx]*(0.2*membranes[idx]-u[idx]);
	}
}


