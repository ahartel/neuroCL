#include <stdlib.h>
#include <iostream>

#define TIMING
#include "helpers.h"

using namespace std;

int main()
{
	INIT_TIMER(complete)

	float membranes[N];
	float u[N];
	float d[N];
	float a[N];
	float I[N];
	unsigned int spikes[N];
	unsigned int k[T];
	k[0] = 0;

	srand(42);
	init_neurons(membranes,u,d,a,I);

	INIT_TIMER(loops)
	for (unsigned int t=0; t<T; t++)
	{
		for (int i=0; i<N; i++)
		{
			if (membranes[i] > v_thresh)
			{
				membranes[i] = v_reset;
				u[i] = d[i];
				spikes[k[t]++] = i;
			}
	/*
		}

		for (int idx=0; idx<N; idx++)
		{
	*/
			membranes[i] += 0.5 * ( ( 0.04 * membranes[i] + 5.0 ) * membranes[i] + 140.0 - u[i] + I[i] );
			membranes[i] += 0.5 * ( ( 0.04 * membranes[i] + 5.0 ) * membranes[i] + 140.0 - u[i] + I[i] );
			u[i]+=a[i]*(0.2*membranes[i]-u[i]);
		}
	}

	STOP_TIMER("loops",loops)

	unsigned int total_spikes = 0;
	for (unsigned int t=0; t<T; t++)
		total_spikes += k[t];
	cout << "found: " << total_spikes << endl;

	STOP_TIMER("complete", complete)
}
