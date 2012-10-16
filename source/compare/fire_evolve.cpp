#include <stdlib.h>
#include <iostream>
#include <fstream>

#define TIMING
#include "helpers.h"

using namespace std;

int main()
{
	INIT_TIMER(complete)

	double membranes[N];
	double u[N];
	float d[N];
	float a[N];
	float I[N];
	unsigned int spikes[N];
	unsigned int k[T];
	for (unsigned int i=0; i<T; i++) k[i] = 0;

	srand(42);
	init_neurons(membranes,u,d,a,I);

	unsigned int total_spikes = 0;

	double watched_membrane[1][T];

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
		}

		for (int i=0; i<N; i++)
		{
			membranes[i] += 0.5 * ( ( 0.04 * membranes[i] + 5.0 ) * membranes[i] + 140.0 - u[i] + I[i] );
			membranes[i] += 0.5 * ( ( 0.04 * membranes[i] + 5.0 ) * membranes[i] + 140.0 - u[i] + I[i] );
			//membranes[i] += 0.25 * ( ( 0.04 * membranes[i] + 5.0 ) * membranes[i] + 140.0 - u[i] + I[i] );
			//membranes[i] += 0.25 * ( ( 0.04 * membranes[i] + 5.0 ) * membranes[i] + 140.0 - u[i] + I[i] );
			u[i]+=a[i]*(0.2*membranes[i]-u[i]);

			watched_membrane[0][t] = membranes[0];
		}

		total_spikes += k[t];
		//cout << "Num spikes after " << t << " ms: " << total_spikes << endl;
	}

	STOP_TIMER("loops",loops)

	cout << "found: " << total_spikes << endl;

	ofstream myfile;
	myfile.open("membrane_compare.txt");
	for (int i=0; i<T; i++) {
		myfile << watched_membrane[0][i] << endl;
	}
	myfile.close();

	STOP_TIMER("complete", complete)
}
