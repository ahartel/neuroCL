#include <stdlib.h>
#include <iostream>
#include <fstream>

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
	for (unsigned int i=0; i<T; i++) k[i] = 0;

	srand(42);
	init_neurons(membranes,u,d,a,I);

	unsigned int total_spikes = 0;

#ifdef WATCH_NEURONS
	float watched_membrane[num_watched_neurons][T];
#ifdef WATCH_ADAPTATION
	float watched_us[num_watched_neurons][T];
#endif
#endif

	INIT_TIMER(loops)
	for (unsigned int t=0; t<T; t++)
	{
		// reset loop
		for (int i=0; i<N; i++)
		{
			if (membranes[i] > v_thresh)
			{
				membranes[i] = v_reset;
				u[i] += d[i];
				spikes[k[t]++] = i;
			}
		}
		// evolve loop
		for (int i=0; i<N; i++)
		{
			membranes[i] += h * 0.5 * (( 0.04 * membranes[i] * membranes[i]) + (5.0 * membranes[i] + (140.0 - u[i] + I[i] )));
			membranes[i] += h * 0.5 * (( 0.04 * membranes[i] * membranes[i]) + (5.0 * membranes[i] + (140.0 - u[i] + I[i] )));

			//float y0_,ya,ya_,yb,yb_,yc,yc_;
			//y0_ = ( 0.04 * membranes[i] + 5.0 ) * membranes[i] + 140.0 - u[i] + I[i];
			//ya = membranes[i] + h * 0.5 * y0_;
			//ya_ = ( 0.04 * ya + 5.0 ) * ya + 140.0 - u[i] + I[i];
			//yb = membranes[i] + h * 0.5 * ya_;
			//yb_ = ( 0.04 * yb + 5.0 ) * yb + 140.0 - u[i] + I[i];
			//yc = membranes[i] + h * yb_;
			//yc_ = ( 0.04 * yc + 5.0 ) * yc + 140.0 - u[i] + I[i];
			//membranes[i] = membranes[i] + h/6.0 * (y0_ + 2 * (ya_+yb_) + yc_);

			u[i] += h * a[i]*(0.2*membranes[i]-u[i]);
		}

#ifdef WATCH_NEURONS
		for (unsigned int i=0; i < num_watched_neurons; i++)
		{
			watched_membrane[i][t] = membranes[neurons_tobe_watched[i]];
#ifdef WATCH_ADAPTATION
			watched_us[i][t] = u[neurons_tobe_watched[i]];
#endif
		}
#endif

		total_spikes += k[t];
		//cout << "Num spikes after " << t << " ms: " << total_spikes << endl;
	}

	STOP_TIMER("loops",loops)

	cout << "found: " << total_spikes << endl;

#ifdef WATCH_NEURONS
	ofstream myfile;
	myfile.open("membrane_compare.txt");
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

	STOP_TIMER("complete", complete)
}
