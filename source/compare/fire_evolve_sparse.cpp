#include <stdlib.h>
#include <iostream>
#include <fstream>

#include "helpers.h"

using namespace std;

void write_watched_membranes(unsigned int sec, float watched_membrane[num_watched_neurons][1000], float watched_us[num_watched_neurons][1000])
{
	ofstream myfile;
	if (sec==0)
		myfile.open("results/membrane_compare.txt");
	else
		myfile.open("results/membrane_compare.txt",ios::app);
	for (int t=0; t<1000; t++) {
		myfile << sec*1000+t;
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
}

void write_spikes(unsigned int sec, unsigned int* spikes, unsigned int* k)
{
	ofstream myfile;
	if (sec==0)
		myfile.open("results/spikes_compare.txt");
	else
		myfile.open("results/spikes_compare.txt",ios::app);
	for (int t=1; t<1000; t++)
	{
		for (unsigned int n=k[t-1]; n < k[t]; n++)
		{
			myfile << sec*1000+t;
			myfile << "," << spikes[n];
			myfile << endl;
		}
	}
	myfile.close();

}

int main()
{
	INIT_TIMER(complete)
	// neuron variables
	float membranes[N];
	float u[N];
	float d[N];
	float a[N];
	// input current
	float I[N][D];
	int delay_pointer = 0;
	// network parameters
	float weights[N][M];
	int delays[N][M];
	int post_neurons[N][M];
	int num_post[N];
	// spike storage
	unsigned int spikes[int(N*1000*h)];
	unsigned int k[1000];
	for (unsigned int i=0; i<1000; i++) k[i] = 0;

	init_neurons(membranes,u,d,a,I,weights,delays,post_neurons,num_post);

	unsigned int total_spikes = 0;

#ifdef WATCH_NEURONS
	float watched_membrane[num_watched_neurons][1000];
	float watched_us[num_watched_neurons][1000];
#endif

	for (unsigned int sec=0; sec<T; sec++)
	{
		INIT_TIMER(loops)
		for (unsigned int t=0;t<1000;t++)
		{
			// random thalamic input for 1 in 1000 neurons
			for (int j=0;j<N/1000;j++)
				I[getrandom(N)][delay_pointer] += 20.0;
			// reset loop
			for (int i=0; i<N; i++)
			{
				if (membranes[i] > v_thresh)
				{
					membranes[i] = v_reset;
					u[i] += d[i];
					spikes[total_spikes++] = i;
				}
			}

			// transmit loop
			int last_spike_count;
			if (t > 0)
			{
				last_spike_count = k[t-1];
			}
			else
				last_spike_count = 0;

			for (int i=last_spike_count; i<total_spikes; i++)
			{
				int neuronID = spikes[i];
				for (int m=0; m<num_post[neuronID]; m++)
				{	
					I[post_neurons[neuronID][m]][(delays[neuronID][m]+delay_pointer+1)%D] += weights[neuronID][m];	
					//if (post_neurons[neuronID][m]==0)
					//{
					//	cout << "Time " << t << ": Injecting current to neuron " << post_neurons[neuronID][m] << " from neuron " << neuronID << " into delay slot " << (delays[neuronID][m]+delay_pointer+1)%D << endl;
					//	cout << I[post_neurons[neuronID][m]][(delays[neuronID][m]+delay_pointer+1)%D] << endl;
					//}
				}
			}

			// evolve loop
			for (int i=0; i<N; i++)
			{

				membranes[i] += h * 0.5 * (( 0.04 * membranes[i] * membranes[i]) + (5.0 * membranes[i] + (140.0 - u[i] + I[i][delay_pointer] )));
				membranes[i] += h * 0.5 * (( 0.04 * membranes[i] * membranes[i]) + (5.0 * membranes[i] + (140.0 - u[i] + I[i][delay_pointer] )));

				u[i] += h * a[i]*(0.2*membranes[i]-u[i]);

				//reset current
				I[i][delay_pointer] = 0;
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
			k[t] = total_spikes;
			//cout << "Num spikes after " << t << " ms: " << k[t] << endl;

			// jump to next delay entry
			delay_pointer++;
			if (delay_pointer == D)
				delay_pointer = 0;
		}

		STOP_TIMER("one second loop",loops)

		cout << "Second " << sec << ": found " << total_spikes << endl;

		write_spikes(sec,spikes,k);

#ifdef WATCH_NEURONS
		write_watched_membranes(sec,watched_membrane,watched_us);
#endif

		total_spikes = 0;
	}
	STOP_TIMER("complete", complete)
}



