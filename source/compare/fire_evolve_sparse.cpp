#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <vector>
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
		myfile.open("results/spikes_compare_sparse.txt");
	else
		myfile.open("results/spikes_compare_sparse.txt",ios::app);
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

int get_delays(int delay_start[D],int synapseID)
{
	int last_delay_start = 0;
	int last_delay = 0;
	for (int d=0; d<D; d++)
	{
		if (delay_start[d] > synapseID && last_delay_start <= synapseID)
		{
			cout << "Found delay " << d << " for synapseID " << synapseID << endl;
			return d;
		}
		else if (delay_start[d] <= synapseID && d==D-1)
		{
			cout << "Found delay " << d << " for synapseID " << synapseID << endl;
			return d;
		}

		if (delay_start[d] > 0)
		{
			last_delay_start = delay_start[d];
			last_delay = d;
		}
	}

	throw "Delay not found";
}


int main()
{
	// neuron variables
	float membranes[N];
	float u[N];
	float d[N];
	float a[N];
	// input current
	float I[N][D];
	int delay_pointer = 0;
	// network parameters
	unsigned int delay_start[N][D];
	unsigned int delay_start_pre[N][D];
	unsigned int delay_count[N][D];
	unsigned int delay_count_pre[N][D];
	vector<float> weights[N];
	vector<unsigned int> post_neurons[N];
	vector<unsigned int> pre_neurons[N];
	unsigned int num_post[N] = {0};
	unsigned int num_pre[N] = {0};
	// STDP functions
	float LTP[N][1001+D], LTD[N];
	vector<float> sd[N];		// matrix of synaptic weights and their derivatives
	vector<float*> sd_pre[N];		// presynaptic weights
	// spike storage
	unsigned int spikes[int(N*1000*h)];
	unsigned int k[1000];
	for (unsigned int i=0; i<1000; i++) k[i] = 0;

	init_neurons_sparse(membranes,u,d,a,I,weights,delay_start,delay_count,post_neurons,num_post,LTP,LTD,sd,sd_pre,num_pre,pre_neurons,delay_start_pre,delay_count_pre);

	if (0)
	{
		for (unsigned int n=0; n<N; n++)
		{
			cout << "num_post[" << n << "]: " << num_post[n] << endl;
			if (num_post[n] > 0)
			{
				cout << "post-neurons: ";
				print_loop(post_neurons[n],post_neurons[n].size());
				cout << "weights: ";
				print_loop(weights[n],post_neurons[n].size());
				cout << "delay_start[" << n << "]: ";
				print_loop(delay_start[n],D);
				cout << "delay_count[" << n << "]: ";
				print_loop(delay_count[n],D);
			}
		}
		//return 0;
	}

	INIT_TIMER(complete)
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
			/* Step 1
			// random thalamic input for 1 in 1000 neurons
			*/
			for (int j=0;j<N/1000+1;j++)
				I[getrandom(N)][delay_pointer] += 20.0;
			// reset loop
			for (int i=0; i<N; i++)
			{
				if (membranes[i] > v_thresh)
				{
					membranes[i] = v_reset;
					u[i] += d[i];
					spikes[total_spikes++] = i;
					LTP[i][t+D]= 0.1;
					LTD[i]=0.12;

					// get first used delay value
					int d = 0;
					int delay_cnt = 0;
					// loop through delay values until a used one is found
					while (delay_count_pre[i][d] == 0 && d < D)
					{
						d++;
					}
					// save the number of neurons that use this delay
					delay_cnt = delay_count_pre[i][d];
					// check if the loop has run through without result
					if (d == D-1 && delay_count_pre[i][d] == 0)
					{
						cout << "Neuron " << i;
						throw "Error: No delays found.";
					}
					for (unsigned int j=0; j<num_pre[i]; ++j)
					{
						*sd_pre[i][j] += LTP[pre_neurons[i][j]][t+D-d-1];// this spike was after pre-synaptic spikes
						// keep track of the number of remaining neurons with current delay
						// if necessary, increase delay value
						if (--delay_cnt == 0)
						{
							d++;
							while (delay_count_pre[i][d] == 0 && d < D)
							{
								d++;
							}
							delay_cnt = delay_count_pre[i][d];
						}
					}
				}
			}

			/* Step 2
			// transmit loop
			*/
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
				//cout << "Found spike in neuron " << neuronID << endl;
				// skip this neuron, if it doesn't connect
				// (this can happen in a randomly connected network)
				if (num_post[neuronID] == 0)
					continue;
				// get first used delay value
				int d = 0;
				int delay_cnt = 0;
				// loop through delay values until a used one is found
				while (delay_count[neuronID][d] == 0 && d < D)
				{
					d++;
				}
				// save the number of neurons that use this delay
				delay_cnt = delay_count[neuronID][d];
				// check if the loop has run through without result
				if (d == D-1 && delay_count[neuronID][d] == 0)
				{
					cout << "Neuron " << neuronID;
					throw "Error: No delays found.";
				}
				// Now, generate the current for all post-syn. neurons
				for (int m=0; m<num_post[neuronID]; m++)
				{
					unsigned int pst_n = post_neurons[neuronID][m];
					if (pst_n < Ne) // this spike is before postsynaptic spikes
						sd[neuronID][m] -= LTD[m];

					I[pst_n][(d+delay_pointer+1)%D] += weights[neuronID][m];
					// keep track of the number of remaining neurons with current delay
					// if necessary, increase delay value
					if (--delay_cnt == 0)
					{
						d++;
						while (delay_count[neuronID][d] == 0 && d < D)
						{
							d++;
						}
						delay_cnt = delay_count[neuronID][d];
					}
				}
			}

			/* Step 3
			// evolve neurons after spike transmission
			*/
			for (int i=0; i<N; i++)
			{

				membranes[i] += h * 0.5 * (( 0.04 * membranes[i] * membranes[i]) + (5.0 * membranes[i] + (140.0 - u[i] + I[i][delay_pointer] )));
				membranes[i] += h * 0.5 * (( 0.04 * membranes[i] * membranes[i]) + (5.0 * membranes[i] + (140.0 - u[i] + I[i][delay_pointer] )));

				u[i] += h * a[i]*(0.2*membranes[i]-u[i]);

				//reset current
				I[i][delay_pointer] = 0;

				LTP[i][t+D+1]=0.95*LTP[i][t+D];
				LTD[i]*=0.95;
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
		} // end of millisecond loop


		for (unsigned int i=0;i<N;i++)		// prepare for the next sec
		{
			for (unsigned int j=0;j<D+1;j++)
				LTP[i][j] = LTP[i][1000+j];

			if (i< Ne)
				for (unsigned int j=0;j<num_post[i];j++)
				{
					weights[i][j] += 0.01+sd[i][j];
					sd[i][j] *= 0.9;
					if (weights[i][j]>sm) weights[i][j]=sm;
					if (weights[i][j]<0) weights[i][j]=0.0;
				}
		}

		STOP_TIMER("one second loop",loops)

		cout << "Second " << sec << ": found " << total_spikes << endl;

		write_spikes(sec,spikes,k);

#ifdef WATCH_NEURONS
		write_watched_membranes(sec,watched_membrane,watched_us);
#endif

		total_spikes = 0;
	} // end of 'second' loop
	STOP_TIMER("complete", complete)
}



