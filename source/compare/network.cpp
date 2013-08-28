#include <stdlib.h>
#include <stdexcept>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include "network.h"

using namespace std;

std::vector<unsigned int> Network::get_last_spikes()
{
	//cout << "Sec : " << sec << endl;
	//cout << "last Sec : " << last_sec << endl;
	//cout << "milli : " << t << endl;
	//if (t>0)
	//	cout << "spikes: " << k[t-1] << endl;
	//if (t>1)
	//	cout << "last spikes: " << k[t-2] << endl;

	//for (size_t ii=0;ii<k[t-1];ii++)
	//	cout << spikes[ii] << ",";
	//cout << endl;

	if (t==0)
	{
		if (sec == 0)
		{
			return std::vector<unsigned int>();
		}
		else
		{
			if (k[999] > k[998])
			{
				std::vector<unsigned int> return_vec(k[999]-k[998]);
				std::vector<unsigned int>::iterator it_start = spikes.begin();
				std::vector<unsigned int>::iterator it_end = spikes.begin();
				it_start = it_start+(k[998]);
				it_end = it_end+(k[999]);

				copy(it_start,it_end,return_vec.begin());
				return return_vec;
			}
			else
				return std::vector<unsigned int>();
		}
	}
	else if (t==1)
	{
		std::vector<unsigned int> return_vec(k[0]);
		copy(spikes.begin(),spikes.begin()+k[0],return_vec.begin());
		return return_vec;
	}
	else
	{
		if (k[t-1] > k[t-2])
		{
			std::vector<unsigned int>::iterator it_start = spikes.begin();
			std::vector<unsigned int>::iterator it_end = spikes.begin();
			it_start = it_start+(k[t-2]);
			it_end = it_end+(k[t-1]);
			std::vector<unsigned int> return_vec(it_start,it_end);

			return return_vec;
		}
		else
			return std::vector<unsigned int>();
	}
}

int Network::get_delays(int delay_start[D],int synapseID)
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

	throw std::logic_error("Delay not found");
}

void Network::add_spikes(vector<pre_spike> const& pre_spikes)
{
	for (pre_spike ps : pre_spikes)
	{
		if (ps.delay > D)
		{
			cout << "Error: Spike with too large delay in add_spikes. Skipping this spike" << endl;
			continue;
		}
		I[ps.neuron][(ps.delay+delay_pointer)%(D+1)] += ps.weight;
	}
}

void Network::eject_dopamine()
{
	DA = 1.0;
	cout << "Dopamine" << endl;
}

Network::Network() :
	name("default_name"),
	delay_pointer(0),
	membranes(N,0),
	u(N,0),
	d(N,0),
	a(N,0),
	sd{N},
	sd_pre{N},
	spikes(1000*N),
	DA(0)
{
	init();
}

Network::Network(std::string const& n) :
	name(n),
	delay_pointer(0),
	membranes(N,0),
	u(N,0),
	d(N,0),
	a(N,0),
	sd{N},
	sd_pre{N},
	spikes(1000*N),
	DA(0)
{
	init();
}

void Network::init()
{
	delay_pointer = 0;
	// spike storage
	for (unsigned int i=0; i<1000; i++) k[i] = 0;

	init_neurons_sparse(membranes,u,d,a,I,weights,delay_start,delay_count,post_neurons,num_post,LTP,LTD,sd,sd_pre,num_pre,pre_neurons,delay_start_pre,delay_count_pre);

#ifdef DEBUG_OUTPUT
	{
		for (unsigned int n=0; n<N; n++)
		{
			cout << "Neuron " << n << endl;
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
			cout << "num_pre[" << n << "]: " << num_pre[n] << endl;
			if (num_pre[n] > 0)
			{
				cout << "pre-neurons: ";
				print_loop(pre_neurons[n],pre_neurons[n].size());
				cout << "delay_start_pre[" << n << "]: ";
				print_loop(delay_start_pre[n],D);
				cout << "delay_count_pre[" << n << "]: ";
				print_loop(delay_count_pre[n],D);
			}
		}
		//return 0;
	}
#endif

	INIT_TIMER(complete)
	last_total_spikes = 0;
	total_spikes = 0;
	sec = 0;
	last_sec = 0;
	t = 0;

#ifdef WATCH_NEURONS
	std::vector<std::vector<float> > watched_LTP {num_watched_neurons};
	std::vector<std::vector<float> > watched_LTD {num_watched_neurons};
	std::vector<std::vector<float> > watched_membrane {num_watched_neurons};
	std::vector<std::vector<float> > watched_us {num_watched_neurons};
#endif
}

void Network::step()
{
	DA *= 0.95;
	for (unsigned int i=0;i<N;i++)		// prepare for the next sec
		if (i< Ne && num_post[i] > 0)
			for (unsigned int j=0;j<num_post[i];j++)
				sd[i][j] *= 0.9;

	//for (unsigned int sec=0; sec<T; sec++)
	//{
		//INIT_TIMER(loops)
		//for (unsigned int t=0;t<1000;t++)
		//{
			/* Step 1
			// random thalamic input for 1 in 1000 neurons
			*/
			for (int j=0;j<N/1000+1;j++)
				I[getrandom(N)][delay_pointer] += 20.0;
			// reset loop
			for (int i=0; i<N; i++)
			{
				// has the neuron fired?
				if (membranes[i] > v_thresh)
				{
					membranes[i] = v_reset;
					u[i] += d[i];
					spikes[k[t]++] = i;
					total_spikes = k[t];
					LTP[i][t+D]= 0.1;
					LTD[i]=0.12;

					// get first used pre-synaptic delay value
					int d = 0;
					int delay_cnt = 0;
					// loop through delay values until a used one is found
					while (delay_count_pre[i][d] == 0 && d < D-1)
					{
						++d;
					}
					// save the number of neurons that use this delay
					delay_cnt = delay_count_pre[i][d];
					// check if the loop has run through without result
					if (d == D-1 && delay_cnt == 0)
					{
						// No problem, neurons may fire without having presynaptic connections
						cout << "Warning: Neuron " << i << " fired without having any presynaptic connections." << endl;
					}
					else
					{
						for (unsigned int j=0; j<num_pre[i]; ++j)
						{
							// if delay_cnt is zero there is an error in the delay table
							if (delay_cnt == 0 && j<num_pre[i]-1)
							{
								stringstream bla;
								bla << "Error, in Neuron " << i << ": No more delays left, but there are still pre-synaptic neurons left.";
								throw std::logic_error(bla.str());
							}

							//cout << "Neuron " << i << " pre " << pre_neurons[i][j] << endl;
							*sd_pre[i][j] += LTP[pre_neurons[i][j]][t+D-d-1];// this spike was after pre-synaptic spikes
							// keep track of the number of remaining neurons with current delay
							// if necessary, increase delay value
							if (--delay_cnt == 0)
							{
								++d;
								while (delay_count_pre[i][d] == 0 && d < D-1)
								{
									++d;
								}
								if (d<D)
									delay_cnt = delay_count_pre[i][d];
							}
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
				while (delay_count[neuronID][d] == 0 && d < D-1)
				{
					++d;
				}
				// save the number of neurons that use this delay
				delay_cnt = delay_count[neuronID][d];
				// check if the loop has run through without result
				if (d == D-1 && delay_cnt == 0)
				{
					// Neurons may fire without having post-synaptic neurons
				}
				else
				{
					// Now, generate the current for all post-syn. neurons
					for (int m=0; m<num_post[neuronID]; m++)
					{
						if (delay_cnt == 0 && m<num_post[neuronID]-1)
						{
							stringstream bla;
							bla << "Error, in Neuron " << neuronID << ": No more delays left, but there are still post-synaptic neurons left.";
							throw std::logic_error(bla.str());
						}
						unsigned int pst_n = post_neurons[neuronID][m];
						if (pst_n < Ne) // this spike is before postsynaptic spikes
							sd[neuronID][m] -= LTD[m];

						I[pst_n][(d+delay_pointer+1)%(D+1)] += weights[neuronID][m];
						// keep track of the number of remaining neurons with current delay
						// if necessary, increase delay value
						if (--delay_cnt == 0)
						{
							++d;
							while (delay_count[neuronID][d] == 0 && d < D-1)
							{
								++d;
							}
							if (d<D)
								delay_cnt = delay_count[neuronID][d];
						}
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
				watched_membrane[i].push_back(membranes[neurons_tobe_watched[i]]);
	#ifdef WATCH_ADAPTATION
				watched_us[i].push_back(u[neurons_tobe_watched[i]]);
	#endif
	#ifdef WATCH_DERIVATIVES
				//watched_sd[i][j].push_back(sd[i][j]);
				watched_LTP[i].push_back(LTP[i][t+D]);
				watched_LTD[i].push_back(LTD[i]);
	#endif
			}
	#endif

			// jump to next delay entry
			delay_pointer++;
			if (delay_pointer == D+1)
				delay_pointer = 0;
		//} // end of millisecond loop

		++t;

		if (t%10==0 && DA > 1e-4)
		{
			cout << "DA level :" << DA << endl;
			for (unsigned int i=0;i<N;i++)		// prepare for the next sec
			{
				if (i< Ne && num_post[i] > 0)
				{
#ifdef DEBUG_OUTPUT
					cout << "Neuron " << i << "'s weights: ";
#endif
					for (unsigned int j=0;j<num_post[i];j++)
					{
						weights[i][j] += DA*(0.01+sd[i][j]);
						if (weights[i][j]>sm) weights[i][j]=sm;
						if (weights[i][j]<0) weights[i][j]=0.0;
#ifdef DEBUG_OUTPUT
						cout << " " << weights[i][j];
#endif
					}
#ifdef DEBUG_OUTPUT
					cout << endl;
					cout << "Neuron " << i << "'s deerivatives: ";
					for (unsigned int j=0;j<num_post[i];j++)
						cout << " " << sd[i][j];
					cout << endl;
#endif
				}
			}
		}

		if (t==1000)
		{
			for (unsigned int i=0;i<N;i++)		// prepare for the next sec
			{
				for (unsigned int j=0;j<D+1;j++)
					LTP[i][j] = LTP[i][1000+j];
			}


			//STOP_TIMER("one second loop",loops)

			cout << "Second " << sec << ": found " << total_spikes << endl;
			cout << "Second " << sec << ": firing rate " << total_spikes/N << endl;

			write_spikes("spikes_"+name+".txt",sec,spikes,k);

#ifdef WATCH_NEURONS
			write_watched_membranes("membrane_"+name+".txt",sec,watched_membrane,watched_us);
		#ifdef WATCH_DERIVATIVES
			write_derivatives("stdp_"+name+".txt",sec,watched_LTP,watched_LTD);
		#endif
#endif

			total_spikes = 0;
			k[0] = 0;
			t = 0;
			last_sec = sec++;
		}
		else
		{
			k[t] = k[t-1];
		}
		// end of 'second' loop
	//STOP_TIMER("complete", complete)

}



