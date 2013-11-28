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

void Network::add_current(vector<pre_spike> const& pre_spikes)
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

#ifdef WITH_DA
void Network::eject_dopamine()
{
	DA = 1.0;
	cout << "Dopamine" << endl;
}
#endif

void Network::enable_random_background()
{
	background_spikes = true;
}

Network::Network(unsigned int e, unsigned int i, unsigned int m) :
	Ne(e),
	Ni(i),
	N(e+i),
	M(m),
	name("default_name"),
	delay_pointer(0),
	membranes(N,0),
	u(N,0),
	d(N,0),
	a(N,0),
	sd{N},
	sd_pre{N},
#ifdef WITH_DA
	DA(0),
#endif
	spikes(1000*N),
#ifdef WATCH_DERIVATIVES
	watched_LTP {std::vector<unsigned int>(WATCHED_NEURONS).size()},
	watched_LTD {std::vector<unsigned int>(WATCHED_NEURONS).size()},
#endif
#ifdef WATCH_NEURONS
	watched_membrane {std::vector<unsigned int>(WATCHED_NEURONS).size()},
	watched_us {std::vector<unsigned int>(WATCHED_NEURONS).size()},
#endif
	background_spikes(false)
{
	init();
}

Network::Network(std::string const& n, NetworkDescription const& spec) :
	N(spec.getNumberOfNeurons()),
	Ne(spec.getNumberOfExcNeurons()),
	Ni(spec.getNumberOfInhNeurons()),
	name(n),
	delay_pointer(0),
	membranes(N,0),
	u(N,0),
	d(N,0),
	a(N,0),
	sd{N},
	sd_pre{N},
#ifdef WITH_DA
	DA(0),
#endif
	spikes(1000*N),
#ifdef WATCH_DERIVATIVES
	watched_LTP {std::vector<unsigned int>(WATCHED_NEURONS).size()},
	watched_LTD {std::vector<unsigned int>(WATCHED_NEURONS).size()},
#endif
#ifdef WATCH_NEURONS
	watched_membrane {std::vector<unsigned int>(WATCHED_NEURONS).size()},
	watched_us {std::vector<unsigned int>(WATCHED_NEURONS).size()},
#endif
	background_spikes(false)
{
	init_from_spec(spec);
}

Network::Network(std::string const& n,unsigned int e, unsigned int i, unsigned int m) :
	Ne(e),
	Ni(i),
	N(e+i),
	M(m),
	name(n),
	delay_pointer(0),
	membranes(N,0),
	u(N,0),
	d(N,0),
	a(N,0),
	sd{N},
	sd_pre{N},
#ifdef WITH_DA
	DA(0),
#endif
	spikes(1000*N),
#ifdef WATCH_DERIVATIVES
	watched_LTP {std::vector<unsigned int>(WATCHED_NEURONS).size()},
	watched_LTD {std::vector<unsigned int>(WATCHED_NEURONS).size()},
#endif
#ifdef WATCH_NEURONS
	watched_membrane {std::vector<unsigned int>(WATCHED_NEURONS).size()},
	watched_us {std::vector<unsigned int>(WATCHED_NEURONS).size()},
#endif
	background_spikes(false)
{
	init();
}

Network::~Network()
{
	for (size_t ii=0; ii<N;++ii)
	{
		delete [] delay_start[ii];
		delete [] delay_start_pre[ii];
		delete [] delay_count[ii];
		delete [] delay_count_pre[ii];
		delete [] LTP[ii];
		delete [] I[ii];
	}
	delete [] delay_start;
	delete [] delay_start_pre;
	delete [] delay_count;
	delete [] delay_count_pre;
	delete [] LTP;
	delete [] I;

	delete[] weights;
	delete[] post_neurons;
	delete[] pre_neurons;
	delete[] num_post;
	delete[] num_pre;
	delete[] LTD;

}

void Network::init()
{
	delay_pointer = 0;
	// spike storage
	for (unsigned int i=0; i<1000; i++) k[i] = 0;

	delay_start = new unsigned int*[N];
	delay_start_pre = new unsigned int*[N];
	delay_count = new unsigned int*[N];
	delay_count_pre = new unsigned int*[N];
	weights = new vector<float>[N];
	post_neurons = new vector<unsigned int>[N];
	pre_neurons = new vector<unsigned int>[N];
	num_post = new unsigned int[N];
	num_pre = new unsigned int[N];
	LTD = new float[N];
	LTP = new float*[N];
	I = new float*[N];

	init_random();

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



}


void Network::init_from_spec(NetworkDescription const& spec)
{
	delay_pointer = 0;
	// spike storage
	for (unsigned int i=0; i<1000; i++) k[i] = 0;

	delay_start = new unsigned int*[N];
	delay_start_pre = new unsigned int*[N];
	delay_count = new unsigned int*[N];
	delay_count_pre = new unsigned int*[N];
	weights = new vector<float>[N];
	post_neurons = new vector<unsigned int>[N];
	pre_neurons = new vector<unsigned int>[N];
	num_post = new unsigned int[N];
	num_pre = new unsigned int[N];
	LTD = new float[N];
	LTP = new float*[N];
	I = new float*[N];

	// This vector stores all pre-synaptic delays temporarily
	vector<unsigned int> pre_delays[N][D];

	for (int i=0; i<N; i++)
	{
		delay_start[i] = new unsigned int[D];
		delay_start_pre[i] = new unsigned int[D];
		delay_count[i] = new unsigned int[D];
		delay_count_pre[i] = new unsigned int[D];

		for (int d=0; d<D; ++d)
			delay_count_pre[i][d] = 0;
		num_pre[i] = 0;
	}

	/*
	 * init neuron parameters
	 */
	// iterate over neurons
	for (int i=0; i<N; i++) {
		I[i] = new float[D+1];
		LTP[i] = new float[1001+D];

		membranes[i] = v_reset;

		for (int del=0; del<=D; del++)
		{
			I[i][del] = 0;//(float)rand()/float(RAND_MAX)*10.0;
		}
		u[i] = 0.2*membranes[i];
		if (i < Ne) {
			a[i] = 0.02;
			d[i] = 8.0;
		}
		else {
			a[i] = 0.1;
			d[i] = 2.0;
		}
	}

	/*
	 * init connectivity
	 */
	// iterate over neurons
	for (int i=0; i<N; i++) {
		// pick how many post-synaptic connections this neuron will have
		num_post[i] = spec.getNumberConnectionsOfNeuron(i);
		// create temporary objects for delays, post-synaptic neurons and weights
		int delays[M];
		// the array synapse_numbers holds all synapse numbers for a given delay
		vector<int> synapse_numbers[D];
		vector<int> pst_n;
		vector<float> wgt_n;
		// reset STDP variables
		LTD[i] = 0.0;
		for (unsigned int j=0;j<1000+D;++j)
			LTP[i][j] = 0.0;
		// set delay_count and delay_start initially to zero
		for (int d=0; d<D; d++)
		{
			delay_start[i][d] = 0;
			delay_count[i][d] = 0;
		}
		// create values for connectivity
		// This will set weights and pick delays for each
		// of the neuron's post-synaptic connections
		for (int m=0; m<num_post[i]; m++)
		{
			// Initialize weights to fixed values
			std::tuple<float, unsigned int> conn = spec.getConnectionOfNeuron(i,m);
			wgt_n.push_back(get<0>(conn));

			// set weight derivatives initially to zero
			sd[i].push_back(0.0);
			// pick a delay and sort it into the delay_count array
			delays[m] = get<1>(conn);
			synapse_numbers[delays[m]].push_back(m);
			delay_count[i][delays[m]]++;
			// pick a post-synaptic neuron number
			unsigned int tmp_post = spec.getPostOfNeuron(i,m);
			// make sure the post-syn. neuron is not the pre-syn. neuron
			if (tmp_post == i)
				throw std::logic_error("Neuron with self-connection found");

			pst_n.push_back(tmp_post);

#ifdef DEBUG_OUTPUT
			cout << "Neuron " << i << ": Added post-neuron " << pst_n.back() << " with delay " << delays[m] << endl;
#endif

			// update pre-synaptic connection arrays
			// add the current neuron to the pre_neuron list
			// of the previously chosen post-synaptic neuron
			pre_delays[tmp_post][delays[m]].push_back(num_pre[tmp_post]);
			delay_count_pre[tmp_post][delays[m]]++;
			//pre_delays[tmp_post].push_back(delays[m]);
			pre_neurons[tmp_post].push_back(i);
			++num_pre[tmp_post];
		}
		// sort connectivity values according to delay values
		vector<unsigned int> pst_n_sort;
		vector<float> wgt_n_sort;
		for (int d=0; d<D; d++)
		{
			vector<int>::iterator delay_it = synapse_numbers[d].begin();
			while (delay_it != synapse_numbers[d].end())
			{
				pst_n_sort.push_back(pst_n[*delay_it]);
				wgt_n_sort.push_back(wgt_n[*delay_it]);
				/*
				   actually, sd would have to be sorted here, too!
				   But since it's initialized with zeros anyway,
				   we can just skip it here
			    */

				delay_it++;
			}
		}

		// add values to global arrays
		weights[i] = wgt_n_sort;
		post_neurons[i] = pst_n_sort;
		/*
		   The next loop calculates the correct values for
		   the delay_count and delay_start arrays, according
		   to Nageswaran et al.
		*/
		int last_delay = 0;
		for (int d=0; d<D; d++)
		{
			if (delay_count[i][d] > 0)
				delay_start[i][d] = last_delay;
			else
				delay_start[i][d] = 0;
			last_delay += delay_count[i][d];
		}
	}
	// iterate over neurons again
	// re-arrange pre-syn. information
	for (int i=0; i<N; i++)
	{
		// get pointers to the derivatives of the weights
		for (int m=0; m<num_post[i]; m++)
			sd_pre[post_neurons[i][m]].push_back(&sd[i][m]);// pointer to the derivative

		// sort connectivity values according to delay values
		vector<unsigned int> pre_n_sort;
#ifdef DEBUG_OUTPUT
		cout << "Neuron " << i << endl;
		cout << "post_neurons: " << post_neurons[i].size() << endl;
		cout << "pre_neurons: " << pre_neurons[i].size() << endl;
#endif
		if (!pre_neurons[i].empty())
		{
			for (int d=0; d<D; d++)
			{
				//cout << "pre_delays[" << d << "]";
				//print_loop(pre_delays[i][d]);
				for (unsigned int delay_it : pre_delays[i][d])
				{
					pre_n_sort.push_back(pre_neurons[i][delay_it]);
				}
			}
			pre_neurons[i] = pre_n_sort;
		}

		int last_delay = 0;
		for (int d=0; d<D; d++)
		{
			if (delay_count_pre[i][d] > 0)
				delay_start_pre[i][d] = last_delay;
			else
				delay_start_pre[i][d] = 0;
			last_delay += delay_count_pre[i][d];
		}
	}
}

void Network::step()
{
#ifdef WITH_DA
	DA *= 0.95;

	for (unsigned int i=0;i<N;i++)		// prepare for the next sec
		if (i< Ne && num_post[i] > 0)
			for (unsigned int j=0;j<num_post[i];j++)
				sd[i][j] *= 0.9;
#endif

	//for (unsigned int sec=0; sec<T; sec++)
	//{
		//INIT_TIMER(loops)
		//for (unsigned int t=0;t<1000;t++)
		//{
			/* Step 1
			// random thalamic input for 1 in 1000 neurons
			*/
			if (background_spikes)
				for (int j=0;j<N/1000+1;j++)
				{
					unsigned int target = getrandom(N);
#ifdef DEBUG_OUTPUT

					cout << "Injecting background spike into neuron " <<  target << endl;
#endif
					I[target][delay_pointer] += 20.0;
				}
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
						//cout << "Warning: Neuron " << i << " fired without having any presynaptic connections." << endl;
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

							cout << "Neuron " << i << " pre " << pre_neurons[i][j] << " adding " << LTP[pre_neurons[i][j]][t+D-d-1] << " to synapse (" << j << "," << i << ")" << endl;
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
						{
							cout << "Neuron " << neuronID << " post " << pst_n << " substracting " << LTD[m] << " from synapse (" << neuronID << "," << pst_n << ")" << endl;
							sd[neuronID][m] -= LTD[m];
							cout << "Resulting sd for this synapse: " << sd[neuronID][m] << endl;	
						}

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

				LTP[i][t+D+1]=0.98*LTP[i][t+D];
				LTD[i]*=0.95;
			}
	#if defined(WATCH_NEURONS) || defined(WATCH_DERIVATIVES)
			for (unsigned int i=0; i < std::vector<unsigned int>(WATCHED_NEURONS).size(); i++)
			{
	#ifdef WATCH_NEURONS
				watched_membrane[i].push_back(membranes[std::vector<unsigned int>(WATCHED_NEURONS)[i]]);
	#ifdef WATCH_ADAPTATION
				watched_us[i].push_back(u[std::vector<unsigned int>(WATCHED_NEURONS)[i]]);
	#endif
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

		if (t==1000)
		{

#ifdef WITH_DA
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
#else
			for (unsigned int i=0;i<N;i++)		// prepare for the next sec
			{
				cout << "Neuron " << i << " num_post " << num_post[i] << endl;
				if (i< Ne && num_post[i] > 0)
				{
#ifdef DEBUG_OUTPUT
					cout << "Neuron " << i << "'s weights: ";
#endif
					for (unsigned int j=0;j<num_post[i];j++)
					{
						weights[i][j] += 0.01+sd[i][j];
						if (weights[i][j]>sm) weights[i][j]=sm;
						if (weights[i][j]<0) weights[i][j]=0.0;
#ifdef DEBUG_OUTPUT
						cout << " " << weights[i][j];
#endif
						sd[i][j] *= 0.9;
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
#endif

			for (unsigned int i=0;i<N;i++)		// prepare for the next sec
			{
				for (unsigned int j=0;j<D+1;j++)
				{
//#ifdef WITH_DA
					//weights[i][j] += 0.01+sd[i][j];
//#endif
					LTP[i][j] = LTP[i][1000+j];
				}
			}


			//STOP_TIMER("one second loop",loops)

			cout << "Second " << sec << ": found " << total_spikes << endl;
			cout << "Second " << sec << ": firing rate " << total_spikes/N << endl;

			write_spikes("spikes_"+name+".txt",sec,spikes,k);

#ifdef WATCH_NEURONS
			write_watched_membranes("membrane_"+name+".txt",sec,watched_membrane,watched_us,WATCHED_NEURONS);
#endif
#ifdef WATCH_DERIVATIVES
			write_derivatives("stdp_"+name+".txt",sec,watched_LTP,watched_LTD,WATCHED_NEURONS);
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

void Network::init_random()
{

	srand(42);
	// This vector stores all pre-synaptic delays temporarily
	vector<unsigned int> pre_delays[N][D];

	for (int i=0; i<N; i++)
	{
		delay_start[i] = new unsigned int[D];
		delay_start_pre[i] = new unsigned int[D];
		delay_count[i] = new unsigned int[D];
		delay_count_pre[i] = new unsigned int[D];

		for (int d=0; d<D; ++d)
			delay_count_pre[i][d] = 0;
		num_pre[i] = 0;
	}

	/*
	 * init neuron parameters
	 */
	// iterate over neurons
	for (int i=0; i<N; i++) {
		I[i] = new float[D+1];
		LTP[i] = new float[1001+D];
		//membranes[i] = (float)rand()/float(RAND_MAX)*50.0;
		membranes[i] = v_reset;

		for (int del=0; del<=D; del++)
		{
			I[i][del] = 0;//(float)rand()/float(RAND_MAX)*10.0;
		}
		u[i] = 0.2*membranes[i];
		if (i < Ne) {
			a[i] = 0.02;
			d[i] = 8.0;
		}
		else {
			a[i] = 0.1;
			d[i] = 2.0;
		}
	}

	/*
	 * init connectivity
	 */
	// iterate over neurons
	for (int i=0; i<N; i++) {
		// pick how many post-synaptic connections this neuron will have
		num_post[i] = getrandom(M);
		// create temporary objects for delays, post-synaptic neurons and weights
		int delays[M];
		// the array synapse_numbers holds all synapse numbers for a given delay
		vector<int> synapse_numbers[D];
		vector<int> pst_n;
		vector<float> wgt_n;
		// reset STDP variables
		LTD[i] = 0.0;
		for (unsigned int j=0;j<1000+D;++j)
			LTP[i][j] = 0.0;
		// set delay_count and delay_start initially to zero
		for (int d=0; d<D; d++)
		{
			delay_start[i][d] = 0;
			delay_count[i][d] = 0;
		}
		// create values for connectivity
		// This will set weights and pick delays for each
		// of the neuron's post-synaptic connections
		for (int m=0; m<num_post[i]; m++)
		{
			// Initialize weights to fixed values
			if (i<Ne)
				wgt_n.push_back(6.);
			else
				wgt_n.push_back(-5.);

			// set weight derivatives initially to zero
			sd[i].push_back(0.0);
			// pick a delay and sort it into the delay_count array
			delays[m] = getrandom(D);
			synapse_numbers[delays[m]].push_back(m);
			delay_count[i][delays[m]]++;
			// pick a post-synaptic neuron number
			unsigned int tmp_post = getrandom(N);
			// make sure the post-syn. neuron is not the pre-syn. neuron
			while (tmp_post == i)
				tmp_post = getrandom(N);

			pst_n.push_back(tmp_post);

#ifdef DEBUG_OUTPUT
			cout << "Neuron " << i << ": Added post-neuron " << pst_n.back() << " with delay " << delays[m] << endl;
#endif

			// update pre-synaptic connection arrays
			// add the current neuron to the pre_neuron list
			// of the previously chosen post-synaptic neuron
			pre_delays[tmp_post][delays[m]].push_back(num_pre[tmp_post]);
			delay_count_pre[tmp_post][delays[m]]++;
			//pre_delays[tmp_post].push_back(delays[m]);
			pre_neurons[tmp_post].push_back(i);
			++num_pre[tmp_post];
		}
		// sort connectivity values according to delay values
		vector<unsigned int> pst_n_sort;
		vector<float> wgt_n_sort;
		for (int d=0; d<D; d++)
		{
			vector<int>::iterator delay_it = synapse_numbers[d].begin();
			while (delay_it != synapse_numbers[d].end())
			{
				pst_n_sort.push_back(pst_n[*delay_it]);
				wgt_n_sort.push_back(wgt_n[*delay_it]);
				/*
				   actually, sd would have to be sorted here, too!
				   But since it's initialized with zeros anyway,
				   we can just skip it here
			    */

				delay_it++;
			}
		}

		// add values to global arrays
		weights[i] = wgt_n_sort;
		post_neurons[i] = pst_n_sort;
		/*
		   The next loop calculates the correct values for
		   the delay_count and delay_start arrays, according
		   to Nageswaran et al.
		*/
		int last_delay = 0;
		for (int d=0; d<D; d++)
		{
			if (delay_count[i][d] > 0)
				delay_start[i][d] = last_delay;
			else
				delay_start[i][d] = 0;
			last_delay += delay_count[i][d];
		}
	}
	// iterate over neurons again
	// re-arrange pre-syn. information
	for (int i=0; i<N; i++)
	{
		// get pointers to the derivatives of the weights
		for (int m=0; m<num_post[i]; m++)
			sd_pre[post_neurons[i][m]].push_back(&sd[i][m]);// pointer to the derivative

		// sort connectivity values according to delay values
		vector<unsigned int> pre_n_sort;
#ifdef DEBUG_OUTPUT
		cout << "Neuron " << i << endl;
		cout << "post_neurons: " << post_neurons[i].size() << endl;
		cout << "pre_neurons: " << pre_neurons[i].size() << endl;
#endif
		if (!pre_neurons[i].empty())
		{
			for (int d=0; d<D; d++)
			{
				//cout << "pre_delays[" << d << "]";
				//print_loop(pre_delays[i][d]);
				for (unsigned int delay_it : pre_delays[i][d])
				{
					pre_n_sort.push_back(pre_neurons[i][delay_it]);
				}
			}
			pre_neurons[i] = pre_n_sort;
		}

		int last_delay = 0;
		for (int d=0; d<D; d++)
		{
			if (delay_count_pre[i][d] > 0)
				delay_start_pre[i][d] = last_delay;
			else
				delay_start_pre[i][d] = 0;
			last_delay += delay_count_pre[i][d];
		}
	}
}

