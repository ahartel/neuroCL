#include <stdlib.h>
#include <iostream>
#include <vector>
/*
 * Constants
 */
const int Ne = 8;//12800;
const int Ni = 2;//3200;
const int N = Ne+Ni; // total number of neurons
const static int M = 2; // number of postsynaptic neurons
const float v_thresh = 30./*mV*/;
const float v_reset = -65/*mV*/;
const int T = 10; // number of seconds to simulate
const float h = 1; // timestep in milliseconds
const static int D = 20; // max. delay in ms
const float	sm = 10.0;		// maximal synaptic strength

using namespace std;

#define getrandom(max1) ((rand()%(int)((max1)))) // random integer between 0 and max-1

#define TIMING
#undef WATCH_NEURONS
#undef WATCH_ADAPTATION

template<typename T>
void print_loop (T array, size_t size)
{
	if (size > 0)
	{
		if (array[0] < 10)
			cout << 0;
		cout << array[0];
		for (unsigned int i=1; i<size; i++)
		{
			cout << ", ";
			if (array[i] < 10)
				cout << 0;
			cout << array[i];
		}
	}
	cout << endl;
}

template<typename T>
void print_loop (T array)
{
	print_loop(array,array.size());
}

void init_neurons_sparse(
	float* membranes,
	float* u,
	float* d,
	float* a,
	float I[N][D],
	vector<float> weights[N],
	unsigned int delay_start[N][D],
	unsigned int delay_count[N][D],
	vector<unsigned int> post_neurons[N],
	unsigned int num_post[N],
	float LTP[N][1001+D],
	float LTD[N],
	vector< vector<float> >& sd,
	vector< vector<float*> >& sd_pre,
	unsigned int num_pre[N],
	vector<unsigned int> pre_neurons[N],
	unsigned int delay_start_pre[N][D],
	unsigned int delay_count_pre[N][D]
	)
{

	srand(42);
	// This vector stores all pre-synaptic delays temporarily
	vector<unsigned int> pre_delays[N][D];

	/*
	 * init neuron parameters
	 */
	// iterate over neurons
	for (int i=0; i<N; i++) {
		//membranes[i] = (float)rand()/float(RAND_MAX)*50.0;
		membranes[i] = v_reset;

		for (int d=0; d<D; d++)
		{
			I[i][d] = 0;//(float)rand()/float(RAND_MAX)*10.0;
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
			// pick a post-synaptic neuron
			unsigned int tmp_post = getrandom(N);
			// make sure the post-syn. neuron is not the pre-syn. neuron
			while (tmp_post == i)
				tmp_post = getrandom(N);

			pst_n.push_back(tmp_post);

			cout << "Neuron " << i << ": Added post-neuron " << pst_n.back() << " with delay " << delays[m] << endl;

			// update pre-synaptic connection arrays
			// add the current neuron to the pre_neuron list
			// of the previously chosen post-synaptic neuron
			pre_delays[tmp_post][delays[m]].push_back(num_pre[tmp_post]);
			//pre_delays[tmp_post].push_back(delays[m]);
			sd_pre[tmp_post].push_back(&sd[i][m]);// pointer to the derivative
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

				delay_it++;
			}
		}

		// add values to global arrays
		weights[i] = wgt_n_sort;
		post_neurons[i] = pst_n_sort;

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
		cout << "Neuron " << i << endl;
		// sort connectivity values according to delay values
		vector<unsigned int> pre_n_sort;
		cout << "post_neurons: " << post_neurons[i].size() << endl;
		cout << "pre_neurons: " << pre_neurons[i].size() << endl;
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

void init_neurons(float* membranes, float* u, float* d, float* a, float I[N][D], float weights[N][M], int delays[N][M], int post_neurons[N][M], int* num_post)
{
	srand(42);
	for (int i=0; i<N; i++) {
		//membranes[i] = (float)rand()/float(RAND_MAX)*50.0;
		membranes[i] = v_reset;

		for (int d=0; d<D; d++)
		{
			I[i][d] = 0;//(float)rand()/float(RAND_MAX)*10.0;
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

		num_post[i] = getrandom(M);

		for (int m=0; m<num_post[i]; m++)
		{
			if (i<Ne)
				weights[i][m] = 6.;
			else
				weights[i][m] = -5.;

			delays[i][m] = getrandom(D);

			post_neurons[i][m] = getrandom(N);
			while (post_neurons[i][m] == i)
				post_neurons[i][m] = getrandom(N);
		}
	}
}

#include <chrono>
#include <iostream>

#ifdef TIMING
#define INIT_TIMER(var) auto var = std::chrono::high_resolution_clock::now();
#define START_TIMER(var)  var = std::chrono::high_resolution_clock::now();
#define STOP_TIMER(name,var)  std::cout << "RUNTIME of " << name << ": " << \
std::chrono::duration_cast<std::chrono::microseconds>( \
		std::chrono::high_resolution_clock::now()-var \
).count() << " us " << std::endl;
#else
#define INIT_TIMER(var)
#define START_TIMER(var)
#define STOP_TIMER(name,var)
#endif

const unsigned int neurons_tobe_watched[] = {0,17};
const unsigned int num_watched_neurons = 2;

