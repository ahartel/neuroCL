#include "helpers.h"
#include <iostream>
#include <fstream>

using namespace std;

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

void write_derivatives(string filename, unsigned int sec, std::vector<std::vector<float> >const& watched_LTP, std::vector<std::vector<float> > const& watched_LTD)
{
	filename = string("./results/")+filename;
	if (sec==0)
	{
		remove( filename.c_str() );
		ofstream myfile( filename, ios::app);
		myfile << "Time";
		for (auto n : neurons_tobe_watched)
		{
			myfile << ",Neuron " << n << " LTP";
			myfile << ",Neuron " << n << " LTD";
		}
		myfile << endl;
		myfile.close();
	}

	ofstream myfile(filename, ios::app);
	for (int t=0; t<1000; t++) {
		myfile << sec*1000+t;
		for (unsigned int n=0; n < num_watched_neurons; n++)
		{
			if (watched_LTP[n][t] > 0)
				myfile << "," << watched_LTP[n][t];
			else
				myfile << ",0.0";
			if (watched_LTD[n][t] > 0)
				myfile << "," << watched_LTD[n][t];
			else
				myfile << ",0.0";
		}
		myfile << endl;
	}
	myfile.close();
}

void write_spikes(string filename, unsigned int sec, vector<unsigned int> const& spikes, unsigned int* k)
{
	filename = string("./results/")+filename;
	if (sec==0)
	{
		ofstream myfile(filename);
		myfile << "time [ms], neuron no." << endl;
		myfile.close();
	}

	ofstream myfile(filename, ios::app);
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

void write_watched_membranes(string filename, unsigned int sec, std::vector<std::vector<float> >const& watched_membrane, std::vector<std::vector<float> > const& watched_us)
{
	filename = string("./results/")+filename;

	if (sec==0)
	{
		remove( filename.c_str() );
		ofstream myfile(filename,ios::app);
		myfile << "Time";
		for (auto n : neurons_tobe_watched)
		{
			myfile << ",Neuron " << n;
		}
		myfile << endl;
		myfile.close();
	}

	ofstream myfile(filename,ios::app);
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

void init_neurons_sparse(
	vector<float> & membranes,
	vector<float> & u,
	vector<float> & d,
	vector<float> & a,
	float I[N][D+1],
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

	for (int i=0; i<N; i++)
	{
		for (int d=0; d<D; ++d)
			delay_count_pre[i][d] = 0;
		num_pre[i] = 0;
	}

	/*
	 * init neuron parameters
	 */
	// iterate over neurons
	for (int i=0; i<N; i++) {
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

void init_neurons(float* membranes, float* u, float* d, float* a, float I[N][D], float weights[N][M], int delays[N][M], int post_neurons[N][M], int* num_post)
{
	srand(42);
	for (int i=0; i<N; i++) {
		//membranes[i] = (float)rand()/float(RAND_MAX)*50.0;
		membranes[i] = v_reset;

		for (int del=0; del<D; del++)
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