#pragma once
#include <string>
#include "helpers.h"
#include "network_description.h"
#include <vector>

using namespace std;

struct pre_spike
{
	unsigned int neuron;
	unsigned int delay;
	float weight;

	pre_spike(unsigned int n, unsigned int d, float w)
	{
		neuron = n;
		delay = d;
		weight = w;
	};
};

class Network
{

public:
	Network(unsigned int e, unsigned int i,unsigned int m);
	Network(std::string const& n, NetworkDescription const&);
	Network(std::string const& n, unsigned int e, unsigned int i, unsigned int m);
	~Network();
	void step();
	void add_current(vector<pre_spike> const&);
	std::vector<unsigned int> get_last_spikes();
#ifdef WITH_DA
	void eject_dopamine();
#endif
	void enable_random_background();

private:
	//network parameters
	unsigned int Ne,Ni,N;
	unsigned int M; // number of postsynaptic neurons
	constexpr static float v_thresh = 30.0/*mV*/;
	constexpr static float v_reset = -65.0/*mV*/;
	constexpr static float h = 1; // timestep in milliseconds
	constexpr static unsigned int D = 20; // max. delay in ms
	constexpr static float sm = 10.0;		// maximal synaptic strength

	int get_delays(int delay_start[D],int synapseID);
	void init();
	void init_from_spec(NetworkDescription const&);
	void init_random();


	string name;

	unsigned int sec;
	unsigned int last_sec;
	unsigned int t;
	unsigned int total_spikes;
	unsigned int last_total_spikes;
	int delay_pointer;

	// neuron variables
	vector<float> membranes;
	vector<float> u;
	vector<float> d;
	vector<float> a;
	// input current
	float** I; //[N][D+1];
	// network parameters
	unsigned int** delay_start; //[N][D];
	unsigned int** delay_start_pre;
	unsigned int** delay_count; //[N][D];
	unsigned int** delay_count_pre;

	vector<float>* weights;//[N];
	vector<unsigned int>* post_neurons;//[N];
	vector<unsigned int>* pre_neurons;//[N];
	unsigned int* num_post;//[N];
	unsigned int* num_pre;//[N];
	// STDP functions
	float** LTP;//[N][1001+D]
	float* LTD;//[N];
	vector< vector<float> > sd;		// matrix of synaptic weights and their derivatives
	vector< vector<float*> > sd_pre;		// presynaptic weights
	// spike storage
	vector<unsigned int> spikes;
	unsigned int k[1000];
#ifdef WITH_DA
	float DA;
#endif

	bool background_spikes;

	std::vector<std::vector<float> > watched_LTP;
	std::vector<std::vector<float> > watched_LTD;
	std::vector<std::vector<float> > watched_membrane;
	std::vector<std::vector<float> > watched_us;
};
