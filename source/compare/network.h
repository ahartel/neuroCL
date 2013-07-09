#pragma once
#include <string>
#include "helpers.h"
#include <vector>

using namespace std;

struct pre_spike
{
	unsigned int neuron;
	unsigned int delay;
	unsigned int weight;

	pre_spike(unsigned int n, unsigned int d, unsigned int w)
	{
		neuron = n;
		delay = d;
		weight = w;
	};
};

class Network
{

public:
	Network();
	Network(std::string const& n);
	void step();
	void add_spikes(vector<pre_spike> const&);
	std::vector<unsigned int> get_last_spikes();

private:
	int get_delays(int delay_start[D],int synapseID);
	void init();

	string name;

	unsigned int sec;
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
	float I[N][D+1];
	// network parameters
	unsigned int delay_start[N][D];
	unsigned int delay_start_pre[N][D];
	unsigned int delay_count[N][D];
	unsigned int delay_count_pre[N][D];
	vector<float> weights[N];
	vector<unsigned int> post_neurons[N];
	vector<unsigned int> pre_neurons[N];
	unsigned int num_post[N];
	unsigned int num_pre[N];
	// STDP functions
	float LTP[N][1001+D], LTD[N];
	vector< vector<float> > sd;		// matrix of synaptic weights and their derivatives
	vector< vector<float*> > sd_pre;		// presynaptic weights
	// spike storage
	vector<unsigned int> spikes;
	unsigned int k[1000];

};
