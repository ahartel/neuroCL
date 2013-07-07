#pragma once

struct pre_spike
{
	unsigned int neuron;
	unsigned int delay;
	unsigned int weight;
};

class Network
{

public:
	Network();
	void step();
	void add_spikes(vector<pre_spike> const&);

private:
	int get_delays(int delay_start[D],int synapseID);

	unsigned int sec;
	unsigned int t;
	unsigned int total_spikes;
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
