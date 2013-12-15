#include <vector>

#include "compare/network.h"
#include "compare/network_description.h"

using namespace std;

int main()
{
	NetworkDescription n_spec;

	vector<unsigned int> exc_neurons;
	vector<unsigned int> inh_neurons;
	vector<unsigned int> inp_neurons;

	// inhibitory population
	for (size_t i=0; i<5; ++i)
		inh_neurons.push_back(n_spec.add_neuron());

	// excitatory population
	for (size_t i=0; i<10; ++i)
	{
		unsigned int neuron = n_spec.add_neuron();
		exc_neurons.push_back(neuron);
		for (auto inh : inh_neurons)
		{
			n_spec.connect_neurons(neuron,inh,6,0,0);
			n_spec.connect_neurons(inh,neuron,-5,0,0);
		}
	}

	// input population
	for (size_t i=0; i<10; ++i)
	{
		unsigned int neuron = n_spec.add_neuron();
		inp_neurons.push_back(neuron);
		for (auto exc : exc_neurons)
		{
			n_spec.connect_neurons(neuron,exc,6,0,1);
		}
	}
	

	Network n("wta",n_spec);
	//n.enable_random_background();

	INIT_TIMER(loops)

/*
	for (unsigned int t=0;t<1000;t++)
	{
		if (t%30 >=0 && t%30 <6)
		{
			n.add_current({{t%30,0,20.}});
			n.add_current({{t%30,0,20.}});
		}
		n.step();
	}
*/

	for (unsigned int t=0;t<50000;t++)
	{
		n.add_current({{15+getrandom(10),0,20.}});
		//if (rand()%2)
		//	n.add_current({{6,getrandom(3),20.}});
		n.step();
	}

	for (unsigned int t=0;t<1000;t++)
	{
		if (t%30 >=0 && t%30 <6)
		{
			n.add_current({{t%30,0,20.}});
			n.add_current({{t%30,0,20.}});
		}
		n.step();
	}

	for (unsigned int t=0;t<1000;t++)
	{
		n.add_current({{getrandom(6),0,20.}});
		n.step();
	}

	STOP_TIMER("all loops",loops)
}
