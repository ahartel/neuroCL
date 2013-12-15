#include <vector>

#include "compare/network.h"
#include "compare/network_description.h"

int main()
{
	NetworkDescription n_spec;

	unsigned int n00 = n_spec.add_neuron();
	unsigned int n01 = n_spec.add_neuron();
	unsigned int n02 = n_spec.add_neuron();
	unsigned int n10 = n_spec.add_neuron();
	unsigned int n11 = n_spec.add_neuron();
	unsigned int n12 = n_spec.add_neuron();
	unsigned int n20 = n_spec.add_neuron();
	unsigned int n21 = n_spec.add_neuron();
	unsigned int n22 = n_spec.add_neuron();
	unsigned int target = n_spec.add_neuron();

	n_spec.connect_neurons(n00,target,6,0,1);
	n_spec.connect_neurons(n01,target,6,0,1);
	n_spec.connect_neurons(n02,target,6,0,1);
	n_spec.connect_neurons(n10,target,6,0,1);
	n_spec.connect_neurons(n11,target,6,0,1);
	n_spec.connect_neurons(n12,target,6,0,1);
	n_spec.connect_neurons(n20,target,6,0,1);
	n_spec.connect_neurons(n21,target,6,0,1);
	n_spec.connect_neurons(n22,target,6,0,1);

	Network n("stdp_test",n_spec);
	n.enable_random_background();

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
		n.add_current({{t%9,0,20.}});
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
