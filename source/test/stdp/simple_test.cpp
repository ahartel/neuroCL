#include <vector>

#include "compare/network.h"
#include "compare/network_description.h"

int main()
{
	NetworkDescription n_spec;

	unsigned int n0 = n_spec.add_neuron();
	unsigned int n1 = n_spec.add_neuron();
	unsigned int n2 = n_spec.add_neuron();
	unsigned int target = n_spec.add_neuron();

	n_spec.connect_neurons(n0,target,6,0);
	n_spec.connect_neurons(n1,target,6,0);
	n_spec.connect_neurons(n2,target,6,0);

	Network n("stdp_test",n_spec);

	INIT_TIMER(loops)

	for (unsigned int t=0;t<2000;t++)
	{
		n.add_current({{t%3,0,20.}});
		if (rand()%2)
			n.add_current({{3,getrandom(3),20.}});
		n.step();
	}

	for (unsigned int t=0;t<1000;t++)
	{
		if (!(t%20))
		{
			n.add_current({{2,0,20.}});
			n.add_current({{1,0,20.}});
			n.add_current({{0,0,20.}});
		}
		//if (rand()%2)
		//	n.add_current({{3,getrandom(3),20.}});
		n.step();
	}

	for (unsigned int t=0;t<1000;t++)
	{
		n.add_current({{t%3,0,20.}});
		//if (rand()%2)
		//	n.add_current({{3,getrandom(3),20.}});
		n.step();
	}

	STOP_TIMER("all loops",loops)
}
