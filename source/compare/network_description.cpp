#include "network_description.h"

Connection::Connection()
{
}

Connection::Connection(Neuron* n, float w, unsigned int d) :
	post_neuron(n),weight(w),delay(d)
{
}


unsigned int Neuron::getPostNeuron(unsigned int i) const { return post_neurons.at(i).getNeuronId(); }
std::tuple<float,unsigned int> Neuron::getConnection(unsigned int i) const { return post_neurons.at(i).data(); }

Neuron::Neuron() :
	connection_index(0)
{
}

Neuron::Neuron(unsigned int i) :
	connection_index(0)
{
	id = i;
}

unsigned int Neuron::add_connection(Neuron* post, float weight, unsigned int delay)
{
	post_neurons[connection_index] = Connection(post,weight,delay);
	return connection_index++;
}

NetworkDescription::NetworkDescription() :
	neuron_index(0)
{

}

unsigned int NetworkDescription::add_neuron()
{
	neurons[neuron_index] = Neuron(neuron_index);
	return neuron_index++;
}

unsigned int NetworkDescription::connect_neurons(unsigned int pre, unsigned int post, float weight, unsigned int delay)
{
	return neurons[pre].add_connection(&(neurons[post]), weight, delay);
}
