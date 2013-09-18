#pragma once

#include <list>
#include <map>
#include <tuple>

class Connection;

class Neuron
{
public:
	Neuron();
	Neuron(unsigned int);
	unsigned int getId() const { return id; };
	unsigned int add_connection(Neuron*, float, unsigned int);
	unsigned int getNumberConnections() const { return connection_index; }
	std::tuple<float,unsigned int> getConnection(unsigned int i) const;
	unsigned int getPostNeuron(unsigned int i) const;
private:
	unsigned int id;
	unsigned int connection_index;
	std::map<unsigned int,Connection> post_neurons;
};

class Connection
{
public:
	Connection();
	Connection(Neuron*,float, unsigned int);
	std::tuple<float,unsigned int> data() const { return std::make_tuple(weight,delay); };
	unsigned int getNeuronId() const { return post_neuron->getId(); }

private:
		unsigned int delay;
		Neuron* post_neuron;
		float weight;
};


class NetworkDescription
{
public:
	NetworkDescription();
	unsigned int add_neuron();
	unsigned int connect_neurons(unsigned int, unsigned int, float, unsigned int);
	unsigned int getNumberOfNeurons() const { return neuron_index; };
	unsigned int getNumberConnectionsOfNeuron(unsigned int i) const { return neurons.at(i).getNumberConnections(); };
	std::tuple<float,unsigned int> getConnectionOfNeuron(unsigned int i, unsigned int j) const { return neurons.at(i).getConnection(j); };
	unsigned int getPostOfNeuron(unsigned int i, unsigned int j) const { return neurons.at(i).getPostNeuron(j); };

private:
	std::map<unsigned int,Neuron> neurons;
	unsigned int neuron_index;
};
