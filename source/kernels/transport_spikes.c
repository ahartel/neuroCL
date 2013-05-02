
__kernel void evolve_neuron (
			__global const float* I,
			const int num,
			const float h
		)
{
	/* Step 2
	// transmit loop
	*/
	int last_spike_count;
	if (t > 0)
	{
		last_spike_count = k[t-1];
	}
	else
		last_spike_count = 0;

	for (int i=last_spike_count; i<total_spikes; i++)
	{
		int neuronID = spikes[i];
		//cout << "Found spike in neuron " << neuronID << endl;
		// skip this neuron, if it doesn't connect
		// (this can happen in a randomly connected network)
		if (num_post[neuronID] == 0)
			continue;
		// get first used delay value
		int d = 0;
		int delay_cnt = 0;
		// loop through delay values until a used one is found
		while (delay_count[neuronID][d] == 0 && d < D)
		{
			d++;
		}
		// save the number of neurons that use this delay
		delay_cnt = delay_count[neuronID][d];	
		// check if the loop has run through without result
		if (d == D-1 && delay_count[neuronID][d] == 0)
		{
			cout << "Neuron " << neuronID;
			throw "Error: No delays found.";
		}
		// Now, generate the current for all post-syn. neurons
		for (int m=0; m<num_post[neuronID]; m++)
		{
			I[post_neurons[neuronID][m]][(d+delay_pointer+1)%D] += weights[neuronID][m];
			// keep track of the number of remaining neurons with current delay
			// if necessary, increase delay value
			if (--delay_cnt == 0)
			{
				d++;
				while (delay_count[neuronID][d] == 0 && d < D)
				{
					d++;
				}
				delay_cnt = delay_count[neuronID][d];	
			}
		}
	}
}
