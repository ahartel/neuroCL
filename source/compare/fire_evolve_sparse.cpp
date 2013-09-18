#include <stdlib.h>
#include <stdexcept>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include "helpers.h"
#include "network.h"

using namespace std;


int main()
{
	Network net(800,200,100);

	for (unsigned int sec=0; sec<100; sec++)
	{
		INIT_TIMER(loops)
		for (unsigned int t=0;t<1000;t++)
		{
			net.step();
		}
		STOP_TIMER("one second loop",loops)
	}
}



