#include <string>
#include <CL/cl.h>

using namespace std;

std::string errorMessage(cl_int error);
void check_error(std::string prefix, cl_int error);
