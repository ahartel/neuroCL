
def options(opt):
	opt.load('compiler_cxx')

def configure(conf):
	conf.load('compiler_cxx')
	conf.check_cxx(lib='OpenCL',use='opencl')
	conf.env.CXXFLAGS = ['--std=c++11','-g']

def build(bld):
	bld.objects(
		name = 'neuroobj',
		source=[
			'source/errorMessage.cpp',
		]
	)
	bld.program(source='source/neurosim.cpp', target='neurosim', use='neuroobj OPENCL')

