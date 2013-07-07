
def options(opt):
	opt.load('compiler_cxx')

def configure(conf):
	conf.load('compiler_cxx')
	conf.check_cxx(lib='OpenCL',use='opencl')
	conf.env.CXXFLAGS = ['-std=c++0x','-g']
	conf.env.INCLUDES = 'source/'

	conf.env.LIB_GL = ['glut','GL','GLU']

def build(bld):
	bld.program(
		target='Pong',
		source='source/Pong/main.cpp',
		use='GL',
		cxxflags='-std=c++0x',
	)

	bld.objects(
		name = 'neuroobj',
		source=[
			'source/errorMessage.cpp',
			'source/compare/network.cpp',
			'source/helpers.cpp',
		]
	)
	#bld.program(source='source/neurosim.cpp', target='neurosim', use='neuroobj OPENCL')
	bld.program(source='source/compare/fire_evolve.cpp', target='loop_fire_evolve', use='neuroobj')
	bld.program(source='source/compare/fire_evolve_sparse.cpp', target='fire_evolve_sparse', use='neuroobj')

