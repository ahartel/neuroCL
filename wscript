import waflib

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
		source=[
            'source/Pong/main.cpp',
			'source/compare/network.cpp',
        ],
		use=['GL','neuroobj'],
        includes='source',
		cxxflags='-std=c++0x -DWITH_DA',
		install_path='bin'
	)

	bld.objects(
		name = 'neuroobj',
		source=[
			'source/errorMessage.cpp',
			'source/compare/network_description.cpp',
			'source/helpers.cpp',
		],
	)

	#bld.program(source='source/neurosim.cpp', target='neurosim', use='neuroobj OPENCL')

	bld.program(
		source=['source/compare/fire_evolve.cpp'],
		target='loop_fire_evolve',
		use='neuroobj',
		install_path='bin'
	)
	bld.program(
		source=['source/compare/fire_evolve_sparse.cpp','source/compare/network.cpp'],
		target='fire_evolve_sparse',
		use='neuroobj',
		install_path='bin'
	)

	for file in bld.path.ant_glob('source/test/**/*.cpp'):
		target = waflib.Node.split_path(file.abspath())[-1].split('.')[0]
		bld.program(
			source=[file,'source/compare/network.cpp'],
			target=target,
			use='neuroobj',
			includes='source',
			cxxflags=['-DDEBUG_OUTPUT', '-DWATCH_DERIVATIVES','-DWATCHED_NEURONS={0,1,2,3}','-DWATCH_NEURONS','-DWATCH_WEIGHTS'],
			install_path='bin/'
		)

