ortools = dependency('ortools',required: false)

benchmark_args = []
if ortools.found()
    benchmark_args += ['-DUSE_ORTOOLS']
endif

executable('benchmark-mcf','benchmark-mcf.cpp',
    dependencies: [dep_mincostflow,ortools],cpp_args: benchmark_args)
