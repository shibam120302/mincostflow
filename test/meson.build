test('Vectorized Map',
    executable('vec-map','vec-map.cpp',
        dependencies: dep_mincostflow))

test('Shortest Path',
    executable('short-path','short-path.cpp',
        dependencies: dep_mincostflow))

test('Max Flow',
    executable('max-flow','max-flow.cpp',
        dependencies: dep_mincostflow))

test('Min Cost Max Flow',
    executable('min-cost-flow','min-cost-flow.cpp',
        dependencies: dep_mincostflow))

python=find_program('python')
python_env = environment()

python_env.set('PYTHONPATH',meson.project_source_root()+'/python')
python_env.set('LD_LIBRARY_PATH',meson.project_build_root())

python_test = files('python_api.py')
test('Python API',
    python,
    args: python_test,
    env: python_env)
