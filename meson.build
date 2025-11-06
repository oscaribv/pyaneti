project('pyaneti', 'c','fortran',
  version: '2',
  meson_version: '>=1.0.1',
)

add_languages('fortran')

# python helper: import the python module and find the installation object
py = import('python')
py_inst = py.find_installation(pure : false)
py_dep = py_inst.dependency()

# program to run f2py
pyprog = find_program('python3', required : true)

# get numpy / f2py include dirs at configure-time
numpy_get_include = run_command(pyprog, ['-c', 'import numpy; print(numpy.get_include())']).stdout().strip()
f2py_get_include = run_command(pyprog, ['-c', 'import numpy.f2py; print(numpy.f2py.get_include())']).stdout().strip()
inc_dirs = include_directories(numpy_get_include, f2py_get_include)

# Fortran sources 
fortran_files = [
  'src/constants.f90',
  'src/todo.f90',
  'src/qpower2.f90',
  'src/quad.f90',
  'src/ftr.f90',
  'src/frv.f90',
  'src/bayesian.f90',
  'src/matrices.f90',
  'src/kernels.f90',
  'src/mcmc.f90',
  'src/multi-gp-routines.f90'
]

pyaneti_link_args = ['-llapack', '-lblas','-fopenmp']

f2py_cmd = [ pyprog, '-m', 'numpy.f2py', '@INPUT@', '-m', 'pyaneti', '--lower' ]

# Path to fortranobject.c in numpy.f2py include dir
# f2py_get_include is the run_command(...) stdout from earlier
fortran_obj_path = f2py_get_include + '/fortranobject.c'
# Create a meson file object for it
fortran_object = files(fortran_obj_path)

# custom_target f2py_gen stays the same
f2py_gen = custom_target('pyaneti_f2py_gen',
  input : fortran_files,
  output : ['pyanetimodule.c','pyaneti-f2pywrappers.f', 'pyaneti-f2pywrappers2.f90'],  # use the names that match your run
  command : f2py_cmd,
)

# Now create the Python extension; include the fortranobject C source so PyFortran_Type is defined.
py_inst.extension_module('pyaneti',
  sources : fortran_files + [f2py_gen] + [fortran_object],
  include_directories : inc_dirs,
  dependencies : [ py_dep ],
  link_args : pyaneti_link_args,
  install : true
)

