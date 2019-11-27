project('xtb4stda', 'fortran', 'c',
        version: '1.0',
        license: 'LGPL3',
        meson_version: '>=0.49')

fc = meson.get_compiler('fortran')

if fc.get_id() == 'gcc'
  add_project_arguments('-ffree-line-length-none', language: 'fortran')
  add_project_arguments('-fbacktrace', language: 'fortran')
elif fc.get_id() == 'intel'
  add_project_arguments('-axAVX2',    language: 'fortran')
  add_project_arguments('-traceback', language: 'fortran')
  add_project_link_arguments('-static', language: 'fortran')
  add_project_link_arguments('-static', language: 'c') # icc will do linking
endif

dependencies = []
dependencies_exe = []
dependencies_sha = []

libmkl = [fc.find_library('pthread')]
libmkl += fc.find_library('m')
libmkl += fc.find_library('dl')
if fc.get_id() == 'gcc'
  libmkl += fc.find_library('mkl_gf_lp64')
  libmkl += fc.find_library('mkl_gnu_thread')
elif fc.get_id() == 'intel'
  libmkl += fc.find_library('mkl_intel_lp64')
  libmkl += fc.find_library('mkl_intel_thread')
endif
libmkl += fc.find_library('mkl_core')
libmkl += fc.find_library('iomp5')
dependencies += libmkl

if fc.get_id() == 'gcc'
  add_project_arguments('-fopenmp', language: 'fortran')
  add_project_link_arguments('-fopenmp', language: 'fortran')
elif fc.get_id() == 'intel'
  add_project_arguments('-qopenmp', language: 'fortran')
  add_project_link_arguments('-qopenmp', language: 'fortran')
endif
dependencies += dependency('threads')

srcs = []
subdir('src')

xtb4stda_exe = executable(meson.project_name(), srcs,
                          dependencies: dependencies,
                          include_directories: include_directories('include'),
                          install: true)

parameter_files = ['.param_stda1.xtb',
                   '.param_stda2.xtb',
                   '.param_gbsa_acetone',
                   '.param_gbsa_acetonitrile',
                   '.param_gbsa_benzene',
                   '.param_gbsa_ch2cl2',
                   '.param_gbsa_chcl3',
                   '.param_gbsa_cs2',
                   '.param_gbsa_dmso',
                   '.param_gbsa_ether',
                   '.param_gbsa_h2o',
                   '.param_gbsa_methanol',
                   '.param_gbsa_thf',
                   '.param_gbsa_toluene']

install_data(parameter_files)
