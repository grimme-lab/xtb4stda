# This file is part of xtb4stda.
#
# Copyright (C) 2019-2020 Sebastian Ehlert
#
# xtb4stda is free software: you can redistribute it and/or modify it under
# the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# xtb4stda is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with xtb4stda.  If not, see <https://www.gnu.org/licenses/>.

project('xtb4stda', 'fortran',
        version: '1.0',
        license: 'LGPL3',
        meson_version: '>=0.49')

if get_option('interface') == '64' and get_option('la_backend') != 'mkl'
  error('64 bit integer interface only supported for MKL backend')
endif

fc = meson.get_compiler('fortran')

if fc.get_id() == 'gcc'
  add_project_arguments('-ffree-line-length-none', language: 'fortran')
  add_project_arguments('-cpp', language: 'fortran')
  add_project_arguments('-fbacktrace', language: 'fortran')
  if get_option('interface') == '64'
    add_project_arguments('-fdefault-integer-8', language: 'fortran')
  endif
elif fc.get_id() == 'intel'
  add_project_arguments('-axAVX2',    language: 'fortran')
  add_project_arguments('-fpp',    language: 'fortran')
  add_project_arguments('-traceback', language: 'fortran')
  if get_option('static')
    add_project_link_arguments('-static', language: 'fortran')
  endif
  if get_option('interface') == '64'
    add_project_arguments('-i8', language: 'fortran')
  endif
elif fc.get_id() == 'intel-cl' # Windows with Intel Fortran
  add_project_arguments('-QaxCORE-AVX2',    language: 'fortran')
  add_project_arguments('-fpp',    language: 'fortran')
  add_project_arguments('-traceback', language: 'fortran')
  if get_option('static')
    add_project_arguments('-MT', language: 'fortran') # link to static Fortran library
  endif
  if get_option('interface') == '64'
    add_project_arguments('-4I8', language: 'fortran')
  endif
endif

dependencies = []

la_backend = get_option('la_backend')
if la_backend == 'mkl'
  if host_machine.system() == 'windows'
    libmkl = [fc.find_library('libiomp5md')]
  else
    libmkl = [fc.find_library('pthread')]
    libmkl += fc.find_library('m')
    libmkl += fc.find_library('dl')
  endif
  if (fc.get_id() == 'intel') or (fc.get_id() == 'intel-cl')
    if get_option('interface') == '64'
      libmkl += fc.find_library('mkl_intel_ilp64')
    else
      libmkl += fc.find_library('mkl_intel_lp64')
    endif
    libmkl += fc.find_library('mkl_intel_thread')
  else
    if get_option('interface') == '64'
      libmkl += fc.find_library('mkl_gf_ilp64')
    else
      libmkl += fc.find_library('mkl_gf_lp64')
    endif
    libmkl += fc.find_library('mkl_gnu_thread')
  endif
  libmkl += fc.find_library('mkl_core')
  if host_machine.system() != 'windows'
    libmkl += fc.find_library('iomp5')
  endif
  dependencies += libmkl
elif la_backend == 'openblas'
  dependencies += fc.find_library('openblas', required : true)
  dependencies += fc.find_library('lapack', required : true)
elif la_backend == 'custom'
  foreach lib: get_option('custom_libraries')
    dependencies += fc.find_library(lib)
  endforeach
else
  dependencies += fc.find_library('blas', required : true)
  dependencies += fc.find_library('lapack', required : true)
endif

if get_option('openmp')
  if fc.get_id() == 'intel'
    add_project_arguments('-qopenmp', language : 'fortran')
    add_project_link_arguments('-qopenmp', language : 'fortran')
  elif fc.get_id() == 'intel-cl'
    add_project_arguments('-Qopenmp', language : 'fortran')
  else
    add_project_arguments('-fopenmp', language : 'fortran')
    add_project_link_arguments('-fopenmp', language : 'fortran')
  endif
endif

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
