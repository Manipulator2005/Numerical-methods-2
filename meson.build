project('chebyshev_methods', 'cpp',
  default_options : ['warning_level=3', 'cpp_std=c++20']
)

project_inc = include_directories('include')
build_type = get_option('buildtype')
iter_src = ['src/chebyshev_one_step.cpp', 'src/chebyshev_two_step.cpp']

if build_type.startswith('release')
  add_project_arguments('-march=native', '-UDEBUG', language: 'cpp')
  lib = shared_library('iter_methods', 
  iter_src, 
  include_directories: project_inc)

  exe = executable('prog', 
    'main.cpp', 
    include_directories: project_inc,
    link_with: lib,
    build_rpath: meson.current_build_dir())
endif