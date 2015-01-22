from distutils.core import setup, Extension

extension_mod = Extension("_breakout_detection",
                          sources = ["breakout_detection.i",
                                     "../src/edm-per.cpp",
                                     "../src/edm-multi.cpp",
                                     "../src/edmTail.cpp",
                                     "../src/edmx.cpp",
                                     "../src/helper.cpp"],
                          swig_opts=['-c++',
                                     '-I../src'],
                          include_dirs=['../src']
                          )

setup(name = "breakout_detection",
      version = '0.1',
      description = 'Breakout Detection',
      ext_modules = [extension_mod],
      py_modules = ["breakout_detection"])
