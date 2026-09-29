"""Install/relocate the C-only package, then compile independent C/Fortran callers."""
import os
from pathlib import Path
import subprocess
import sys
import tempfile

cmake, build, source, fixtures, backend, gcc, gfortran = sys.argv[1:]
source = Path(source)


def run(*args, clean_loader=False):
    env = os.environ.copy()
    if clean_loader:
        # The relocated package must supply its own runtime search paths.
        for name in ('LD_LIBRARY_PATH', 'DYLD_LIBRARY_PATH', 'DYLD_FALLBACK_LIBRARY_PATH'):
            env.pop(name, None)
    subprocess.run(args, check=True, env=env)


with tempfile.TemporaryDirectory(prefix='accelnet-c-consumer-') as directory:
    root = Path(directory)
    run(cmake, '--install', build, '--prefix', str(root / 'install'), '--component', 'TargetC')
    prefix = root / 'relocated'
    (root / 'install').rename(prefix)
    project = root / 'project'
    project.mkdir()
    (project / 'CMakeLists.txt').write_text('''cmake_minimum_required(VERSION 3.20)
project(CConsumer LANGUAGES C)
find_package(AccelNetC CONFIG REQUIRED)
get_target_property(_options AccelNet::TargetC INTERFACE_LINK_OPTIONS)
if(_options)
  message(FATAL_ERROR "C ABI target leaked compiler link options: ${_options}")
endif()
add_executable(c_client "${CLIENT_SOURCE}/AccelNetPredictor/test/test_target_c_loading.c")
target_link_libraries(c_client PRIVATE AccelNet::TargetC)
if(CLIENT_FORTRAN)
  enable_language(Fortran)
  add_executable(fortran_client "${CLIENT_SOURCE}/docs/validation/target-c-loading-2026-09-28/gnu_client.f90")
  target_link_libraries(fortran_client PRIVATE AccelNet::TargetC)
endif()
''')
    out = root / 'build'
    run(cmake, '-S', str(project), '-B', str(out), '-DCMAKE_PREFIX_PATH=' + str(prefix),
        '-DCMAKE_C_COMPILER=' + gcc, '-DCLIENT_SOURCE=' + str(source))
    run(cmake, '--build', str(out), '--parallel', '2')
    run(str(out / 'c_client'), fixtures, str(source / 'AccelNetPredictor/test/data'), backend,
        clean_loader=True)
    # Only this second configuration enables Fortran; the C-only consumer above
    # must not need a Fortran compiler or native AccelNet module files at all.
    run(cmake, '-S', str(project), '-B', str(out), '-DCLIENT_FORTRAN=ON',
        '-DCMAKE_Fortran_COMPILER=' + gfortran)
    run(cmake, '--build', str(out), '--parallel', '2')
    run(str(out / 'fortran_client'), str(source / 'AccelNetPredictor/test/data/n2p2-virial-angular'),
        str(Path(fixtures) / 'n2p2-virial-angular.ref'), backend, clean_loader=True)
print('Relocated C-only package and independent GNU Fortran caller passed')
