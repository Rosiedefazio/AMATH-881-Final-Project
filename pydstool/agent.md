# Agent Guide

This repository is a legacy scientific Python package for modeling, simulating, and analyzing dynamical systems. Read this file before making changes so you do not trip over the native build paths, import side effects, or broad public API surface.

## Repository Map

- `PyDSTool/`: main package.
- `PyDSTool/__init__.py`: heavy package entry point that eagerly imports and re-exports much of the public API.
- `PyDSTool/Generator/`: solver wrappers and generator classes. `Dopri_ODEsystem.py` and `Radau_ODEsystem.py` use generated C code and native builds. `Vode_ODEsystem.py` and `Euler_ODEsystem.py` are important pure-Python/SciPy paths.
- `PyDSTool/core/codegenerators/`: code generation backends with focused tests.
- `PyDSTool/integrator/`: C, Fortran, SWIG, and Makefile assets used by compiled integrators.
- `PyDSTool/PyCont/`: continuation and bifurcation code. `PyDSTool/PyCont/auto/` contains vendored AUTO sources and should be treated as high risk.
- `PyDSTool/Toolbox/`: optional higher-level analysis utilities.
- `tests/`: pytest suite, roughly organized by package area.
- `examples/`: manual examples and smoke scripts; several require plotting or native toolchains.
- `PyCont_docs/`: historical PyCont documentation artifacts.

## Environment Expectations

- Package metadata targets Python 3.6+.
- `tox.ini` and `.travis.yml` are centered on Python 3.6 and 3.7-era environments.
- Required Python deps: `numpy`, `scipy`.
- Test deps: `pytest`, `pytest-mock`, `pytest-xdist`.
- Native build deps for compiled interfaces: `swig`, a C compiler, and `gfortran`.
- CI also installs LAPACK/Atlas development packages for compiled solver paths.

## Common Commands

- Install editable: `python setup.py develop`
- Run full test suite: `pytest`
- Alternate full test path: `python setup.py test`
- Run one test file: `pytest tests/test_common.py`
- Run codegen-focused tests: `pytest tests/core/codegenerators/test_c.py tests/generator/test_codegen.py`
- Run examples: `cd examples && python run_all_tests.py`
- Clean generated artifacts: `make clean`

`setup.cfg` defines pytest defaults: `-rsxX -q --boxed`, `testpaths = tests`, and warnings are treated as errors except for a short allowlist. Because of `--boxed`, `pytest-xdist` must be installed even for normal test runs.

## Working Rules

- Preserve public API compatibility unless the change explicitly intends to break it. `from PyDSTool import *` is a supported pattern here.
- Prefer small, local fixes over broad refactors. This codebase uses legacy module structure and older conventions.
- When changing imports or initialization code, always consider `import PyDSTool` behavior. Package import is broad and side-effectful.
- Keep pure-Python changes separate from native-build changes when possible. It makes breakage easier to isolate.
- Avoid editing vendored AUTO code under `PyDSTool/PyCont/auto/` unless the defect is clearly inside that subtree.
- Do not commit generated shared objects, temp build directories, or generated vector-field files.

## Native Build Notes

- Compiled generators create artifacts in the current working directory, not in a dedicated cache. Common outputs include `dop853_temp/`, `radau5_temp/`, `auto_temp/`, generated `*_vf.py`, `*_vf_wrap.c`, shared libraries, and build logs.
- `PyDSTool/Generator/mixins.py` documents a real limitation: once a compiled shared library exists, rebuilding cleanly in the same Python session is unreliable. If you change compiled solver code, delete generated artifacts and restart the Python process before retesting.
- For code generation work that does not need compilation, prefer existing `nobuild=True` paths used throughout `tests/generator/test_codegen.py` and related tests.

## Testing Strategy

- Start with the narrowest relevant tests.
- Good first checks for symbolic or code-generation changes:
  - `tests/core/codegenerators/`
  - `tests/generator/test_codegen.py`
  - `tests/test_funcspec.py`
- Native compilation coverage lives in places like:
  - `tests/generator/test_compiled_interfaces.py`
  - `tests/pycont/`
- The examples are useful smoke tests but are slower, may open matplotlib windows, and can require manual interaction depending on platform/toolchain.

## High-Risk Areas

- `PyDSTool/__init__.py`: changing exports or import order can break large parts of the package.
- `PyDSTool/Generator/` plus `PyDSTool/integrator/`: mixed Python/native boundary.
- `PyDSTool/PyCont/ContClass.py` and `PyDSTool/PyCont/auto/`: compiled continuation path with vendored code.
- Cleanup commands in `setup.py clean` and `Makefile`: they remove generated artifacts aggressively. Read them before extending cleanup behavior.

## Practical Advice

- If a test only needs generated source, avoid invoking the compiler.
- If a failure appears only after rebuilding a compiled module, suspect stale temp artifacts first.
- If you touch a public object or re-export, add at least one test that exercises the import path the way users do.
- When in doubt, prefer a targeted fix plus a focused pytest file over a full-repo rewrite.
