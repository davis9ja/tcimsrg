# Tensorized C++ IM-SRG

To build program, run

    mkdir build; cd build
    cmake .. ../build/

To execute program, run

   ./solve_imsrg <args>

## Executable arguments

 - numStates=<int> (specify number of S.P. states in reference)
 - nHoles=<int> (specify number of holes states in reference)
 - nParticles=<int> (specify number of particle states in reference)
 - d=<double> (energy spacing)
 - g=<double> (pairing strength)
 - pb=<double> (pair-breaking strength)
 - t0=<double> (starting point for ODE solver)
 - t1=<double> (ending point for ODE solver)
 - dt=<double> (step width for ODE solver)
 - generator_id=<0,1> (0: White generator with MP denominators, 1: Atan-regulated White)
 - reference_type=<int> (unimplemented)
 - solver=<0,1> (0:Adams-Bashforth-Moulton, 1:Runge-Kutta 4 point)

`CMakeLists.txt` will automatically link the main.cpp executable against the TACO library.

Tensor Algebra Compiler code at [https://github.com/tensor-compiler/taco]()

### TODO:

- [x] Implement TACO backend
  - pass occupation tensors directly so no conversion?
- [x] Add density matrix normal-ordering functionality
  - [ ] optimize density matrix normal-ordering
    - DMCPP code is way too slow
- [x] Remove particle-hole distinction
  - This makes density matrix process simpler
- [ ] Come up with user input to specify density matrix weights
- [x] Come up with user input to specify reference state type (single ref, or ensemble)
- [ ] Optimize loops and parallelize TACO code
- [x] Overhaul the system logging
  - want two files: IMSRG coefficients and vacuum coefficients (including s vals and ||eta1b||, ||eta2b||)
- [x] Add White Atan generator
  - seems to work better for floating point reference states, without including irreducible DMs in flow equations (should investigate this)
- [ ] Add the multi-reference IMSRG equations
  - the flow seems unstable for fractional density matrices without including irreducible components in flow equations

## Purpose
Build system for solving the pairing-plus-particle-hole model Hamiltonian using IMSRG(2). 

## Code Logic
Modularized so that new pieces can be added to code, including density matrix normal-ordering, new model Hamiltonians, and new backends for computing the flow equations. Logic is still WIP. 

`main.cpp`\
User-defined flow setup, handling, and solution. 

`pairinghamiltonian.cpp`\
Defines the pairing-plus-particle-hole model Hamiltonian and computes vacuum coefficients and normal-ordered coefficients (with reference state provided as input).

`white.cpp`\
Computes the White generator for input IMSRG coefficients. Inherits from `generator.hpp`.

`flow_imsrg2.cpp`\
Handles computation of the flow equations. Transforms occupation tensors to UBLAS vectors before passing to backend.

`system.cpp`\
ODE system class for `boost::numeric::odeint`. Track E, f, Gamma and dE, df, dGamma while passing through solver.

`system_observer.cpp`\
Singleton class that observes and outputs flow data. Passed to ODEINT solver.

`occupation_factors.cpp`\
Write to file tensor-train decomposed occupation factors that appear in the flow equations. Also handles reading back into memory pointer.

`BACKEND_ublas.cpp`\
Compute the IMSRG flow using only vector loops. Primarily for performance benchmarking. Inherits from `BACKEND.hpp`.

`BACKEND_taco.cpp`\
Compute the IMSRG flow using tensor library TACO. Convert UBLAS vectors into TACO tensors and run contractions. Inherits from `BACKEND.hpp`.
