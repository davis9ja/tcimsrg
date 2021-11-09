# Tensorized C++ IM-SRG

To build program, run

    mkdir build; cd build
    cmake .. ../build/
    ./main

`CMakeLists.txt` will automatically link the main.cpp executable against the TACO library.

Tensor Algebra Compiler code at [https://github.com/tensor-compiler/taco]()

### TODO:

- [ ] Make SystemObserver a singleton class
  - necessary to reconcile flow data across all systems
- [ ] Implement TACO backend
  - pass occupation tensors directly so no conversion?
- [ ] Add density matrix normal-ordering functionality
- [ ] Remove particle-hole distinction
  - This makes density matrix process simpler

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
