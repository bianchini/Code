Code
====

Instructions

0. Download the code to `$CMSSW_BASE/src/TTH/MEIntegratorStandalone`
1. copy the libraries from `lib/*.so` to `$CMSSW_BASE/lib/$SCRAM_ARCH/`
2. Make sure you have configured lhapdf: `scram setup lhapdf`
3. compile the MEM code with `scram b`
4. Check the example of how to call the MEM in `bin/integrator.cpp`
