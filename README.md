# GEANT4_HPGe_Simulation

## Introduction
This program is based on GEANT4 simulation for HPGe detector
Initially Cs-137 source is defined by default.

## Directories


Directory | Contents
----------|-----------
src       | source code
include   | Header Files


## Get started

1. Make sure that [ROOT][] (version 6 and above) and GEANT4 (version 10 or above) is installed.


2. Execute the following commands in a terminal:
~~~sh
git clone https://github.com/SubhroD/GEANT4_HPGe_Simulation.git
cd GEANT4_HPGe_Simulation-main/
mkdir build
cd build/
cmake ../
make -jN
./exampleB4a
~~~

[ROOT]:https://root.cern.ch
[GEANT4]:https://geant4.web.cern.ch/
