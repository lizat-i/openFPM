# openFPM
Backbone of the implementation is the openFPM Framework:
OpenFPM is a scalable and open C++ framework for particles and mesh simulation.
see also : http://openfpm.mpi-cbg.de/

Implementation of the multiphase pcisph according to 
https://diglib.tugraz.at/download.php?id=5c4a48e8ceb60&location=browse
# Files 
main.cpp        contains main loop
main.h          contais functions from the openFPM SPH example implementation
                -domain setup, kernel etc.
constants.h         contains all the paramters for the simulation
PCISPHFunctions.h   contains the Functions used to perfrom the PCISPH loop

# how to get it running, first time setup:
this repo contains a 
.vscode and .devcontainer file.
It runs everything inside a docker container.

#1 docker pull openfpm/ubuntu:install20.04
#2 open vscode from within PCISPH folder 
#3 open inside Dev Container 
#4 cd /openfpm/openfpm_pdata
#5 ./install && sudo make install 
#6 follow installation instruction


# Compilation and run simulation 

#1 cd /workspace
#2 source /root/openfpm_vars
#3 make clean   cleans ouput and binaries
#4 make -j2     compile the main.cpp
#5 make run     runs te binarie using one core
#6 output will be written into output folder as vtk file
#7 logfile.txt will capture simulation log
#8  in order to open vtk file in paraview the 
    Intdebugger.sh file needs to be executes to rename datatype in the vrk file
