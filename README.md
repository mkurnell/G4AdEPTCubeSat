# G4AdEPTCubeSat
G4AdEPTCubeSat is a simulation of the Advanced Energetic Pair Telescope's (AdEPT) imaging instrument scaled down for a cubesat mission. The application is based on the Geant4 Monte Carlo toolkit which allows simulation of the response of the instrument to different radiation environments.

# Build Notes
To build this application on your computer, ensure that you have a working version of CMake and Geant 4.10 or higher. You will also need OpenCV, OpenGL and Qt to be working with your Geant4 install if you wish to see a visualization of the detector.

For simplicity, I recommend using the following directory structure:

- $G4WORKDIR/G4AdEPTCubeSat : This directory will hold the information regarding the simulation and is a clone of this GitHub repository
- $G4WORKDIR/G4AdEPTCubeSat-build : This directory houses the build files for the simulation. This folder also contains the executable file and macros for running various radiation sources

## Steps to compile and run:
##
### Step 1: Source the Geant4 environment setup script
source /opt/local/libexec/Geant4/Geant4.10.1/geant4.sh

### Step 2: Create a local build directory and enter it
in $G4WORKDIR: mkdir G4AdEPTCubeSat-build && cd $_

### Step 3: Compile the code, make the build and run the executable
cmake -DGeant4_DIR=/opt/local/libexec/Geant4/Geant4.10.1
~/$G4WORKDIR/G4AdEPTCubeSat; make clean; make; ./AdEPTCubeSat

If you have setup Geant4 with multithreading, you may substitute "make" with "make -j8" to speed up this process.
After making changes to the code, simply repeat step three to run the code with the changes made.
You may also choose to run a macro file at this point (e.g. cmake -DGeant4_DIR=/opt/local/libexec/Geant4/Geant4.10.1
~/$G4WORKDIR/G4AdEPTCubeSat; make clean; make; ./AdEPTCubeSat runGammas_ISO.mac)

## Scoring of Events

This simulation will track the energy deposited in the Sensitive Gas Volume (SGV) by electrons, positrons and the sum of these two energies. In addition, it will also count the number of secondary electrons, positrons and photons created by these energy deposition events in the SGV. With the number of secondaries created tracked you will be able to determine how the photons interacted within the gas volume (Photoelectric Effect, Compton Scattering or Pair Production).

All of this information will be printed out in the build directory in the form of CSV files.
