source ~/geant4/geant4.10.2-install/bin/geant4.sh

cmake -DGeant4_DIR=~/geant4/geant4.10.2-install/lib/Geant4-10.2.1 ~/git/G4AdE
PTCubeSat/G4AdEPTCubeSat-4U/; make clean; make -j8; ./AdEPTCubeSat runGammas_ISO.mac



