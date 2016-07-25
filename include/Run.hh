#ifndef Run_h
#define Run_h 1

#include "G4Run.hh"
#include "globals.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class Run : public G4Run
{
	public:
		// Constructor
  		Run();
  		// Destructor
  		virtual ~Run();
		
		// Methods
		virtual void RecordEvent(const G4Event*);
// 		virtual void Merge(const G4Run*);

	private:
		G4int ID_PVSensitiveGas_eDep;
		G4int ID_PVSensitiveGas_eDep_Positron;
		G4int ID_PVSensitiveGas_eDep_Electron;
		G4int ID_PVSensitiveGas_eDep_Triton;
		G4int ID_PVSensitiveGas_eDep_Proton;
		G4int ID_PVSensitiveGas_trackLengthPassage;
		G4int ID_PVSensitiveGas_secondaryElectrons;
		G4int ID_PVSensitiveGas_secondaryPositrons;
		G4int ID_PVSensitiveGas_secondaryPhotons;
		G4int ID_PVSensitiveGas_secondaryTritons;
		G4int ID_PVSensitiveGas_secondaryProtons;
};

#endif