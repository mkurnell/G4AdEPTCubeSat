#ifndef PhysicsList_h
#define PhysicsList_h 1

#include "G4VModularPhysicsList.hh"
#include "globals.hh"

class G4VPhysicsConstructor;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class PhysicsList: public G4VModularPhysicsList
{
public:

	PhysicsList();
  	virtual ~PhysicsList();

  	virtual void ConstructParticle();
  	
  	virtual void ConstructProcess();
        
  	void AddPhysicsList(const G4String& name);
  	
  	virtual void SetCuts();

private:

	void AddPAIModel(const G4String&);
	void NewPAIModel(const G4ParticleDefinition* part, 
                     const G4String& modname,
                     const G4String& procname);

 	G4VPhysicsConstructor*  fEmPhysicsList;
  	G4VPhysicsConstructor*  fDecayPhysicsList;
  	std::vector<G4VPhysicsConstructor*> fHadronPhys;
  	G4String fEmName;

};

#endif