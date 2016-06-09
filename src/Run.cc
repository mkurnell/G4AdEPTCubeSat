#include "Run.hh"
#include "G4Event.hh"
#include "G4Run.hh"
#include "G4SDManager.hh"
#include "G4HCofThisEvent.hh"
#include "G4THitsMap.hh"
#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"

// Select output format for Analysis Manager
#include "Analysis.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

Run::Run():G4Run()
{
	G4SDManager* SDMan = G4SDManager::GetSDMpointer(); 
    	
    ID_PVSensitiveGas_eDep = SDMan->GetCollectionID("PVSensitiveGas/eDep");
	ID_PVSensitiveGas_trackLengthPassage = SDMan->GetCollectionID("PVSensitiveGas/trackLengthPassage");
	ID_PVSensitiveGas_ionizations = SDMan->GetCollectionID("PVSensitiveGas/ionizations");
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

Run::~Run()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Run::RecordEvent(const G4Event* event)
{ 	
  	// Get hits collections
  	G4HCofThisEvent* HCE = event->GetHCofThisEvent();
  	if(!HCE) { 
    	G4ExceptionDescription msg; 
    	msg << "No hits collection of this event found.\n"; 
    	G4Exception("Run::RecordEvent()","Code001", JustWarning, msg); 
    	return; 
  	} 
  	
	// Zero out the variables
	G4double PVSensitiveGas_eDep = 0.;
	G4double PVSensitiveGas_trackLengthPassage = 0.;
	G4double PVSensitiveGas_ionizations = 0.;
	
	// Get the HitMaps for this event
	G4THitsMap<G4double>* event_PVSensitiveGas_eDep = (G4THitsMap<G4double>*)(HCE->GetHC(ID_PVSensitiveGas_eDep));
	G4THitsMap<G4double>* event_PVSensitiveGas_trackLengthPassage = (G4THitsMap<G4double>*)(HCE->GetHC(ID_PVSensitiveGas_trackLengthPassage));
	G4THitsMap<G4double>* event_PVSensitiveGas_ionizations = (G4THitsMap<G4double>*)(HCE->GetHC(ID_PVSensitiveGas_ionizations));
	
	std::map<G4int,G4double*>::iterator itr;
	
	// Get the total energy deposited for this event in the Sensitive Gas Volume
	for (itr = event_PVSensitiveGas_eDep->GetMap()->begin(); itr != event_PVSensitiveGas_eDep->GetMap()->end(); itr++) {
		PVSensitiveGas_eDep += *(itr->second);
	}
	
	// Get the passage track length for this event in the Sensitive Gas Volume
	for (itr = event_PVSensitiveGas_trackLengthPassage->GetMap()->begin(); itr != event_PVSensitiveGas_trackLengthPassage->GetMap()->end(); itr++) {
		PVSensitiveGas_trackLengthPassage += *(itr->second);
	} 

	// Get the number of secondary ionizations in the Sensitive Gas Volume
	for (itr = event_PVSensitiveGas_ionizations->GetMap()->begin(); itr != event_PVSensitiveGas_ionizations->GetMap()->end(); itr++) {
		PVSensitiveGas_ionizations += *(itr->second);
	}
	// Record Sensitive Gas events with non-zero deposited energy
	if (PVSensitiveGas_eDep > 0) {
		// Get analysis manager
  		G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
  		
  	// Fill ntuple
  		analysisManager->FillNtupleDColumn(0, PVSensitiveGas_eDep/eV);
  		analysisManager->FillNtupleDColumn(1, PVSensitiveGas_trackLengthPassage/mm);
  		analysisManager->FillNtupleDColumn(2, PVSensitiveGas_ionizations);
  		analysisManager->AddNtupleRow();
	}
	
	// Invoke base class method
  	G4Run::RecordEvent(event); 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
// 
// void Run::Merge(const G4Run* aRun)
// {
//   	const Run* localRun = static_cast<const Run*>(aRun);
//   	
//   	//  Invoke base class method
//   	G4Run::Merge(aRun); 
// } 