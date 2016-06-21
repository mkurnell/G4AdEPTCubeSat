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
    ID_PVSensitiveGas_eDep_Positron = SDMan->GetCollectionID("PVSensitiveGas/eDepP");
	ID_PVSensitiveGas_eDep_Electron = SDMan->GetCollectionID("PVSensitiveGas/eDepE");
	ID_PVSensitiveGas_eDep_Triton = SDMan->GetCollectionID("PVSensitiveGas/eDepT");
	ID_PVSensitiveGas_trackLengthPassage = SDMan->GetCollectionID("PVSensitiveGas/trackLengthPassage");
	ID_PVSensitiveGas_secondaryElectrons = SDMan->GetCollectionID("PVSensitiveGas/secondaryElectrons");
	ID_PVSensitiveGas_secondaryPhotons = SDMan->GetCollectionID("PVSensitiveGas/secondaryPhotons");
	ID_PVSensitiveGas_secondaryPositrons = SDMan->GetCollectionID("PVSensitiveGas/secondaryPositrons");
	ID_PVSensitiveGas_secondaryTritons = SDMan->GetCollectionID("PVSensitiveGas/secondaryTritons");	
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
	G4double PVSensitiveGas_eDep_Positron = 0.;
	G4double PVSensitiveGas_eDep_Electron = 0.;
	G4double PVSensitiveGas_eDep_Triton = 0.;
	G4double PVSensitiveGas_trackLengthPassage = 0.;
	G4double PVSensitiveGas_secondaryElectrons = 0.;
	G4double PVSensitiveGas_secondaryPhotons = 0.;
	G4double PVSensitiveGas_secondaryPositrons = 0.;
	G4double PVSensitiveGas_secondaryTritons = 0.;
	
	// Get the HitMaps for this event
	G4THitsMap<G4double>* event_PVSensitiveGas_eDep = (G4THitsMap<G4double>*)(HCE->GetHC(ID_PVSensitiveGas_eDep));
	G4THitsMap<G4double>* event_PVSensitiveGas_eDep_Positron = (G4THitsMap<G4double>*)(HCE->GetHC(ID_PVSensitiveGas_eDep_Positron));
	G4THitsMap<G4double>* event_PVSensitiveGas_eDep_Electron = (G4THitsMap<G4double>*)(HCE->GetHC(ID_PVSensitiveGas_eDep_Electron));
	G4THitsMap<G4double>* event_PVSensitiveGas_eDep_Triton = (G4THitsMap<G4double>*)(HCE->GetHC(ID_PVSensitiveGas_eDep_Triton));
	G4THitsMap<G4double>* event_PVSensitiveGas_trackLengthPassage = (G4THitsMap<G4double>*)(HCE->GetHC(ID_PVSensitiveGas_trackLengthPassage));
	G4THitsMap<G4double>* event_PVSensitiveGas_secondaryElectrons = (G4THitsMap<G4double>*)(HCE->GetHC(ID_PVSensitiveGas_secondaryElectrons));
	G4THitsMap<G4double>* event_PVSensitiveGas_secondaryPhotons = (G4THitsMap<G4double>*)(HCE->GetHC(ID_PVSensitiveGas_secondaryPhotons));
	G4THitsMap<G4double>* event_PVSensitiveGas_secondaryPositrons = (G4THitsMap<G4double>*)(HCE->GetHC(ID_PVSensitiveGas_secondaryPositrons));
	G4THitsMap<G4double>* event_PVSensitiveGas_secondaryTritons = (G4THitsMap<G4double>*)(HCE->GetHC(ID_PVSensitiveGas_secondaryTritons));
	
	std::map<G4int,G4double*>::iterator itr;
	
	// Get the total energy deposited for this event in the Sensitive Gas Volume
	for (itr = event_PVSensitiveGas_eDep->GetMap()->begin(); itr != event_PVSensitiveGas_eDep->GetMap()->end(); itr++) {
		PVSensitiveGas_eDep += *(itr->second);
	}
	
	// Get the positron energy deposited for this event in the Sensitive Gas Volume
	for (itr = event_PVSensitiveGas_eDep_Positron->GetMap()->begin(); itr != event_PVSensitiveGas_eDep_Positron->GetMap()->end(); itr++) {
		PVSensitiveGas_eDep_Positron += *(itr->second);
	}
	
	// Get the electron energy deposited for this event in the Sensitive Gas Volume
	for (itr = event_PVSensitiveGas_eDep_Electron->GetMap()->begin(); itr != event_PVSensitiveGas_eDep_Electron->GetMap()->end(); itr++) {
		PVSensitiveGas_eDep_Electron += *(itr->second);
	}
	
	// Get the triton energy deposited for this event in the Sensitive Gas Volume
	for (itr = event_PVSensitiveGas_eDep_Triton->GetMap()->begin(); itr != event_PVSensitiveGas_eDep_Triton->GetMap()->end(); itr++) {
		PVSensitiveGas_eDep_Triton += *(itr->second);
	}
	
	// Get the passage track length for this event in the Sensitive Gas Volume
	for (itr = event_PVSensitiveGas_trackLengthPassage->GetMap()->begin(); itr != event_PVSensitiveGas_trackLengthPassage->GetMap()->end(); itr++) {
		PVSensitiveGas_trackLengthPassage += *(itr->second);
	} 

	// Get the number of secondary electrons in the Sensitive Gas Volume
	for (itr = event_PVSensitiveGas_secondaryElectrons->GetMap()->begin(); itr != event_PVSensitiveGas_secondaryElectrons->GetMap()->end(); itr++) {
		PVSensitiveGas_secondaryElectrons += *(itr->second);
	}
	
	// Get the number of secondary photons in the Sensitive Gas Volume
	for (itr = event_PVSensitiveGas_secondaryPhotons->GetMap()->begin(); itr != event_PVSensitiveGas_secondaryPhotons->GetMap()->end(); itr++) {
		PVSensitiveGas_secondaryPhotons += *(itr->second);
	}
	
	// Get the number of secondary positrons in the Sensitive Gas Volume
	for (itr = event_PVSensitiveGas_secondaryPositrons->GetMap()->begin(); itr != event_PVSensitiveGas_secondaryPositrons->GetMap()->end(); itr++) {
		PVSensitiveGas_secondaryPositrons += *(itr->second);
	}
	
	// Get the number of secondary tritons in the Sensitive Gas Volume
	for (itr = event_PVSensitiveGas_secondaryTritons->GetMap()->begin(); itr != event_PVSensitiveGas_secondaryTritons->GetMap()->end(); itr++) {
		PVSensitiveGas_secondaryTritons += *(itr->second);
	}
	
	// Record Sensitive Gas events with non-zero deposited energy
	if (PVSensitiveGas_eDep > 0) {
		// Get analysis manager
  		G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
  		
  	// Fill ntuple
	  	analysisManager->FillNtupleDColumn(0, PVSensitiveGas_eDep/eV);
  		analysisManager->FillNtupleDColumn(1, PVSensitiveGas_eDep_Positron/eV);
		analysisManager->FillNtupleDColumn(2, PVSensitiveGas_eDep_Electron/eV);
		analysisManager->FillNtupleDColumn(3, PVSensitiveGas_eDep_Triton/eV);
  		analysisManager->FillNtupleDColumn(4, PVSensitiveGas_trackLengthPassage/mm);
  		analysisManager->FillNtupleDColumn(5, PVSensitiveGas_secondaryElectrons);
		analysisManager->FillNtupleDColumn(6, PVSensitiveGas_secondaryPhotons);
		analysisManager->FillNtupleDColumn(7, PVSensitiveGas_secondaryPositrons);
		analysisManager->FillNtupleDColumn(8, PVSensitiveGas_secondaryTritons);
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