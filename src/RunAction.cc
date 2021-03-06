#include "RunAction.hh"
#include "Run.hh"
#include "PrimaryGeneratorAction.hh"
#include "DetectorConstruction.hh"
#include "G4Run.hh"
#include "G4RunManager.hh"
#include "G4UImanager.hh"
#include "G4VVisManager.hh"
#include "G4SystemOfUnits.hh"

// Select output format for Analysis Manager
#include "Analysis.hh"

#include <stdio.h>
#include <time.h>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

RunAction::RunAction(DetectorConstruction* det, PrimaryGeneratorAction* primary):G4UserRunAction(),
detector(det), particleGun(primary)
{
	// Set printing event number per each event
  	G4RunManager::GetRunManager()->SetPrintProgress(1E5);  
  	
	// Set starting seed for the Random Number Generator
	long seeds[2];
	time_t systime = time(NULL);
	seeds[0] = (long) systime;
	seeds[1] = (long) (systime*G4UniformRand());  
    G4Random::setTheSeeds(seeds);
  	
  	// Create analysis manager 
  	// The choice of analysis technology is done via selection of an appropriate namespace
  	// (g4root, g4xml, or g4csv)
  	G4AnalysisManager* analysisManager = G4AnalysisManager::Instance(); 
  	G4cout << "Using " << analysisManager->GetType() << G4endl; 
  	
	// Default settings 
  	analysisManager->SetVerboseLevel(0); 
  	
  	// Create ntuple 
  	// Pressure Sensitive Gas Volume
 	analysisManager->CreateNtuple("G4AdEPTCubeSat", "Edep and TrackLength");
	analysisManager->CreateNtupleDColumn("eDep_PVSensitiveGas");
 	analysisManager->CreateNtupleDColumn("eDep_PVSensitiveGas_Positron");
	analysisManager->CreateNtupleDColumn("eDep_PVSensitiveGas_Electron");
	analysisManager->CreateNtupleDColumn("eDep_PVSensitiveGas_Triton");
	analysisManager->CreateNtupleDColumn("eDep_PVSensitiveGas_Proton");
 	analysisManager->CreateNtupleDColumn("trackLength_PVSensitiveGas");
	analysisManager->CreateNtupleDColumn("Secondary Electrons");
	analysisManager->CreateNtupleDColumn("Secondary Photons");
	analysisManager->CreateNtupleDColumn("Secondary Positrons");
	analysisManager->CreateNtupleDColumn("Secondary Tritons");
	analysisManager->CreateNtupleDColumn("Secondary Protons");
 	analysisManager->FinishNtuple();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

RunAction::~RunAction()
{
	// Delete analysis manager
	delete G4AnalysisManager::Instance();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4Run* RunAction::GenerateRun()
{ 
	return new Run; 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::BeginOfRunAction(const G4Run* aRun)
{
	// Get analysis manager
  	G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();

  	// Open an AnalysisManager data file for the worker threads
  	if (!IsMaster()){
		// Filename for AnalysisManager is provided in the macro file using
		// /analysis/setFileName command
		analysisManager->OpenFile();
  	}
  	
  	// For the master let's create an info file
  	if (IsMaster()){
		// Get the local time at the start of the simulation
		time_t now = time(0);
		
		// Create an information file for the run using the same filename as the Analysis Manager
    	outputFile_INFO = analysisManager->GetFileName() + ".info";
    	
   		// Open the Information File
		pFile_INFO = fopen(outputFile_INFO,"w+");
		std::ofstream outFile_INFO(outputFile_INFO);
		
		// Export Source Information
		outFile_INFO << "============================    Simulation Information    ============================" << G4endl;
		outFile_INFO << "Start Time: \t\t" <<  ctime(&now);
  	}
  	
}

// //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::EndOfRunAction(const G4Run* aRun)
{ 	
  	// Output & close analysis file 
	G4AnalysisManager* analysisManager = G4AnalysisManager::Instance(); 
	if (!IsMaster()){
  		analysisManager->CloseFile(); 
  	}
  	
  	// Append Source Information to the INFO file
  	if (IsMaster()){
		// Open the Information File
		std::ofstream outFile_INFO(outputFile_INFO,std::ios::out|std::ios::app);
		
		//Get the local time at the end of the simulation
		time_t now = time(0);
    	
    	//Export Source Information
    	outFile_INFO << "End Time: \t\t\t" <<  ctime(&now);
		outFile_INFO << "============================    Source Information    ============================" << G4endl;
		outFile_INFO <<  "Number of Events: \t" << aRun->GetNumberOfEvent() << G4endl;	
		//outFile_INFO << "============================    Detector Information    ============================" << G4endl;
		//outFile_INFO <<  "Number of Ionizations: \t" << detector->GetDetectorAngle()/degree << " deg" << G4endl;	
		outFile_INFO << "==================================================================================" << G4endl; 
		
		//Close file
		fclose(pFile_INFO);
  	}
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
