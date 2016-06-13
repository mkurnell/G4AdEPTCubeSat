// ********************************************************************
// DetectorConstruction.cc
//
// Description: Definition of the AdEPT CubeSat detector construction.
//				Includes the detector's pressure vessel (PV), Ar+CS2 gas
//				volume, various PCBs and the AdEPT Micro Well Detector (MWD)
//
// Created: Mitchell D. Kurnell
// 
// ********************************************************************

#include "DetectorConstruction.hh"
#include <cmath>

// Units and constants
#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"

// Manager classes
#include "G4RunManager.hh"
#include "G4NistManager.hh"
#include "G4GeometryManager.hh"
#include "G4SDManager.hh"

// Store classes
#include "G4PhysicalVolumeStore.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4SolidStore.hh"

// Geometry classes
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4RotationMatrix.hh"
#include "G4Transform3D.hh"
#include "G4VisAttributes.hh"

// Primitive geometry types
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4Sphere.hh"

// Boolean operations on volumes
#include "G4UnionSolid.hh"
#include "G4SubtractionSolid.hh"

// Regions
#include "G4Region.hh"

// Messenger classes
#include "G4GenericMessenger.hh"

// Scoring Components
#include "G4MultiFunctionalDetector.hh"
#include "G4PSEnergyDeposit.hh"
#include "G4PSTrackLength.hh"
#include "G4PSPassageTrackLength.hh"
#include "G4PSPassageCellCurrent.hh"
#include "G4PSNofSecondary.hh"
#include "G4VSDFilter.hh"
#include "G4SDParticleFilter.hh"
#include "G4ProductionCuts.hh"


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::DetectorConstruction(): G4VUserDetectorConstruction(), fCheckOverlaps(true),
WorldPhysical(0)
{	
	// Geometry Parameters (Default)
	// Pressure Vessel Top
	PV_length = 100.*mm;
	PV_width = 100.*mm;
	PV_height = 109.*mm;
	
	// Pressure Vessel Top Side Cutouts
	PV_sidecut_length = 83.*mm;
	PV_sidecut_width = 83.*mm;
	PV_sidecut_depth = 5.5*mm;
	
	// Pressure Vessel Top Cap Cutouts
	PV_topcut_length = 83.*mm;
	PV_topcut_width = 83.*mm;
	PV_topcut_depth = 1.*mm;
	
	// Pressure Vessel Bottom
	PV_bottom_length = 83.*mm;
	PV_bottom_width = 83.*mm;
	PV_bottom_height = 12.1*mm;
	
	// Pressure Vessel Bottom Cutouts
	PV_bottom_cut_length = 66.*mm;
	PV_bottom_cut_width = 66.*mm;
	PV_bottom_cut_depth = 1.*mm;
	
	// Pressure Vessel Top Gas Volume
	PV_gas_length = 85.*mm;
	PV_gas_width = 85.*mm;
	PV_gas_height = 103.5*mm;
	
	// Pressure Vessel Intermediate Gas Volume
	PV_mid_gas_length = 73.*mm;
	PV_mid_gas_width = 73*mm;
	PV_mid_gas_height = 2.5*mm;
	
	// Pressure Vessel Bottom Gas Volume
	PV_bottom_gas_length = 73.*mm;
	PV_bottom_gas_width = 73.*mm;
	PV_bottom_gas_height = 11.6*mm;
	
	// Presure Vessel Sensitive Gas Volume
	PV_sensitive_gas_length = 78.5*mm;
	PV_sensitive_gas_width = 78.5*mm;
	PV_sensitive_gas_height = 95.*mm;
	
	// Inner Electric Field Cage
	Cage_length = 85.*mm;
	Cage_width = 85.*mm;
	Cage_height = 95.*mm;
	
	// Pressure Vessel Top PCBs
	Top_PCB_length = 85.*mm;
	Top_PCB_width = 85.*mm;
	Top_PCB_thickness = 1.6*mm;
	
	// Bottom PCB
	Bottom_PCB_length = 73.*mm;
	Bottom_PCB_width = 73.*mm;
	Bottom_PCB_thickness = 1.6*mm;
	
	// Micro-Well Detector
	MWD_length = 77.5*mm;
	MWD_width = 77.5*mm;
	MWD_thickness = 1.*mm;
			
	// Rotation Angle
	rotX = 0.0*deg;		
			 
	// Define Materials
	DefineMaterials();
	
	// Define commands to control the geometry
   	DefineCommands();
	   
		// Production cuts for secondary particle generation
	G4double cut = 50.*um;
  	fTrackerCuts = new G4ProductionCuts();
  	fTrackerCuts->SetProductionCut(cut,"gamma");
  	fTrackerCuts->SetProductionCut(cut,"e-");
  	fTrackerCuts->SetProductionCut(cut,"e+");
  	fTrackerCuts->SetProductionCut(cut,"proton");
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::~DetectorConstruction()
{
	delete fTrackerCuts; 
}



//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

// DetectorConstruction::~DetectorConstruction()
// {}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::DefineMaterials()
{
	// NIST Manager
	G4NistManager* nistManager = G4NistManager::Instance();
	nistManager->SetVerbose(0);
	
	// NIST elements
  	G4Element* elH = nistManager->FindOrBuildElement(1);
  	G4Element* elC = nistManager->FindOrBuildElement(6);
  	G4Element* elO = nistManager->FindOrBuildElement(8);
  	G4Element* elSi = nistManager->FindOrBuildElement(14);

  	// NIST materials
  	G4Material* galactic = nistManager->FindOrBuildMaterial("G4_Galactic");
  	G4Material* Al = nistManager->FindOrBuildMaterial("G4_Al");
	G4Material* Si = nistManager->FindOrBuildMaterial("G4_Si");
  	
	// G10 material for the PCB
	G4double density = 1.70*g/cm3;
	G4Material* G10 = new G4Material("G10", density, 4);
	G10->AddElement(elH, 3);
	G10->AddElement(elC, 3);
	G10->AddElement(elO, 2);
	G10->AddElement(elSi, 1);
	
	// Detector Gas
	G4Material* Ar_293K_1p5atm = nistManager->ConstructNewGasMaterial("Ar_293K_1p5atm","G4_Ar",293.15*kelvin,1.5*atmosphere);
	
  	// Set the materials for the Geometry
  	fMatWorld = galactic;
  	fMatPressureVessel = Al;
  	fMatGas = Ar_293K_1p5atm;
	fMatPCB = G10;
	fMatMWD = Si; 
  	
  	// Print materials
// 	G4cout << *(G4Material::GetMaterialTable()) << G4endl;

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* DetectorConstruction::Construct()
{ 	
	// Cleanup
	// G4GeometryManager::GetInstance()->OpenGeometry();	
  	// if(fRegGasDet) { delete fRegGasDet; }
  	// fRegGasDet = new G4Region("Region_Sensitive_Gas");
  	// fRegGasDet->SetProductionCuts(fTrackerCuts);
	
	// Cleanup old geometry
  	G4GeometryManager::GetInstance()->OpenGeometry();
  	G4PhysicalVolumeStore::GetInstance()->Clean();
  	G4LogicalVolumeStore::GetInstance()->Clean();
  	G4SolidStore::GetInstance()->Clean();	
  	
		
	////////////////////////////////////////////////////////////////////////
	// Construct The World Volume (Vacuum)

	G4double world_X = 2*PV_gas_length;
	G4double world_Y = 2*PV_gas_width;
	G4double world_Z = 2*PV_gas_height;
	
	G4Box* WorldSolid = new G4Box("World", world_X/2, world_Y/2, world_Z/2);
  
	WorldLogical = 
		new G4LogicalVolume(WorldSolid,						// The Solid
							fMatWorld,					// Material
							"World");						// Name
  
	WorldPhysical = 
		new G4PVPlacement(	0,								// Rotation
							G4ThreeVector(),				// Translation vector
							WorldLogical,					// Logical volume
							"World",						// Name
							0,								// Mother volume
							false,							// Unused boolean parameter
							0,								// Copy number
							fCheckOverlaps);				// Overlap Check
							
	////////////////////////////////////////////////////////////////////////
	// Presure Vessel Construct
	
	G4Box* PVtopSolid = new G4Box("PressureVesselTop", PV_length/2, PV_width/2, PV_height/2);
	G4Box* PVbottomSolid = new G4Box("PressureVesselBottom", PV_bottom_length/2, PV_bottom_width/2, PV_bottom_height/2);
	
	G4UnionSolid* PVSolidBody = new G4UnionSolid("PressureVessel", PVtopSolid, PVbottomSolid, 0, G4ThreeVector(0,0,(PV_height+PV_bottom_height)/2 - 0.001*mm));

	
	////////////////////////////////////////////////////////////////////////
	// Pressure Vessel Side Cutouts
	
	G4Box* PVSideCut = new G4Box("PressureVesselSideCut", PV_sidecut_depth/2, PV_sidecut_length/2, PV_sidecut_width/2);
	
	// PVSideCutLogical_1 = 
	// 	new G4LogicalVolume(PVSideCut_1,
	// 						fMatWorld,
	// 						"PressureVesselSideCut");
	
	// PVSideCutPhysical_1 = 
	// 	new G4PVPlacement(	0,
	// 						G4ThreeVector((PV_width-PV_sidecut_depth)/2,0,0),
	// 						PVSideCutLogical_1,
	// 						"PressureVesselSideCut",
	// 						WorldLogical,
	// 						false,
	// 						0,
	// 						fCheckOverlaps);
	
	// //////////						
							
	// G4Box* PVSideCut_2 = new G4Box("PressureVesselSideCut", PV_sidecut_depth/2, PV_sidecut_length/2, PV_sidecut_width/2);
	
	// PVSideCutLogical_2 = 
	// 	new G4LogicalVolume(PVSideCut_2,
	// 						fMatWorld,
	// 						"PressureVesselSideCut");
	
	// PVSideCutPhysical_2 = 
	// 	new G4PVPlacement(	0,
	// 						G4ThreeVector(-(PV_width-PV_sidecut_depth)/2,0,0),
	// 						PVSideCutLogical_2,
	// 						"PressureVesselSideCut",
	// 						WorldLogical,
	// 						false,
	// 						0,
	// 						fCheckOverlaps);
							
	// //////////		
							
	G4Box* PVSideCut_3 = new G4Box("PressureVesselSideCut", PV_sidecut_length/2, PV_sidecut_depth/2, PV_sidecut_width/2);
	
	// PVSideCutLogical_3 = 
	// 	new G4LogicalVolume(PVSideCut_3,
	// 						fMatWorld,
	// 						"PressureVesselSideCut");
	
	// PVSideCutPhysical_3 = 
	// 	new G4PVPlacement(	0,
	// 						G4ThreeVector(0,(PV_width-PV_sidecut_depth)/2,0),
	// 						PVSideCutLogical_3,
	// 						"PressureVesselSideCut",
	// 						WorldLogical,
	// 						false,
	// 						0,
	// 						fCheckOverlaps);
							
	// //////////		

	// G4Box* PVSideCut_4 = new G4Box("PressureVesselSideCut", PV_sidecut_length/2, PV_sidecut_depth/2, PV_sidecut_width/2);
	
	// PVSideCutLogical_4 = 
	// 	new G4LogicalVolume(PVSideCut_4,
	// 						fMatWorld,
	// 						"PressureVesselSideCut");
	
	// PVSideCutPhysical_4 = 
	// 	new G4PVPlacement(	0,
	// 						G4ThreeVector(0,-(PV_width-PV_sidecut_depth)/2,0),
	// 						PVSideCutLogical_4,
	// 						"PressureVesselSideCut",
	// 						WorldLogical,
	// 						false,
	// 						0,
	// 						fCheckOverlaps);
	
	// ////////////////////////////////////////////////////////////////////////
	// // Pressure Vessel Top Cutouts
	
	G4Box* PVTopCut = new G4Box("PressureVesselTopCut", PV_topcut_length/2, PV_topcut_width/2, PV_topcut_depth/2);
	
	// PVTopCutLogical = 
	// 	new G4LogicalVolume(PVTopCut,
	// 						fMatWorld,
	// 						"PressureVesselTopCut");
							
	// PVTopCutPhysical = 
	// 	new G4PVPlacement(	0,
	// 						G4ThreeVector(0,0,-(PV_height-PV_topcut_depth)/2),
	// 						PVTopCutLogical,
	// 						"PressureVesselTopCut",
	// 						WorldLogical,
	// 						false,
	// 						0,
	// 						fCheckOverlaps);
								
	////////////////////////////////////////////////////////////////////////
	// Pressure Vessel Bottom Cutouts

	// G4Box* PVBottomCut = new G4Box("PressureVesselBottomCut", PV_bottom_cut_length/2, PV_bottom_cut_width/2, PV_bottom_cut_depth/2);
	
	// PVBottomCutLogical = 
	// 	new G4LogicalVolume(PVBottomCut,
	// 						fMatWorld,
	// 						"PressureVesselBottomCut");
							
	// PVBottomCutPhysical =
	// 	new G4PVPlacement(	0,
	// 						G4ThreeVector(0,0,(PV_height+PV_bottom_height*2-PV_bottom_cut_depth*2)/2),
	// 						PVBottomCutLogical,
	// 						"PressureVesselBottomCut",
	// 						WorldLogical,
	// 						false,
	// 						0,
	// 						fCheckOverlaps);
	
	////////////////////////////////////////////////////////////////////////
	// Subtraction of cutouts from Pressure Vessel
	
	G4SubtractionSolid* PVSolid_1 = 
		new G4SubtractionSolid("PressureVesselSolid", PVSolidBody, PVSideCut, 0, G4ThreeVector((PV_width-PV_sidecut_depth)/2,0,0));
	
	G4SubtractionSolid* PVSolid_2 = 
		new G4SubtractionSolid("PressureVesselSolid", PVSolid_1, PVSideCut, 0, G4ThreeVector(-(PV_width-PV_sidecut_depth)/2,0,0));
	
	G4SubtractionSolid* PVSolid_3 = 
		new G4SubtractionSolid("PressureVesselSolid", PVSolid_2, PVSideCut_3, 0, G4ThreeVector(0,(PV_width-PV_sidecut_depth)/2,0));
	
	G4SubtractionSolid* PVSolid_4 = 
		new G4SubtractionSolid("PressureVesselSolid", PVSolid_3, PVSideCut_3, 0, G4ThreeVector(0,-(PV_width-PV_sidecut_depth)/2,0));
	
	G4SubtractionSolid* PVSolid =
		new G4SubtractionSolid("PressureVesselSolid", PVSolid_4, PVTopCut, 0, G4ThreeVector(0,0,-(PV_height-PV_topcut_depth)/2));
	
	PVLogical =
		new G4LogicalVolume(PVSolid,
							fMatPressureVessel,
							"PressureVessel");
							
	PVPhysical = 
		new G4PVPlacement(	0,
							G4ThreeVector(),
							PVLogical,
							"PressureVessel",
							WorldLogical,
							false,
							0,
							fCheckOverlaps);

	////////////////////////////////////////////////////////////////////////
	// Pressure Vessel Detector Gas
	
	G4Box* PVGasTop = new G4Box("PressureVesselGasTop", PV_gas_length/2, PV_gas_width/2, PV_gas_height/2);
	G4Box* PVGasMid = new G4Box("PressureVesselGasMid", PV_mid_gas_length/2, PV_mid_gas_width/2, PV_mid_gas_height/2);
	G4Box* PVGasBottom = new G4Box("PressureVesselGasBottom", PV_bottom_gas_length/2, PV_bottom_gas_width/2, PV_bottom_gas_height/2);
	
	G4UnionSolid* PVGasTemp = new G4UnionSolid("PressureVesselGasTemp", PVGasTop, PVGasMid, 0, G4ThreeVector(0,0,(PV_gas_height+PV_mid_gas_height)/2));
	G4UnionSolid* PVGas = new G4UnionSolid("PressureVesselGas", PVGasTemp, PVGasBottom, 0, G4ThreeVector(0,0,(PV_mid_gas_height+PV_gas_height+PV_bottom_gas_height)/2));
	
	PVGasLogical = 
		new G4LogicalVolume(PVGas,
							fMatGas,
							"PressureVesselGas");
							
	PVGasPhysical = 
		new G4PVPlacement(	0,
							G4ThreeVector(),
							PVGasLogical,
							"PressureVesselGas",
							PVLogical,
							false,
							0,
							fCheckOverlaps);

	G4Region* regPVGas = new G4Region("Region_PV_Gas");
  	PVGasLogical->SetRegion(regPVGas);
  	regPVGas->AddRootLogicalVolume(PVGasLogical);
	  
	////////////////////////////////////////////////////////////////////////
	// Pressure Vessel Sensitive Detector Gas					
	
	G4Box* PVGasSensitive = new G4Box("SensitiveGas", PV_sensitive_gas_length/2, PV_sensitive_gas_width/2, PV_sensitive_gas_height/2);
	
	PVSensitiveGasLogical =
		new G4LogicalVolume(PVGasSensitive,
							fMatGas,
							"SensitiveGas");
							
	PVSensitiveGasPhysical =
		new G4PVPlacement(	0,
							G4ThreeVector(0,0,-0.25),
							PVSensitiveGasLogical,
							"SensitiveGas",
							PVGasLogical,
							false,
							0,
							fCheckOverlaps);
							
	G4Region* regSensitiveGas = new G4Region("Region_Sensitive_Gas");
  	PVSensitiveGasLogical->SetRegion(regSensitiveGas);
  	regSensitiveGas->AddRootLogicalVolume(PVSensitiveGasLogical);
							
	////////////////////////////////////////////////////////////////////////
	// Top PCB in Top Pressure Vessel
	
	G4Box* TopPCB = new G4Box("TopPCB", Top_PCB_length/2, Top_PCB_width/2, Top_PCB_thickness/2);
	
	TopPCBLogical =
		new G4LogicalVolume(TopPCB,
							fMatPCB,
							"TopPCB");
	
	TopPCBPhysical =
		new G4PVPlacement(	0,
							G4ThreeVector(0,0,-47.75-0.8),//(-PV_gas_height+8.*mm+Top_PCB_thickness)/2),
							TopPCBLogical,
							"TopPCB",
							PVGasLogical,
							false,
							0,
							fCheckOverlaps);
							
	////////////////////////////////////////////////////////////////////////
	// Bottom PCB in Top Pressure Vessel
	
	G4Box* BottomPCB = new G4Box("BottomPCB", Top_PCB_length/2, Top_PCB_width/2, Top_PCB_thickness/2);
	
	BottomPCBLogical =
		new G4LogicalVolume(BottomPCB,
							fMatPCB,
							"BottomPCB");
	
	BottomPCBPhysical =
		new G4PVPlacement(	0,
							G4ThreeVector(0,0,47.25+0.8),//(PV_gas_height-Top_PCB_thickness-9*mm)/2),
							BottomPCBLogical,
							"BottomPCB",
							PVGasLogical,
							false,
							0,
							fCheckOverlaps);
							
	////////////////////////////////////////////////////////////////////////
	// PCB in Bottom Pressure Vessel
	
	G4Box* CapPCB = new G4Box("CapPCB", Bottom_PCB_length/2, Bottom_PCB_width/2, Bottom_PCB_thickness/2);
	
	CapPCBLogical =
		new G4LogicalVolume(CapPCB,
							fMatPCB,
							"CapPCB");
	
	CapPCBPhysical = 
		new G4PVPlacement(	0,
							G4ThreeVector(0,0,((PV_gas_height-Bottom_PCB_thickness)/2+PV_bottom_gas_height+PV_mid_gas_height/2)),
							CapPCBLogical,
							"CapPCB",
							PVGasLogical,
							false,
							0,
							fCheckOverlaps);
	
	////////////////////////////////////////////////////////////////////////
	// Micro-Well Detector
	
	G4Box* MWD = new G4Box("MicrowellDetector", MWD_length/2, MWD_width/2, MWD_thickness/2);
	
	MWDLogical =
		new G4LogicalVolume(MWD,
							fMatMWD,
							"MicrowellDetector");
							
	MWDPhysical =
		new G4PVPlacement(	0,
							G4ThreeVector(0,0,(-0.8+0.5)),//(PV_gas_height/2-5.6*mm-MWD_thickness)),
							MWDLogical,
							"MicrowellDetector",
							BottomPCBLogical,
							false,
							0,
							fCheckOverlaps);
	
	////////////////////////////////////////////////////////////////////////
	// Field Shaping Cage
	
	G4Box* CageSolid = new G4Box("FieldCageSolid", Cage_length/2, Cage_width/2, Cage_height/2);
	G4Box* CageCut = new G4Box("FieldCageCut", (Cage_length-1.*mm)/2, (Cage_width-1.*mm)/2, Cage_height/2);
	
	G4SubtractionSolid* Cage = 
		new G4SubtractionSolid("Cage", CageSolid, CageCut, 0, G4ThreeVector());
	
	CageLogical =
		new G4LogicalVolume(Cage,
							fMatPCB,
							"Cage");
							
	CagePhysical = 
		new G4PVPlacement(	0,
							G4ThreeVector(0,0,-0.25*mm),
							CageLogical,
							"Cage",
							PVGasLogical,
							false,
							0,
							fCheckOverlaps);
	
							
	////////////////////////////////////////////////////////////////////////
  	// Visualisation attributes
  	
  	// ** World Volume
  	G4VisAttributes* Vis_World = new G4VisAttributes(G4Colour(1.,1.,1.,0.));
  	Vis_World->SetForceWireframe(true);
  	WorldLogical->SetVisAttributes(Vis_World);
  	
  	// // ** Pressure Vessel Volume
  	G4VisAttributes* Vis_PV = new G4VisAttributes(G4Colour(1.,1.,1.,0.2));
  	Vis_PV->SetForceWireframe(false);
  	PVLogical->SetVisAttributes(Vis_PV);
	  
  	// ** Pressure Vessel Cutouts
  	// G4VisAttributes* Vis_Cutout = new G4VisAttributes(G4Colour(0.,0.,0.,0.2));
  	// Vis_Cutout->SetForceWireframe(false);
  	// PVSideCutLogical_1->SetVisAttributes(Vis_Cutout);
	// PVSideCutLogical_2->SetVisAttributes(Vis_Cutout);
	// PVSideCutLogical_3->SetVisAttributes(Vis_Cutout);
	// PVSideCutLogical_4->SetVisAttributes(Vis_Cutout);
  	// PVTopCutLogical->SetVisAttributes(Vis_Cutout);
	// PVBottomCutLogical->SetVisAttributes(Vis_Cutout);
  	
  	// ** Pressure Vessel Gas Volume
  	G4VisAttributes* Vis_Gas = new G4VisAttributes(G4Colour(0.,0.,1.,0.3));
  	Vis_Gas->SetForceWireframe(false);
  	PVGasLogical->SetVisAttributes(Vis_Gas);
	  
  	// ** Pressure Vessel Sensitive Gas Volume
  	G4VisAttributes* Vis_Sensitive_Gas = new G4VisAttributes(G4Colour(0.5,0.,0.5,0.32));
  	Vis_Sensitive_Gas->SetForceWireframe(false);
  	PVSensitiveGasLogical->SetVisAttributes(Vis_Sensitive_Gas);
  	
  	// ** PCBs
  	G4VisAttributes* Vis_PCB = new G4VisAttributes(G4Colour(0.,1.,0.,0.5));
  	Vis_PCB->SetForceWireframe(false);
  	TopPCBLogical->SetVisAttributes(Vis_PCB);
	BottomPCBLogical->SetVisAttributes(Vis_PCB);
	CapPCBLogical->SetVisAttributes(Vis_PCB);
	  
  	// ** Micro-Well Detector
  	G4VisAttributes* Vis_MWD = new G4VisAttributes(G4Colour(0.94,0.9,0.55,0.6));
  	Vis_MWD->SetForceWireframe(false);
  	MWDLogical->SetVisAttributes(Vis_MWD);
  	
  	// // ** Field Shaping Cage
  	G4VisAttributes* Vis_Cage = new G4VisAttributes(G4Colour(0.8,0.52,0.25,0.4));
  	Vis_Cage->SetForceWireframe(false);
  	CageLogical->SetVisAttributes(Vis_Cage);

	////////////////////////////////////////////////////////////////////////
	// Return world volume
	return WorldPhysical; 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::ConstructSDandField()
{ 	
  	////////////////////////////////////////////////////////////////////////
	// Construct the Multi Functional Detector for the Sensitive Gas Volume
	
	G4MultiFunctionalDetector* PVGasScorer = new G4MultiFunctionalDetector("PVSensitiveGas");
	G4SDManager::GetSDMpointer()->AddNewDetector(PVGasScorer);	
	G4SDManager::GetSDMpointer()->SetVerboseLevel(0);
	PVSensitiveGasLogical->SetSensitiveDetector(PVGasScorer);
	
	// MFD Fitlers
	G4VSDFilter* ElectronFilter = new G4SDParticleFilter("electronFilter","e-");
	G4VSDFilter* PositronFilter = new G4SDParticleFilter("positronFilter","e+");
	G4VSDFilter* PhotonFilter = new G4SDParticleFilter("gammaFilter","gamma");

	// Total Edep/, Count secondary photons/, count electrons/, positrons
	
	// Total Energy deposited in the Sensitive Gas Volume
	G4VPrimitiveScorer* eDep_Gas = new G4PSEnergyDeposit("eDep");
    PVGasScorer->RegisterPrimitive(eDep_Gas);
	
	// Energy deposited from positrons in the Sensitive Gas Volume;
	G4VPrimitiveScorer* eDep_Gas_Positron = new G4PSEnergyDeposit("eDepP");
	eDep_Gas_Positron->SetFilter(PositronFilter);
    PVGasScorer->RegisterPrimitive(eDep_Gas_Positron);
	
	// Energy deposited from Electrons in the Sensitive Gas Volume
	G4VPrimitiveScorer* eDep_Gas_Electron = new G4PSEnergyDeposit("eDepE");
	eDep_Gas_Electron->SetFilter(ElectronFilter);
	PVGasScorer->RegisterPrimitive(eDep_Gas_Electron);

    G4PSPassageTrackLength* trackLengthPassage_Gas = new G4PSPassageTrackLength("trackLengthPassage");
 	PVGasScorer->RegisterPrimitive(trackLengthPassage_Gas);
	 
	// Number of secondary electrons
	G4PSNofSecondary* secondaryElectrons = new G4PSNofSecondary("secondaryElectrons");
	secondaryElectrons->SetFilter(ElectronFilter);
	PVGasScorer->RegisterPrimitive(secondaryElectrons);
	
	// Number of secondary photons
	G4PSNofSecondary* secondaryPhotons = new G4PSNofSecondary("secondaryPhotons");
	secondaryPhotons->SetFilter(PhotonFilter);
	PVGasScorer->RegisterPrimitive(secondaryPhotons);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::DefineCommands()
{
    // Define /AdEPTCubeSat/ command directory using generic messenger class
    fMessenger = new G4GenericMessenger(this, "/AdEPTCubeSat/", "Geometry control");
}
