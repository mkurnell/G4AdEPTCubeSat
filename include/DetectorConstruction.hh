#ifndef DetectorConstruction_h
#define DetectorConstruction_h 1

#include "G4VUserDetectorConstruction.hh"
#include "globals.hh"

class G4VPhysicalVolume;
class G4LogicalVolume;
class G4Material;
class G4GenericMessenger;
class G4ProductionCuts;

class DetectorConstruction : public G4VUserDetectorConstruction
{
  public:
  	// Constructor
    DetectorConstruction();
    // Destructor
    virtual ~DetectorConstruction();
    
  public:
    // Defines the detector geometry and returns a pointer to the physical World Volume
    virtual G4VPhysicalVolume* Construct();
    
    // Sensitive Detector 
	virtual void ConstructSDandField();
    
  public:
    // Set Methods
    void SetDetectorAngle(G4double val);
    
    // Get Methods
    G4double GetDetectorAngle();
    
  private:
    // Defines all the detector materials
    void DefineMaterials();
    
    // Define commands to change the geometry
    void DefineCommands();
    
    G4GenericMessenger* fMessenger;
    G4bool  fCheckOverlaps;
    
    // Standard Materials
    G4Material* fMatWorld;
    G4Material* fMatPressureVessel;
    G4Material* fMatGas;
    G4Material* fMatPCB;
    G4Material* fMatMWD;
    
    // Logical Volumes
    G4LogicalVolume* WorldLogical;
    G4LogicalVolume* PVLogical;
    G4LogicalVolume* PVGasLogical;
    G4LogicalVolume* PVSensitiveGasLogical;
    G4LogicalVolume* PVSideCutLogical_1;
    G4LogicalVolume* PVSideCutLogical_2;
    G4LogicalVolume* PVSideCutLogical_3;
    G4LogicalVolume* PVSideCutLogical_4;
    G4LogicalVolume* PVTopCutLogical;
    G4LogicalVolume* PVBottomCutLogical;
    G4LogicalVolume* TopPCBLogical;
    G4LogicalVolume* BottomPCBLogical;
    G4LogicalVolume* CapPCBLogical;
    G4LogicalVolume* MWDLogical;
    G4LogicalVolume* CageLogical;
    
    // Physical Volumes
    G4VPhysicalVolume* WorldPhysical;
    G4VPhysicalVolume* PVPhysical;
    G4VPhysicalVolume* PVGasPhysical;
    G4VPhysicalVolume* PVSensitiveGasPhysical;
    G4VPhysicalVolume* PVSideCutPhysical_1;
    G4VPhysicalVolume* PVSideCutPhysical_2;
    G4VPhysicalVolume* PVSideCutPhysical_3;
    G4VPhysicalVolume* PVSideCutPhysical_4;
    G4VPhysicalVolume* PVTopCutPhysical;
    G4VPhysicalVolume* PVBottomCutPhysical;
    G4VPhysicalVolume* TopPCBPhysical;
    G4VPhysicalVolume* BottomPCBPhysical;
    G4VPhysicalVolume* CapPCBPhysical;
    G4VPhysicalVolume* MWDPhysical;
    G4VPhysicalVolume* CagePhysical;
    
    // Geometry Parameters
    
    // Detector Dimensions
    G4double PV_length;
    G4double PV_width;
    G4double PV_height;
    
    G4double PV_sidecut_length;
    G4double PV_sidecut_width;
    G4double PV_sidecut_depth;
    
    G4double PV_topcut_length;
    G4double PV_topcut_width;
    G4double PV_topcut_depth;
    
    G4double PV_bottom_length;
    G4double PV_bottom_width;
    G4double PV_bottom_height;
    
    G4double PV_bottom_cut_length;
    G4double PV_bottom_cut_width;
    G4double PV_bottom_cut_depth;
    
    G4double PV_gas_length;
    G4double PV_gas_width;
    G4double PV_gas_height;
    
    G4double PV_mid_gas_length;
    G4double PV_mid_gas_width;
    G4double PV_mid_gas_height;
    
    G4double PV_bottom_gas_length;
    G4double PV_bottom_gas_width;
    G4double PV_bottom_gas_height;
    
    G4double PV_sensitive_gas_length;
    G4double PV_sensitive_gas_width;
    G4double PV_sensitive_gas_height;
    
    G4double Cage_length;
    G4double Cage_width;
    G4double Cage_height;
    
    G4double Top_PCB_length;
    G4double Top_PCB_width;
    G4double Top_PCB_thickness;
    
    G4double Bottom_PCB_length;
    G4double Bottom_PCB_width;
    G4double Bottom_PCB_thickness;
    
    G4double MWD_length;
    G4double MWD_width;
    G4double MWD_thickness;
    
	  // Rotation Angles
	  G4double rotX;
    G4ProductionCuts*  fTrackerCuts;
    
};

#endif