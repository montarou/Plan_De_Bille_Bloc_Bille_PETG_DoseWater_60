#ifndef DetectorConstruction_h
#define DetectorConstruction_h 1

#include "G4VUserDetectorConstruction.hh"
#include "globals.hh"
#include <vector>

class G4VPhysicalVolume;
class G4LogicalVolume;
class G4Material;

class DetectorConstruction : public G4VUserDetectorConstruction
{
public:
    DetectorConstruction();
    virtual ~DetectorConstruction();
    
    virtual G4VPhysicalVolume* Construct();

private:
    G4LogicalVolume* fDetecteurRP_log;
    
    // ═══════════════════════════════════════════════════════════════════
    // PARAMÈTRES DE LA PLAQUE W/PETG AVEC CAVITÉ ET BILLES
    // ═══════════════════════════════════════════════════════════════════
    
    // Dimensions de la plaque W/PETG
    G4double fPlaqueX;          // Largeur en X (10 cm)
    G4double fPlaqueY;          // Largeur en Y (10 cm)
    G4double fPlaqueZ;          // Épaisseur totale en Z (18 mm)
    G4double fPlaquePosZ;       // Position du centre en Z
    
    // ═══════════════════════════════════════════════════════════════════
    // CAVITÉ D'AIR POUR L'EMPILEMENT DE BILLES
    // ═══════════════════════════════════════════════════════════════════
    G4double fCaviteX;          // Largeur cavité en X (50 mm)
    G4double fCaviteY;          // Largeur cavité en Y (50 mm)
    G4double fCaviteZ;          // Hauteur cavité en Z (épaisseur empilement)
    
    // ═══════════════════════════════════════════════════════════════════
    // EMPILEMENT DE BILLES DE BISMUTH
    // ═══════════════════════════════════════════════════════════════════
    G4double fBilleDiameter;    // Diamètre des billes (2 mm)
    G4int fNombrePlans;         // Nombre de plans (5)
    
    // Volumes logiques pour les billes
    G4LogicalVolume* fBilleLog;
    std::vector<G4VPhysicalVolume*> fBillePhys;

    // ═══════════════════════════════════════════════════════════════════
    // DÉTECTEUR DOSE
    // ═══════════════════════════════════════════════════════════════════
    G4double fSourcePosZ;               // Position Z de la source (2 cm)
    G4double fDistanceSourceDetecteur;  // Distance source-détecteur (18 cm)

    // ═══════════════════════════════════════════════════════════════════
    // MATÉRIAUX
    // ═══════════════════════════════════════════════════════════════════
    G4Material* fPETG;          // Matériau PETG
    G4Material* fTungsten;      // Matériau Tungstène (W)
    G4Material* fW_PETG;        // Mélange W/PETG 75%/25%
    G4Material* fBismuth;       // Matériau Bismuth solide (pour les billes)
    G4Material* fAir;           // Air (pour la cavité)
    
    // ═══════════════════════════════════════════════════════════════════
    // MÉTHODE PRIVÉE POUR CRÉER L'EMPILEMENT DE BILLES
    // ═══════════════════════════════════════════════════════════════════
    G4int CreateBillesEmpilement(G4LogicalVolume* motherVolume,
                                  G4Material* billeMaterial,
                                  G4double containerLx,
                                  G4double containerLy,
                                  G4double baseZ);
};

#endif
