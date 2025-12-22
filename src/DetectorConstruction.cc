#include "DetectorConstruction.hh"

#include "G4RunManager.hh"
#include "G4NistManager.hh"
#include "G4Box.hh"
#include "G4Orb.hh"
#include "G4Sphere.hh"
#include "G4SubtractionSolid.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4SystemOfUnits.hh"
#include "G4VisAttributes.hh"
#include "G4Colour.hh"
#include "G4UserLimits.hh"
#include "G4UnitsTable.hh"

#include <cmath>
#include <vector>
#include <iomanip>

DetectorConstruction::DetectorConstruction()
: G4VUserDetectorConstruction(),
  fDetecteurRP_log(nullptr),
  // ═══════════════════════════════════════════════════════════════════
  // PLAQUE W/PETG DE 18 mm AVEC CAVITÉ CENTRALE
  // ═══════════════════════════════════════════════════════════════════
  fPlaqueX(10.0*cm),              // Largeur en X : 10 cm
  fPlaqueY(10.0*cm),              // Largeur en Y : 10 cm
  fPlaqueZ(18.0*mm),              // Épaisseur totale : 18 mm
  fPlaquePosZ(4.65*cm),           // Centre à z = 4.65 cm (INCHANGÉ)
  // ═══════════════════════════════════════════════════════════════════
  // CAVITÉ D'AIR POUR L'EMPILEMENT DE BILLES
  // ═══════════════════════════════════════════════════════════════════
  fCaviteX(50.0*mm),              // Largeur cavité : 50 mm
  fCaviteY(50.0*mm),              // Largeur cavité : 50 mm
  fCaviteZ(0.0*mm),               // Sera calculé selon l'empilement
  // ═══════════════════════════════════════════════════════════════════
  // EMPILEMENT DE BILLES
  // ═══════════════════════════════════════════════════════════════════
  fBilleDiameter(2.0*mm),         // Diamètre billes : 2 mm
  fNombrePlans(5),                // 5 plans de billes
  fBilleLog(nullptr),
  // ═══════════════════════════════════════════════════════════════════
  // DÉTECTEUR DOSE - DISTANCE DE 18 cm (détecteur à z = 20 cm)
  // ═══════════════════════════════════════════════════════════════════
  fSourcePosZ(2.0*cm),            // Position de la source : z = 2 cm
  fDistanceSourceDetecteur(18.0*cm), // Distance source-détecteur : 18 cm
  // ═══════════════════════════════════════════════════════════════════
  // MATÉRIAUX
  // ═══════════════════════════════════════════════════════════════════
  fPETG(nullptr),
  fTungsten(nullptr),
  fW_PETG(nullptr),
  fBismuth(nullptr),
  fAir(nullptr)
{ 
    // Calcul de la hauteur de la cavité (empilement de 5 plans)
    G4double bille_radius = fBilleDiameter / 2.0;
    G4double dz_compact = fBilleDiameter * std::sqrt(2.0/3.0);  // 1.633 mm
    // Hauteur = rayon + 4*dz + rayon = diamètre + 4*dz
    fCaviteZ = fBilleDiameter + 4.0 * dz_compact;  // ~8.53 mm
}

DetectorConstruction::~DetectorConstruction()
{ }

// =============================================================================
// MÉTHODE POUR CRÉER L'EMPILEMENT DE BILLES
// =============================================================================
G4int DetectorConstruction::CreateBillesEmpilement(
    G4LogicalVolume* motherVolume,
    G4Material* billeMaterial,
    G4double containerLx,
    G4double containerLy,
    G4double baseZ)
{
    // =========================================================================
    // PARAMÈTRES GÉOMÉTRIQUES
    // =========================================================================
    
    G4double d = fBilleDiameter;
    G4double r = d / 2.0;
    
    // Constantes d'empilement hexagonal compact
    G4double dy_hex = d * std::sqrt(3.0) / 2.0;   // 1.732 mm
    G4double dz = d * std::sqrt(2.0 / 3.0);       // 1.633 mm
    G4double dx_B = d / 2.0;                      // 1.0 mm
    G4double dy_B = d / (2.0 * std::sqrt(3.0));   // 0.577 mm
    
    // Positions Z des 5 plans (relatives à baseZ)
    G4double z[5];
    for (G4int i = 0; i < 5; i++) {
        z[i] = baseZ + r + i * dz;
    }
    
    // =========================================================================
    // CRÉATION DU VOLUME BILLE
    // =========================================================================
    
    G4Orb* solidBille = new G4Orb("Bille", r);
    fBilleLog = new G4LogicalVolume(solidBille, billeMaterial, "BilleLog");
    
    // Visualisation (doré pour le bismuth)
    G4VisAttributes* billeVis = new G4VisAttributes(G4Colour(0.8, 0.6, 0.2, 0.9));
    billeVis->SetForceSolid(true);
    fBilleLog->SetVisAttributes(billeVis);
    
    // =========================================================================
    // GÉNÉRATION DES BILLES
    // =========================================================================
    
    G4int totalBilles = 0;
    G4int billesParPlan[5] = {0, 0, 0, 0, 0};
    
    for (G4int plan = 0; plan < 5; plan++) {
        
        // Type de plan : A (0,2,4) ou B (1,3)
        G4bool typeB = (plan == 1 || plan == 3);
        
        // Décalages XY pour ce plan
        G4double offX = typeB ? dx_B : 0.0;
        G4double offY = typeB ? dy_B : 0.0;
        
        // Nombre de rangées
        G4double y0 = r + offY;
        G4int nRows = (G4int)std::floor((containerLy - r - y0) / dy_hex) + 1;
        
        for (G4int row = 0; row < nRows; row++) {
            
            G4double yc = y0 + row * dy_hex;
            if (yc + r > containerLy) continue;
            
            // Décalage X alterné dans le plan hexagonal
            G4double rowOffX = (row % 2 == 1) ? r : 0.0;
            G4double x0 = r + offX + rowOffX;
            
            // Parcours des positions X
            for (G4double xc = x0; xc + r <= containerLx; xc += d) {
                
                if (xc - r < 0) continue;
                
                // Position relative au centre de la cavité
                // La cavité est centrée en (0,0), donc on décale
                G4double xPos = xc - containerLx/2.0;
                G4double yPos = yc - containerLy/2.0;
                G4double zPos = z[plan];
                
                G4ThreeVector pos(xPos, yPos, zPos);
                
                G4String name = "Bille_P" + std::to_string(plan+1) + 
                               "_" + std::to_string(totalBilles);
                
                G4VPhysicalVolume* physBille = new G4PVPlacement(
                    0,              // Rotation
                    pos,            // Position
                    fBilleLog,      // Volume logique
                    name,           // Nom
                    motherVolume,   // Volume mère (cavité d'air)
                    false,          // Boolean
                    totalBilles,    // Copy number
                    false           // Check overlaps (false pour vitesse)
                );
                
                fBillePhys.push_back(physBille);
                totalBilles++;
                billesParPlan[plan]++;
            }
        }
    }
    
    // Affichage récapitulatif
    G4cout << "\n+===========================================================+" << G4endl;
    G4cout << "|      EMPILEMENT DE 5 PLANS DE BILLES DE BISMUTH           |" << G4endl;
    G4cout << "+===========================================================+" << G4endl;
    G4cout << "|  Diametre billes : " << d/mm << " mm" << G4endl;
    G4cout << "|  Container       : " << containerLx/mm << " x " << containerLy/mm << " mm2" << G4endl;
    G4cout << "|  Materiau        : " << billeMaterial->GetName() << G4endl;
    G4cout << "+-----------------------------------------------------------+" << G4endl;
    
    for (G4int p = 0; p < 5; p++) {
        G4bool tB = (p == 1 || p == 3);
        G4cout << "|  Plan " << p+1 << " (Type " << (tB ? "B" : "A") 
               << ") : " << std::setw(4) << billesParPlan[p] 
               << " billes, z = " << std::fixed << std::setprecision(3) 
               << z[p]/mm << " mm (relatif)" << G4endl;
    }
    
    G4cout << "+===========================================================+" << G4endl;
    G4cout << "|  TOTAL           : " << totalBilles << " billes" << G4endl;
    G4cout << "|  Epaisseur       : " << (z[4] + r - baseZ)/mm << " mm" << G4endl;
    G4cout << "+===========================================================+\n" << G4endl;
    
    return totalBilles;
}

// =============================================================================
// MÉTHODE CONSTRUCT
// =============================================================================
G4VPhysicalVolume* DetectorConstruction::Construct()
{
    G4NistManager* nist = G4NistManager::Instance();

    // =============================================================================
    // MATÉRIAUX DE BASE
    // =============================================================================

    fAir = nist->FindOrBuildMaterial("G4_AIR");
    G4Material* water = nist->FindOrBuildMaterial("G4_WATER");

    // =============================================================================
    // TUNGSTÈNE (W)
    // =============================================================================
    fTungsten = nist->FindOrBuildMaterial("G4_W");

    // =============================================================================
    // PETG (approximation PET)
    // =============================================================================
    G4Element* elC = nist->FindOrBuildElement("C");
    G4Element* elH = nist->FindOrBuildElement("H");
    G4Element* elO = nist->FindOrBuildElement("O");

    fPETG = new G4Material("PETG", 1.27*g/cm3, 3, kStateSolid);
    fPETG->AddElement(elC, 10);
    fPETG->AddElement(elH,  8);
    fPETG->AddElement(elO,  4);

    // =============================================================================
    // Mélange W/PETG : 75%/25% (fractions massiques)
    // =============================================================================
    G4double rhoW = fTungsten->GetDensity();
    G4double rhoPETG = fPETG->GetDensity();
    G4double massFracW = 0.75;
    G4double massFracPETG = 0.25;

    G4double rhoMix = 1.0 / (massFracW/rhoW + massFracPETG/rhoPETG);

    fW_PETG = new G4Material("W_PETG_75_25", rhoMix, 2, kStateSolid);
    fW_PETG->AddMaterial(fTungsten, massFracW);
    fW_PETG->AddMaterial(fPETG, massFracPETG);

    // =============================================================================
    // BISMUTH SOLIDE (pour les billes)
    // =============================================================================
    fBismuth = nist->FindOrBuildMaterial("G4_Bi");

    G4cout << "\n=== Proprietes des Materiaux ===" << G4endl;
    G4cout << "W/PETG 75%/25%:" << G4endl;
    G4cout << "  W density      = " << G4BestUnit(rhoW, "Volumic Mass") << G4endl;
    G4cout << "  PETG density   = " << G4BestUnit(rhoPETG, "Volumic Mass") << G4endl;
    G4cout << "  Mix density    = " << G4BestUnit(rhoMix, "Volumic Mass") << G4endl;
    G4cout << "Bismuth (billes):" << G4endl;
    G4cout << "  Bi density     = " << G4BestUnit(fBismuth->GetDensity(), "Volumic Mass") << G4endl;
    G4cout << "================================\n" << G4endl;

    // =============================================================================
    // WORLD (MONDE)
    // =============================================================================

    G4double world_size = 400*cm;
    G4Box* solidWorld = new G4Box("World",
                                  world_size/2,
                                  world_size/2,
                                  world_size/2);

    G4LogicalVolume* logicWorld = new G4LogicalVolume(solidWorld, fAir, "World");

    G4VPhysicalVolume* physWorld = new G4PVPlacement(0,
                                                     G4ThreeVector(),
                                                     logicWorld,
                                                     "World",
                                                     0,
                                                     false,
                                                     0);

    logicWorld->SetVisAttributes(G4VisAttributes::GetInvisible());

    // ===================================================================
    // ENVELOPPE
    // ===================================================================
    G4Box* solidEnveloppe = new G4Box("Enveloppe", 40*cm, 40*cm, 40*cm);

    G4LogicalVolume* logicEnveloppe = new G4LogicalVolume(solidEnveloppe, fAir, "Enveloppe");

    new G4PVPlacement(nullptr,
                      G4ThreeVector(),
                      logicEnveloppe,
                      "Enveloppe",
                      logicWorld,
                      false,
                      0,
                      true);

    G4VisAttributes* enveloppeVis = new G4VisAttributes(G4Colour(1.0, 1.0, 1.0, 0.1));
    enveloppeVis->SetVisibility(true);
    logicEnveloppe->SetVisAttributes(enveloppeVis);

    // ═══════════════════════════════════════════════════════════════════
    // PLAQUE W/PETG DE 18 mm AVEC CAVITÉ CENTRALE CONTENANT LES BILLES
    // ═══════════════════════════════════════════════════════════════════

    G4double plaque_front_z = fPlaquePosZ - fPlaqueZ/2;  // 3.75 cm
    G4double plaque_back_z = fPlaquePosZ + fPlaqueZ/2;   // 5.55 cm

    // ───────────────────────────────────────────────────────────────
    // CRÉATION DE LA PLAQUE W/PETG AVEC CAVITÉ (soustraction booléenne)
    // ───────────────────────────────────────────────────────────────
    
    // Solide plein de la plaque
    G4Box* solidPlaqueFull = new G4Box("PlaqueFull",
                                        fPlaqueX/2,
                                        fPlaqueY/2,
                                        fPlaqueZ/2);
    
    // Solide de la cavité (centrée dans la plaque)
    G4Box* solidCavite = new G4Box("Cavite",
                                    fCaviteX/2,
                                    fCaviteY/2,
                                    fCaviteZ/2);
    
    // Soustraction : plaque - cavité (cavité centrée)
    G4SubtractionSolid* solidPlaque = new G4SubtractionSolid(
        "PlaqueWPETG",
        solidPlaqueFull,
        solidCavite,
        0,                      // Pas de rotation
        G4ThreeVector(0, 0, 0)  // Cavité centrée
    );
    
    G4LogicalVolume* logicPlaque = new G4LogicalVolume(solidPlaque,
                                                        fW_PETG,  // fAir
                                                        "PlaqueWPETGLog");
    
    new G4PVPlacement(nullptr,
                      G4ThreeVector(0, 0, fPlaquePosZ),
                      logicPlaque,
                      "PlaqueWPETG",
                      logicEnveloppe,
                      false,
                      0,
                      true);
    
    // Attributs visuels (vert pour W/PETG)
    G4VisAttributes* wpetgVis = new G4VisAttributes(G4Colour(0.0, 0.8, 0.0, 0.5));
    wpetgVis->SetForceSolid(true);
    logicPlaque->SetVisAttributes(wpetgVis);

    // ───────────────────────────────────────────────────────────────
    // CAVITÉ D'AIR (volume fille de l'enveloppe, positionnée au centre de la plaque)
    // ───────────────────────────────────────────────────────────────
    
    G4LogicalVolume* logicCavite = new G4LogicalVolume(solidCavite,
                                                        fAir,
                                                        "CaviteAirLog");
    
    new G4PVPlacement(nullptr,
                      G4ThreeVector(0, 0, fPlaquePosZ),  // Même centre que la plaque
                      logicCavite,
                      "CaviteAir",
                      logicEnveloppe,
                      false,
                      0,
                      true);
    
    // Attributs visuels (bleu clair transparent pour l'air)
    G4VisAttributes* caviteVis = new G4VisAttributes(G4Colour(0.7, 0.9, 1.0, 0.3));
    caviteVis->SetForceSolid(true);
    logicCavite->SetVisAttributes(caviteVis);

    // ───────────────────────────────────────────────────────────────
    // EMPILEMENT DE BILLES DE BISMUTH DANS LA CAVITÉ
    // ───────────────────────────────────────────────────────────────
    
    G4double empilement_base_z = -fCaviteZ/2;
    
    G4int nombreBilles = CreateBillesEmpilement(
        logicCavite,
        fBismuth,  // fAir
        fCaviteX,
        fCaviteY,
        empilement_base_z
    );

    // =============================================================================
    // PLANS DE COMPTAGE (UPSTREAM ET DOWNSTREAM)
    // =============================================================================

    G4double detector_thickness = 2.0*mm;
    G4double detector_gap = 2.0*mm;
    G4double detector_size = 12.0*cm;

    G4double upstream_detector_z = plaque_front_z - detector_gap - detector_thickness/2.0;
    G4double downstream_detector_z = plaque_back_z + detector_gap + detector_thickness/2.0;

    G4Box* solidDetector = new G4Box("Detector",
                                     detector_size/2.0,
                                     detector_size/2.0,
                                     detector_thickness/2.0);

    // ───────────────────────────────────────────────────────────────
    // UPSTREAM DETECTOR
    // ───────────────────────────────────────────────────────────────

    G4LogicalVolume* logicUpstreamDetector = new G4LogicalVolume(solidDetector, fAir, "UpstreamDetector");

    G4VisAttributes* upstreamVisAtt = new G4VisAttributes(G4Colour(0.0, 0.0, 1.0, 0.3));
    upstreamVisAtt->SetForceSolid(true);
    logicUpstreamDetector->SetVisAttributes(upstreamVisAtt);

    new G4PVPlacement(0,
                      G4ThreeVector(0., 0., upstream_detector_z),
                      logicUpstreamDetector,
                      "UpstreamDetector",
                      logicEnveloppe,
                      false,
                      0,
                      false);

    // ───────────────────────────────────────────────────────────────
    // DOWNSTREAM DETECTOR
    // ───────────────────────────────────────────────────────────────

    G4LogicalVolume* logicDownstreamDetector = new G4LogicalVolume(solidDetector, fAir, "DownstreamDetector");

    G4VisAttributes* downstreamVisAtt = new G4VisAttributes(G4Colour(1.0, 0.0, 0.0, 0.3));
    downstreamVisAtt->SetForceSolid(true);
    logicDownstreamDetector->SetVisAttributes(downstreamVisAtt);

    new G4PVPlacement(0,
                      G4ThreeVector(0., 0., downstream_detector_z),
                      logicDownstreamDetector,
                      "DownstreamDetector",
                      logicEnveloppe,
                      false,
                      0,
                      false);

    // =============================================================================
    // DÉTECTEUR DOSE (sphère d'EAU) - DISTANCE DE 18 cm (détecteur à z = 20 cm)
    // =============================================================================

    G4double dose_position_z = fSourcePosZ + fDistanceSourceDetecteur;  // 2 cm + 18 cm = 20 cm
    G4double dose_radius = 2.0*cm;

    G4Orb* doseDetector_solid = new G4Orb("DoseDetector", dose_radius);

    G4LogicalVolume* doseDetector_log = new G4LogicalVolume(
        doseDetector_solid,
        water,
        "DoseDetectorLog");

    G4double maxStep = 0.1*mm;
    G4UserLimits* doseLimits = new G4UserLimits(maxStep);
    doseDetector_log->SetUserLimits(doseLimits);

    G4VisAttributes* doseVis = new G4VisAttributes(G4Colour(0.0, 0.5, 1.0, 0.5));
    doseVis->SetForceSolid(true);
    doseDetector_log->SetVisAttributes(doseVis);

    G4ThreeVector dose_position = G4ThreeVector(0, 0, dose_position_z);

    new G4PVPlacement(
        0,
        dose_position,
        doseDetector_log,
        "DoseDetector",
        logicEnveloppe,
        false,
        0,
        true);

    // =============================================================================
    // CALCULS POUR L'AFFICHAGE
    // =============================================================================
    
    G4double dose_volume = (4.0/3.0) * M_PI * std::pow(dose_radius, 3);
    G4double water_density = water->GetDensity();
    G4double dose_mass = dose_volume * water_density;

    // Volume et masse de la plaque W/PETG (avec cavité)
    G4double plaque_volume_full = fPlaqueX * fPlaqueY * fPlaqueZ;
    G4double cavite_volume = fCaviteX * fCaviteY * fCaviteZ;
    G4double plaque_volume = plaque_volume_full - cavite_volume;
    G4double plaque_mass = plaque_volume * fW_PETG->GetDensity();

    // Volume et masse des billes
    G4double bille_radius = fBilleDiameter / 2.0;
    G4double bille_volume = (4.0/3.0) * M_PI * std::pow(bille_radius, 3);
    G4double billes_volume_total = nombreBilles * bille_volume;
    G4double billes_mass = billes_volume_total * fBismuth->GetDensity();

    // Masse totale
    G4double total_mass = plaque_mass + billes_mass;

    // Masse surfacique
    G4double surface = fPlaqueX * fPlaqueY;
    G4double plaque_surfacic_mass = plaque_mass / surface;
    G4double billes_surfacic_mass = billes_mass / surface;
    G4double total_surfacic_mass = total_mass / surface;

    // Distance plaque - détecteur
    G4double distance_plaque_detecteur = dose_position_z - plaque_back_z;

    // =============================================================================
    // AFFICHAGE RÉCAPITULATIF
    // =============================================================================

    G4cout << "\n+=====================================================================+" << G4endl;
    G4cout << "|      PROGRAMME : MESURE DE DOSE DANS L'EAU                          |" << G4endl;
    G4cout << "|      (Version PLAQUE W/PETG + BILLES DE BISMUTH)                    |" << G4endl;
    G4cout << "|      >>> DISTANCE SOURCE-DETECTEUR = 18 cm <<<                      |" << G4endl;
    G4cout << "+=====================================================================+" << G4endl;
    G4cout << "|                                                                     |" << G4endl;
    G4cout << "| STRUCTURE :                                                         |" << G4endl;
    G4cout << "| +-------------------------------------------------------------------+" << G4endl;
    G4cout << "| |                      PLAQUE W/PETG 18 mm                          |" << G4endl;
    G4cout << "| |  +-----------------------------------------------------+          |" << G4endl;
    G4cout << "| |  |     CAVITE D'AIR avec BILLES DE BISMUTH             |          |" << G4endl;
    G4cout << "| |  |     (50x50 mm, " << std::fixed << std::setprecision(2) << fCaviteZ/mm << " mm hauteur)                    |          |" << G4endl;
    G4cout << "| |  |     " << nombreBilles << " billes de " << fBilleDiameter/mm << " mm (5 plans HCP)            |          |" << G4endl;
    G4cout << "| |  +-----------------------------------------------------+          |" << G4endl;
    G4cout << "| +-------------------------------------------------------------------+" << G4endl;
    G4cout << "|                                                                     |" << G4endl;
    G4cout << "+---------------------------------------------------------------------+" << G4endl;
    G4cout << "| POSITIONS SUR L'AXE Z :                                             |" << G4endl;
    G4cout << "| * Source             : z = " << fSourcePosZ/cm << " cm                                 |" << G4endl;
    G4cout << "| * Face avant plaque  : z = " << plaque_front_z/cm << " cm                              |" << G4endl;
    G4cout << "| * Centre plaque      : z = " << fPlaquePosZ/cm << " cm                              |" << G4endl;
    G4cout << "| * Face arriere plaque: z = " << plaque_back_z/cm << " cm                              |" << G4endl;
    G4cout << "| * Centre detecteur   : z = " << dose_position_z/cm << " cm                                |" << G4endl;
    G4cout << "|                                                                     |" << G4endl;
    G4cout << "+---------------------------------------------------------------------+" << G4endl;
    G4cout << "| DISTANCES :                                                         |" << G4endl;
    G4cout << "| * Source -> Detecteur    : " << fDistanceSourceDetecteur/cm << " cm                              |" << G4endl;
    G4cout << "| * Plaque -> Detecteur    : " << distance_plaque_detecteur/cm << " cm                            |" << G4endl;
    G4cout << "|                                                                     |" << G4endl;
    G4cout << "+---------------------------------------------------------------------+" << G4endl;
    G4cout << "| PLAQUE W/PETG :                                                     |" << G4endl;
    G4cout << "| * Dimensions X x Y   : " << fPlaqueX/cm << " x " << fPlaqueY/cm << " cm                            |" << G4endl;
    G4cout << "| * Epaisseur totale   : " << fPlaqueZ/mm << " mm                                     |" << G4endl;
    G4cout << "| * Densite W/PETG     : " << fW_PETG->GetDensity()/(g/cm3) << " g/cm3                         |" << G4endl;
    G4cout << "| * Masse W/PETG       : " << plaque_mass/g << " g                               |" << G4endl;
    G4cout << "| * Masse surf. W/PETG : " << plaque_surfacic_mass/(g/cm2) << " g/cm2                          |" << G4endl;
    G4cout << "|                                                                     |" << G4endl;
    G4cout << "+---------------------------------------------------------------------+" << G4endl;
    G4cout << "| CAVITE ET BILLES :                                                  |" << G4endl;
    G4cout << "| * Dimensions cavite  : " << fCaviteX/mm << " x " << fCaviteY/mm << " x " << fCaviteZ/mm << " mm3                  |" << G4endl;
    G4cout << "| * Nombre de billes   : " << nombreBilles << "                                        |" << G4endl;
    G4cout << "| * Diametre billes    : " << fBilleDiameter/mm << " mm                                       |" << G4endl;
    G4cout << "| * Materiau billes    : Bismuth (Z=83)                               |" << G4endl;
    G4cout << "| * Densite Bi         : " << fBismuth->GetDensity()/(g/cm3) << " g/cm3                            |" << G4endl;
    G4cout << "| * Masse billes       : " << billes_mass/g << " g                                |" << G4endl;
    G4cout << "| * Masse surf. billes : " << billes_surfacic_mass/(g/cm2) << " g/cm2                           |" << G4endl;
    G4cout << "|                                                                     |" << G4endl;
    G4cout << "+---------------------------------------------------------------------+" << G4endl;
    G4cout << "| TOTAL :                                                             |" << G4endl;
    G4cout << "| * Masse totale       : " << total_mass/g << " g                               |" << G4endl;
    G4cout << "| * Masse surfacique   : " << total_surfacic_mass/(g/cm2) << " g/cm2                          |" << G4endl;
    G4cout << "|                                                                     |" << G4endl;
    G4cout << "+---------------------------------------------------------------------+" << G4endl;
    G4cout << "| PLANS DE COMPTAGE :                                                 |" << G4endl;
    G4cout << "| * Upstream           : z = " << upstream_detector_z/cm << " cm                             |" << G4endl;
    G4cout << "| * Downstream         : z = " << downstream_detector_z/cm << " cm                             |" << G4endl;
    G4cout << "|                                                                     |" << G4endl;
    G4cout << "+---------------------------------------------------------------------+" << G4endl;
    G4cout << "| DETECTEUR DOSE DANS L'EAU :                                         |" << G4endl;
    G4cout << "| * Position centre    : z = " << dose_position_z/cm << " cm                               |" << G4endl;
    G4cout << "| * Distance source    : " << fDistanceSourceDetecteur/cm << " cm                               |" << G4endl;
    G4cout << "| * Rayon              : " << dose_radius/cm << " cm                                     |" << G4endl;
    G4cout << "| * Masse              : " << dose_mass/g << " g                                  |" << G4endl;
    G4cout << "|                                                                     |" << G4endl;
    G4cout << "+=====================================================================+\n" << G4endl;

    return physWorld;
}
