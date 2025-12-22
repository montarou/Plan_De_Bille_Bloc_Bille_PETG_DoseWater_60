#include "PhysicsList.hh"

#include "G4EmLivermorePhysics.hh"
#include "G4EmStandardPhysics_option4.hh"
#include "G4DecayPhysics.hh"
#include "G4StepLimiterPhysics.hh"
#include "G4SystemOfUnits.hh"

PhysicsList::PhysicsList() : FTFP_BERT() {
  

  // Modèles basse énergie optimisés pour photons/électrons de keV à quelques MeV
  ReplacePhysics(new G4EmLivermorePhysics());

  // Physique de décroissance
  RegisterPhysics(new G4DecayPhysics());
  
  // Step Limiter : permet d'appliquer les G4UserLimits dans les volumes
  // Nécessaire pour forcer des steps courts dans le détecteur Kerma
  RegisterPhysics(new G4StepLimiterPhysics());

  G4cout << "\n========== PHYSIQUE ==========\n"
  << "Liste de physique : FTFP_BERT\n"
  << "EM Physics : Livermore (optimisée basse énergie)\n"
  << "  - Photoélectrique avec couches atomiques\n"
  << "  - Compton avec fonction de diffusion\n"
  << "  - Diffusion Rayleigh\n"
  << "  - Production de paires\n"
  << "  - Bremsstrahlung\n"
  << "  - Ionisation\n"
  << "Step Limiter : activé (pour UserLimits)\n"
  << "==============================\n" << G4endl;
}

void PhysicsList::SetCuts() {
  
  // Cuts de production (distances minimales)
  SetCutValue(0.1*mm, "gamma");
  SetCutValue(0.1*mm, "e-");
  SetCutValue(0.1*mm, "e+");
  SetCutValue(0.1*mm, "proton");
  
  // ALTERNATIVE : Cuts encore plus fins pour précision maximale (mais plus lent)
  // Décommenter les lignes suivantes pour des cuts de 0.01 mm
  // SetCutValue(0.01*mm, "gamma");
  // SetCutValue(0.01*mm, "e-");
  // SetCutValue(0.01*mm, "e+");
  // SetCutValue(0.01*mm, "proton");

  if (verboseLevel > 0) {
    G4cout << "\n========== CUTS DE PRODUCTION ==========\n"
           << "Gamma : 0.1 mm\n"
           << "e-    : 0.1 mm\n"
           << "e+    : 0.1 mm\n"
           << "Proton: 0.1 mm\n"
           << "========================================\n" << G4endl;
  }
}
