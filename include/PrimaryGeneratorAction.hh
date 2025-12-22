#ifndef PrimaryGeneratorAction_h
#define PrimaryGeneratorAction_h 1

#include "G4VUserPrimaryGeneratorAction.hh"
#include "G4ParticleGun.hh"
#include "G4ThreeVector.hh"
#include "globals.hh"
#include <vector>

class G4ParticleGun;
class G4Event;

/// @brief Génération des particules primaires selon le spectre Eu-152
///
/// Cette classe génère des gammas selon le spectre de l'Europium-152.
/// Plusieurs gammas peuvent être émis par événement (désintégration).
/// Les informations sont stockées automatiquement dans G4Event et
/// récupérées par EventAction::BeginOfEventAction().

class PrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction
{
public:
    PrimaryGeneratorAction();    
    virtual ~PrimaryGeneratorAction();
    
    virtual void GeneratePrimaries(G4Event*);
    
    const G4ParticleGun* GetParticleGun() const { return fParticleGun; }

    // ═══════════════════════════════════════════════════════════════
    // ACCESSEURS POUR DIAGNOSTIC
    // ═══════════════════════════════════════════════════════════════
    G4int GetLastEventGammaCount() const { return fLastEventGammaCount; }

    // Accès au spectre (pour vérification)
    const std::vector<G4double>& GetGammaEnergies() const { return fGammaEnergies; }
    const std::vector<G4double>& GetGammaProbabilities() const { return fGammaProbabilities; }
    
private:
    G4ParticleGun* fParticleGun;

    // ═══════════════════════════════════════════════════════════════
    // SPECTRE GAMMA Europium-152
    // ═══════════════════════════════════════════════════════════════
    std::vector<G4double> fGammaEnergies;        // Énergies des raies (keV)
    std::vector<G4double> fGammaIntensities;     // Intensités relatives (%)
    std::vector<G4double> fGammaProbabilities;   // Probabilités d'émission (fraction 0-1)

    // ═══════════════════════════════════════════════════════════════
    // COMPTEUR POUR LE DERNIER ÉVÉNEMENT
    // ═══════════════════════════════════════════════════════════════
    G4int fLastEventGammaCount;

    // ═══════════════════════════════════════════════════════════════
    // PARAMÈTRES DE LA SOURCE
    // ═══════════════════════════════════════════════════════════════
    G4double fConeAngle;        // Demi-angle du cône d'émission
    G4ThreeVector fSourcePosition;  // Position de la source

    // ═══════════════════════════════════════════════════════════════
    // MÉTHODE POUR GÉNÉRER UNE DIRECTION DANS LE CÔNE
    // ═══════════════════════════════════════════════════════════════
    void GenerateDirectionInCone(G4double coneAngle,
                                 G4double& theta,
                                 G4double& phi,
                                 G4ThreeVector& direction);
};

#endif
