#ifndef RunAction_h
#define RunAction_h 1

#include "G4UserRunAction.hh"
#include "globals.hh"
#include <vector>

class G4Run;

/// @brief Gestion des runs avec histogrammes et statistiques
///
/// Cette classe crée les histogrammes et ntuples pour enregistrer
/// les statistiques de la simulation. Elle accumule les données
/// envoyées par EventAction à chaque événement.

class RunAction : public G4UserRunAction
{
public:
    RunAction();
    virtual ~RunAction();
    
    virtual void BeginOfRunAction(const G4Run*);
    virtual void EndOfRunAction(const G4Run*);
    
    // ═══════════════════════════════════════════════════════════════
    // MÉTHODES POUR ACCUMULER LES DONNÉES (appelées par EventAction)
    // ═══════════════════════════════════════════════════════════════
    
    /// Ajoute l'énergie déposée dans le détecteur kerma (Méthode 1)
    void AddKermaEnergy(G4double edep);
    
    /// Ajoute la contribution fluence au kerma (Méthode 2)
    void AddKermaFluence(G4double fluenceContrib);
    
    /// Ajoute l'énergie forcée au kerma (Méthode 1bis)
    void AddKermaEnergyForced(G4double forcedDeposit);
    
    /// Enregistre les statistiques d'un événement
    void RecordEventStatistics(G4int nPrimaries, 
                               const std::vector<G4double>& primaryEnergies,
                               G4int nTransmitted,
                               G4int nAbsorbed,
                               G4double kermaDeposit);

private:
    // ═══════════════════════════════════════════════════════════════
    // DONNÉES POUR LE CALCUL DU KERMA
    // ═══════════════════════════════════════════════════════════════
    G4double fKermaTotalEnergy;     // Somme des énergies déposées (MeV) - Méthode 1
    G4double fKermaTotalEnergy2;    // Somme des énergies² (pour variance)
    G4double fKermaMass;            // Masse du détecteur d'air (kg)
    G4double fKermaRadius;          // Rayon du détecteur (cm)
    G4double fKermaPosition;        // Position z du détecteur (cm)
    G4int fKermaEventCount;         // Nombre d'événements avec dépôt
    
    // Données pour la méthode 2 (fluence)
    G4double fKermaFluenceTotal;    // Somme de E × μ_en/ρ (keV·cm²/g)
    G4int fKermaFluenceCount;       // Nombre de gammas traversant le détecteur
    
    // Données pour la méthode 1bis (forçage d'interaction)
    G4double fKermaForcedTotal;     // Somme des énergies forcées déposées (MeV)
    G4int fKermaForcedCount;        // Nombre de gammas avec forçage
    
    // Paramètres source pour normalisation du débit de Kerma
    G4double fActivity4pi;          // Activité source 4π (Bq)
    G4double fConeAngle;            // Demi-angle du cône d'émission
    
    // Paramètres géométriques pour normalisation améliorée
    G4double fSourcePosZ;           // Position Z de la source
    G4double fDetectorPosZ;         // Position Z du détecteur
    G4double fMeanGammasPerDecay;   // Nombre moyen de gammas par désintégration

    // ═══════════════════════════════════════════════════════════════
    // COMPTEURS POUR STATISTIQUES GLOBALES
    // ═══════════════════════════════════════════════════════════════
    G4int fTotalPrimariesGenerated;   // Total de gammas primaires générés
    G4int fTotalEventsWithZeroGamma;  // Événements sans gamma (désintégration vide)
    G4int fTotalTransmitted;          // Total de gammas transmis
    G4int fTotalAbsorbed;             // Total de gammas absorbés
    
    // ═══════════════════════════════════════════════════════════════
    // NOM DU FICHIER DE SORTIE
    // ═══════════════════════════════════════════════════════════════
    G4String fOutputFileName;
};

#endif
