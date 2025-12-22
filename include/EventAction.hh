#ifndef EventAction_h
#define EventAction_h 1

#include "G4UserEventAction.hh"
#include "globals.hh"
#include <vector>
#include <map>

class RunAction;

/// @brief Gestion des événements avec suivi des gammas primaires multiples
///
/// Cette classe gère le cycle de vie d'un événement et stocke les informations
/// de tous les gammas primaires générés. Elle récupère les infos depuis G4Event
/// dans BeginOfEventAction et les met à jour via SteppingAction.

class EventAction : public G4UserEventAction
{
public:
    EventAction(RunAction* runAction);
    virtual ~EventAction();
    
    virtual void BeginOfEventAction(const G4Event*);
    virtual void EndOfEventAction(const G4Event*);
    
    // ═══════════════════════════════════════════════════════════════
    // STRUCTURE POUR LES GAMMAS PRIMAIRES
    // ═══════════════════════════════════════════════════════════════
    struct PrimaryGammaInfo {
        G4int trackID;              // Identifiant unique de la trace (1, 2, 3, ...)
        G4double energyInitial;     // Énergie à la génération (MeV)
        G4double energyUpstream;    // Énergie au plan upstream (0 si pas détecté)
        G4double energyDownstream;  // Énergie au plan downstream (0 si pas détecté)
        G4double theta;             // Angle polaire initial (rad)
        G4double phi;               // Angle azimutal initial (rad)
        G4bool detectedUpstream;    // A traversé le plan upstream ?
        G4bool detectedDownstream;  // A traversé le plan downstream ?
        G4bool transmitted;         // Transmis sans perte d'énergie significative ?
    };

    // ═══════════════════════════════════════════════════════════════
    // STRUCTURE POUR LES PARTICULES SECONDAIRES DOWNSTREAM
    // ═══════════════════════════════════════════════════════════════
    struct SecondaryParticleInfo {
        G4int trackID;              // Identifiant de la trace
        G4int parentID;             // Identifiant du parent
        G4int pdgCode;              // Code PDG de la particule
        G4double energy;            // Énergie cinétique (MeV)
        G4String creatorProcess;    // Processus créateur (compt, phot, etc.)
    };

    // ═══════════════════════════════════════════════════════════════
    // MÉTHODES POUR ENREGISTRER LES PASSAGES (appelées par SteppingAction)
    // ═══════════════════════════════════════════════════════════════

    /// Enregistre le passage d'un gamma primaire au plan upstream
    void RecordPrimaryUpstream(G4int trackID, G4double energy);

    /// Enregistre le passage d'un gamma primaire au plan downstream
    void RecordPrimaryDownstream(G4int trackID, G4double energy);

    /// Enregistre une particule secondaire au plan downstream
    void RecordSecondaryDownstream(G4int trackID, G4int parentID,
                                   G4int pdgCode, G4double energy,
                                   const G4String& process);

    // ═══════════════════════════════════════════════════════════════
    // MÉTHODES POUR LE KERMA DANS L'AIR
    // ═══════════════════════════════════════════════════════════════
    void AddKermaEnergy(G4double edep);     // Méthode 1 : dépôt
    void AddKermaFluence(G4double energy, G4double chordLength);  // Méthode 2 : fluence
    void AddKermaEnergyForced(G4double energy, G4double chordLength);  // Méthode 1bis : forçage

    // ═══════════════════════════════════════════════════════════════
    // ACCESSEURS
    // ═══════════════════════════════════════════════════════════════
    /// Retourne le vecteur des gammas primaires
    const std::vector<PrimaryGammaInfo>& GetPrimaryGammas() const {
        return fPrimaryGammas;
    }

    /// Retourne le vecteur des secondaires downstream
    const std::vector<SecondaryParticleInfo>& GetSecondariesDownstream() const {
        return fSecondariesDownstream;
    }

    /// Nombre de gammas primaires dans cet événement
    G4int GetNumberOfPrimaries() const { return fPrimaryGammas.size(); }

    /// Nombre de gammas transmis
    G4int GetNumberTransmitted() const;

    /// Nombre de gammas absorbés (détectés upstream mais pas downstream)
    G4int GetNumberAbsorbed() const;

    /// Vérifie si un trackID correspond à un primaire
    G4bool IsPrimaryTrack(G4int trackID) const;

private:
    RunAction* fRunAction;

    // ═══════════════════════════════════════════════════════════════
    // STOCKAGE DES INFORMATIONS
    // ═══════════════════════════════════════════════════════════════

    /// Vecteur des gammas primaires de cet événement
    std::vector<PrimaryGammaInfo> fPrimaryGammas;

    /// Vecteur des particules secondaires détectées downstream
    std::vector<SecondaryParticleInfo> fSecondariesDownstream;

    /// Map pour accès rapide : trackID → index dans fPrimaryGammas
    std::map<G4int, size_t> fTrackIDtoIndex;

    /// Énergie déposée dans le détecteur kerma (Méthode 1)
    G4double fKermaEnergyDeposit;

    /// Contribution fluence au kerma (Méthode 2) : Σ(E × μ_en/ρ)
    G4double fKermaFluence;

    /// Énergie forcée déposée (Méthode 1bis) : Σ(E × L × μ_en × ρ)
    G4double fKermaForcedDeposit;

    // ═══════════════════════════════════════════════════════════════
    // PARAMÈTRES
    // ═══════════════════════════════════════════════════════════════

    /// Tolérance pour considérer un gamma comme "transmis" (pas de perte)
    G4double fTransmissionTolerance;

    /// Niveau de verbosité (0=silencieux, 1=résumé, 2=détaillé)
    G4int fVerboseLevel;
};

#endif
