#ifndef SteppingAction_h
#define SteppingAction_h 1

#include "G4UserSteppingAction.hh"
#include "globals.hh"
#include <vector>

class EventAction;

/// @brief Suivi pas-à-pas des particules
///
/// Cette classe détecte les passages des particules dans les volumes
/// de détection (upstream, downstream, kerma) et met à jour les
/// structures de EventAction.
///
/// Identification des primaires : parentID == 0

class SteppingAction : public G4UserSteppingAction
{
public:
    SteppingAction(EventAction* eventAction);
    virtual ~SteppingAction();
    
    virtual void UserSteppingAction(const G4Step*);
    
private:
    EventAction* fEventAction;

    // ═══════════════════════════════════════════════════════════════
    // PARAMÈTRES DE DEBUG
    // ═══════════════════════════════════════════════════════════════
    G4bool fVerbose;           // Affichage détaillé
    G4int fVerboseMaxEvents;   // Nombre max d'événements verbeux
};

#endif
