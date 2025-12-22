#include "ActionInitialization.hh"
#include "PrimaryGeneratorAction.hh"
#include "RunAction.hh"
#include "EventAction.hh"
#include "SteppingAction.hh"

ActionInitialization::ActionInitialization()
: G4VUserActionInitialization()
{}

ActionInitialization::~ActionInitialization()
{}

void ActionInitialization::Build() const
{
    // Set primary generator action
    SetUserAction(new PrimaryGeneratorAction);

    // ═══════════════════════════════════════════════════════════════
    // Set run action - STOCKER le pointeur pour le passer à EventAction
    // ═══════════════════════════════════════════════════════════════
    RunAction* runAction = new RunAction();
    SetUserAction(runAction);

    // ═══════════════════════════════════════════════════════════════
    // Set event action - MODIFIÉ : Passer runAction au constructeur
    // ═══════════════════════════════════════════════════════════════
    EventAction* eventAction = new EventAction(runAction);
    SetUserAction(eventAction);

    // Set stepping action (needs EventAction pointer)
    SetUserAction(new SteppingAction(eventAction));

}
