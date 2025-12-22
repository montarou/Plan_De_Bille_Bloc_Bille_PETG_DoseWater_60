#include "SteppingAction.hh"
#include "EventAction.hh"

#include "G4Step.hh"
#include "G4Track.hh"
#include "G4StepPoint.hh"
#include "G4VPhysicalVolume.hh"
#include "G4SystemOfUnits.hh"
#include "G4RunManager.hh"
#include <cmath>

SteppingAction::SteppingAction(EventAction* eventAction)
: G4UserSteppingAction(),
  fEventAction(eventAction),
  fVerbose(false),
  fVerboseMaxEvents(5)
{}

SteppingAction::~SteppingAction()
{}

void SteppingAction::UserSteppingAction(const G4Step* step)
{
    // ═══════════════════════════════════════════════════════════════
    // RÉCUPÉRATION DES INFORMATIONS DE BASE
    // ═══════════════════════════════════════════════════════════════

    G4StepPoint* preStepPoint = step->GetPreStepPoint();
    G4StepPoint* postStepPoint = step->GetPostStepPoint();

    // Vérifier que les volumes existent
    if (!preStepPoint->GetPhysicalVolume() || !postStepPoint->GetPhysicalVolume()) {
        return;
    }

    G4String preVolumeName = preStepPoint->GetPhysicalVolume()->GetName();
    G4String postVolumeName = postStepPoint->GetPhysicalVolume()->GetName();

    // Informations sur la trace
    G4Track* track = step->GetTrack();
    G4int trackID = track->GetTrackID();
    G4int parentID = track->GetParentID();
    G4String particleName = track->GetDefinition()->GetParticleName();
    G4double kineticEnergy = preStepPoint->GetKineticEnergy();

    // ID de l'événement pour le debug
    G4int eventID = G4RunManager::GetRunManager()->GetCurrentEvent()->GetEventID();

    // ═══════════════════════════════════════════════════════════════
    // DÉTECTION DANS LE DÉTECTEUR DOSE (Méthode 1 : dépôt d'énergie)
    // ═══════════════════════════════════════════════════════════════

    G4LogicalVolume* volume = preStepPoint->GetTouchableHandle()
                                        ->GetVolume()
                                        ->GetLogicalVolume();
    G4String logicalVolumeName = volume->GetName();

    // CHANGEMENT : DoseDetectorLog au lieu de KermaDetectorLog
    if (logicalVolumeName == "DoseDetectorLog") {
        G4double edep = step->GetTotalEnergyDeposit();
        if (edep > 0.) {
            fEventAction->AddKermaEnergy(edep);
        }
    }

    // ═══════════════════════════════════════════════════════════════
    // DÉTECTION DES GAMMAS ENTRANT DANS LE DÉTECTEUR DOSE (Méthode 2 : fluence)
    // ═══════════════════════════════════════════════════════════════

    G4String preLogVolName = "OutOfWorld";
    G4String postLogVolName = "OutOfWorld";

    if (preStepPoint->GetPhysicalVolume()) {
        preLogVolName = preStepPoint->GetPhysicalVolume()->GetLogicalVolume()->GetName();
    }
    if (postStepPoint->GetPhysicalVolume()) {
        postLogVolName = postStepPoint->GetPhysicalVolume()->GetLogicalVolume()->GetName();
    }

    // CHANGEMENT : DoseDetectorLog au lieu de KermaDetectorLog
    if (postLogVolName == "DoseDetectorLog" && preLogVolName != "DoseDetectorLog") {
        if (particleName == "gamma") {
            G4double gammaEnergy = postStepPoint->GetKineticEnergy();

            G4ThreeVector position = postStepPoint->GetPosition();
            G4ThreeVector direction = track->GetMomentumDirection();

            // ───────────────────────────────────────────────────────
            // CALCUL DE LA LONGUEUR DE CORDE DANS LA SPHÈRE
            // ───────────────────────────────────────────────────────
            G4ThreeVector center(0., 0., 20.*cm);
            G4double radius = 2.0*cm;

            G4ThreeVector oc = position - center;
            G4double b = 2.0 * oc.dot(direction);
            G4double c = oc.mag2() - radius*radius;
            G4double discriminant = b*b - 4.0*c;

            G4double chordLength = 0.;
            if (discriminant >= 0.) {
                chordLength = std::sqrt(discriminant);
            }

            // ───────────────────────────────────────────────────────
            // MÉTHODE 2 : FLUENCE × μ_en/ρ × L (longueur de corde)
            // ───────────────────────────────────────────────────────
            if (chordLength > 0.) {
                fEventAction->AddKermaFluence(gammaEnergy, chordLength);
            }

            // ───────────────────────────────────────────────────────
            // MÉTHODE 1bis : FORÇAGE D'INTERACTION
            // E_déposé = E × L × μ_en × ρ_eau
            // ───────────────────────────────────────────────────────
            if (chordLength > 0.) {
                fEventAction->AddKermaEnergyForced(gammaEnergy, chordLength);
            }
        }
    }

    // ═══════════════════════════════════════════════════════════════
    // DÉTECTION AU PLAN UPSTREAM
    // ═══════════════════════════════════════════════════════════════

    if (postVolumeName == "UpstreamDetector" && preVolumeName != "UpstreamDetector") {

        G4ThreeVector momentum = track->GetMomentumDirection();
        G4double pz = momentum.z();

        if (pz > 0) {

            if (parentID == 0 && particleName == "gamma") {
                fEventAction->RecordPrimaryUpstream(trackID, kineticEnergy);

                if (fVerbose && eventID < fVerboseMaxEvents) {
                    G4cout << "  UPSTREAM | Event " << eventID
                    << " | PRIMARY gamma trackID=" << trackID
                    << " E=" << kineticEnergy/keV << " keV"
                    << G4endl;
                }
            }
        }
    }

    // ═══════════════════════════════════════════════════════════════
    // DÉTECTION AU PLAN DOWNSTREAM
    // ═══════════════════════════════════════════════════════════════

    if (postVolumeName == "DownstreamDetector" && preVolumeName != "DownstreamDetector") {

        G4ThreeVector momentum = track->GetMomentumDirection();
        G4double pz = momentum.z();

        if (pz > 0) {

            if (parentID == 0 && particleName == "gamma") {
                fEventAction->RecordPrimaryDownstream(trackID, kineticEnergy);

                if (fVerbose && eventID < fVerboseMaxEvents) {
                    G4cout << "  DOWNSTREAM | Event " << eventID
                    << " | PRIMARY gamma trackID=" << trackID
                    << " E=" << kineticEnergy/keV << " keV"
                    << G4endl;
                }

            } else {
                G4String processName = "Unknown";
                const G4VProcess* creatorProcess = track->GetCreatorProcess();
                if (creatorProcess) {
                    processName = creatorProcess->GetProcessName();
                }

                G4int pdgCode = track->GetDefinition()->GetPDGEncoding();

                fEventAction->RecordSecondaryDownstream(
                    trackID,
                    parentID,
                    pdgCode,
                    kineticEnergy,
                    processName
                );

                if (fVerbose && eventID < fVerboseMaxEvents) {
                    G4cout << "  DOWNSTREAM | Event " << eventID
                    << " | SECONDARY " << particleName
                    << " trackID=" << trackID
                    << " parentID=" << parentID
                    << " E=" << kineticEnergy/keV << " keV"
                    << " process=" << processName
                    << G4endl;
                }
            }
        }
    }

    // ═══════════════════════════════════════════════════════════════
    // DEBUG DÉTAILLÉ (optionnel)
    // ═══════════════════════════════════════════════════════════════

    if (fVerbose && eventID < fVerboseMaxEvents) {
        G4ThreeVector prePos = preStepPoint->GetPosition();
        G4ThreeVector postPos = postStepPoint->GetPosition();

        if (parentID == 0 && particleName == "gamma") {
            G4cout << "  Step | trackID=" << trackID
            << " " << preVolumeName << " → " << postVolumeName
            << " z: " << prePos.z()/mm << " → " << postPos.z()/mm << " mm"
            << " E=" << kineticEnergy/keV << " keV"
            << G4endl;
        }
    }
}
