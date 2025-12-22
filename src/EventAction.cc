#include "EventAction.hh"
#include "RunAction.hh"

#include "G4Event.hh"
#include "G4PrimaryVertex.hh"
#include "G4PrimaryParticle.hh"
#include "G4SystemOfUnits.hh"
#include "G4AnalysisManager.hh"
#include <cmath>

EventAction::EventAction(RunAction* runAction)
: G4UserEventAction(),
  fRunAction(runAction),
  fKermaEnergyDeposit(0.),
  fKermaFluence(0.),
  fKermaForcedDeposit(0.),
  fTransmissionTolerance(1.0 * keV),
  fVerboseLevel(1)
{}

EventAction::~EventAction()
{}

void EventAction::BeginOfEventAction(const G4Event* event)
{

    // ═══════════════════════════════════════════════════════════════
    // RESET DE TOUTES LES STRUCTURES AU DÉBUT DE CHAQUE ÉVÉNEMENT
    // ═══════════════════════════════════════════════════════════════
    fPrimaryGammas.clear();
    fSecondariesDownstream.clear();
    fTrackIDtoIndex.clear();
    fKermaEnergyDeposit = 0.;
    fKermaFluence = 0.;
    fKermaForcedDeposit = 0.;

    // ═══════════════════════════════════════════════════════════════
    // RÉCUPÉRATION DES INFORMATIONS DES GAMMAS PRIMAIRES
    // DIRECTEMENT DEPUIS G4Event
    // ═══════════════════════════════════════════════════════════════

    G4int nVertices = event->GetNumberOfPrimaryVertex();
    G4int eventID = event->GetEventID();

    G4int expectedTrackID = 1;

    for (G4int iVertex = 0; iVertex < nVertices; ++iVertex) {

        G4PrimaryVertex* vertex = event->GetPrimaryVertex(iVertex);
        if (!vertex) continue;

        G4PrimaryParticle* primary = vertex->GetPrimary();

        while (primary) {

            G4double energy = primary->GetKineticEnergy();
            G4ThreeVector momentum = primary->GetMomentumDirection();

            G4double theta = std::acos(momentum.z());
            G4double phi = std::atan2(momentum.y(), momentum.x());

            PrimaryGammaInfo info;
            info.trackID = expectedTrackID;
            info.energyInitial = energy;
            info.energyUpstream = 0.;
            info.energyDownstream = 0.;
            info.theta = theta;
            info.phi = phi;
            info.detectedUpstream = false;
            info.detectedDownstream = false;
            info.transmitted = false;

            fTrackIDtoIndex[expectedTrackID] = fPrimaryGammas.size();
            fPrimaryGammas.push_back(info);

            if (fVerboseLevel >= 2 && eventID < 5) {
                G4cout << "  BeginOfEvent | Registered primary gamma: "
                << "trackID=" << expectedTrackID
                << ", E=" << energy/keV << " keV"
                << ", theta=" << theta/deg << " deg"
                << ", phi=" << phi/deg << " deg"
                << G4endl;
            }

            expectedTrackID++;
            primary = primary->GetNext();
        }
    }

    if (fVerboseLevel >= 1 && (eventID < 10 || eventID % 10000 == 0)) {
        G4cout << "BeginOfEvent " << eventID
        << " | " << fPrimaryGammas.size() << " primary gamma(s) registered"
        << G4endl;
    }
}

void EventAction::EndOfEventAction(const G4Event* event)
{

    G4int eventID = event->GetEventID();
    G4int nPrimaries = fPrimaryGammas.size();

    std::vector<G4double> primaryEnergies;
    G4double totalEnergy = 0.;
    G4int nTransmitted = 0;
    G4int nAbsorbed = 0;
    G4int nScattered = 0;

    auto analysisManager = G4AnalysisManager::Instance();

    for (size_t i = 0; i < fPrimaryGammas.size(); ++i) {
        const auto& g = fPrimaryGammas[i];

        primaryEnergies.push_back(g.energyInitial);
        totalEnergy += g.energyInitial;

        if (g.transmitted) {
            nTransmitted++;
        } else if (g.detectedUpstream && g.detectedDownstream) {
            nScattered++;
        } else if (g.detectedUpstream && !g.detectedDownstream) {
            nAbsorbed++;
        }

        // Remplir le ntuple GammaData
        analysisManager->FillNtupleIColumn(1, 0, eventID);
        analysisManager->FillNtupleIColumn(1, 1, i);
        analysisManager->FillNtupleDColumn(1, 2, g.energyInitial/keV);
        analysisManager->FillNtupleDColumn(1, 3, g.energyUpstream/keV);
        analysisManager->FillNtupleDColumn(1, 4, g.energyDownstream/keV);
        analysisManager->FillNtupleDColumn(1, 5, g.theta/deg);
        analysisManager->FillNtupleDColumn(1, 6, g.phi/deg);
        analysisManager->FillNtupleIColumn(1, 7, g.detectedUpstream ? 1 : 0);
        analysisManager->FillNtupleIColumn(1, 8, g.detectedDownstream ? 1 : 0);
        analysisManager->FillNtupleIColumn(1, 9, g.transmitted ? 1 : 0);
        analysisManager->AddNtupleRow(1);
    }

    // Remplir le ntuple EventData
    analysisManager->FillNtupleIColumn(0, 0, eventID);
    analysisManager->FillNtupleIColumn(0, 1, nPrimaries);
    analysisManager->FillNtupleDColumn(0, 2, totalEnergy/keV);
    analysisManager->FillNtupleIColumn(0, 3, nTransmitted);
    analysisManager->FillNtupleIColumn(0, 4, nAbsorbed);
    analysisManager->FillNtupleIColumn(0, 5, nScattered);
    analysisManager->FillNtupleIColumn(0, 6, (G4int)fSecondariesDownstream.size());
    analysisManager->FillNtupleDColumn(0, 7, fKermaEnergyDeposit/keV);
    analysisManager->AddNtupleRow(0);

    // Envoyer les statistiques à RunAction
    if (fRunAction) {
        fRunAction->RecordEventStatistics(nPrimaries, primaryEnergies,
                                          nTransmitted, nAbsorbed,
                                          fKermaEnergyDeposit);
        // ═══════════════════════════════════════════════════════════════
        // MÉTHODE 1 : Dépôt d'énergie Monte Carlo direct
        // ═══════════════════════════════════════════════════════════════
        if (fKermaEnergyDeposit > 0.) {
            fRunAction->AddKermaEnergy(fKermaEnergyDeposit);
        }
        // MÉTHODE 2 : Fluence × μ_en/ρ
        if (fKermaFluence > 0.) {
            fRunAction->AddKermaFluence(fKermaFluence);
        }
        // MÉTHODE 1bis : Forçage d'interaction
        if (fKermaForcedDeposit > 0.) {
            fRunAction->AddKermaEnergyForced(fKermaForcedDeposit);
        }
    }

    // Affichage de diagnostic
    if (fVerboseLevel >= 1 && (eventID < 10 || eventID % 10000 == 0)) {
        G4cout << "\n══════════════════════════════════════════════════" << G4endl;
        G4cout << "EVENT " << eventID << " SUMMARY" << G4endl;
        G4cout << "══════════════════════════════════════════════════" << G4endl;
        G4cout << "Primary gammas: " << nPrimaries << " | Total E: " << totalEnergy/keV << " keV" << G4endl;

        for (size_t i = 0; i < fPrimaryGammas.size(); ++i) {
            const auto& g = fPrimaryGammas[i];
            G4String status = "UNKNOWN";
            if (g.transmitted) status = "TRANSMITTED";
            else if (g.detectedUpstream && g.detectedDownstream) status = "SCATTERED";
            else if (g.detectedUpstream && !g.detectedDownstream) status = "ABSORBED";
            else if (!g.detectedUpstream) status = "MISSED_UPSTREAM";

            G4cout << "  [" << i << "] trackID=" << g.trackID
            << " E_init=" << g.energyInitial/keV << " keV"
            << " → Up(" << g.energyUpstream/keV << ")"
            << " → Down(" << g.energyDownstream/keV << ")"
            << " [" << status << "]" << G4endl;
        }

        if (fSecondariesDownstream.size() > 0) {
            G4cout << "Secondaries downstream: " << fSecondariesDownstream.size() << G4endl;
            for (const auto& s : fSecondariesDownstream) {
                G4cout << "  - PDG=" << s.pdgCode
                << " E=" << s.energy/keV << " keV"
                << " parent=" << s.parentID
                << " process=" << s.creatorProcess << G4endl;
            }
        }

        G4cout << "Dose deposit: " << fKermaEnergyDeposit/keV << " keV" << G4endl;
        G4cout << "Statistics: Transmitted=" << nTransmitted
        << " Absorbed=" << nAbsorbed
        << " Scattered=" << nScattered << G4endl;
        G4cout << "══════════════════════════════════════════════════\n" << G4endl;
    }
}

void EventAction::RecordPrimaryUpstream(G4int trackID, G4double energy)
{
    auto it = fTrackIDtoIndex.find(trackID);
    if (it != fTrackIDtoIndex.end()) {
        size_t index = it->second;
        fPrimaryGammas[index].energyUpstream = energy;
        fPrimaryGammas[index].detectedUpstream = true;
    }
}

void EventAction::RecordPrimaryDownstream(G4int trackID, G4double energy)
{
    auto it = fTrackIDtoIndex.find(trackID);
    if (it != fTrackIDtoIndex.end()) {
        size_t index = it->second;
        fPrimaryGammas[index].energyDownstream = energy;
        fPrimaryGammas[index].detectedDownstream = true;

        G4double deltaE = std::abs(fPrimaryGammas[index].energyUpstream - energy);
        if (deltaE < fTransmissionTolerance) {
            fPrimaryGammas[index].transmitted = true;
        }
    }
}

void EventAction::RecordSecondaryDownstream(G4int trackID, G4int parentID,
                                            G4int pdgCode, G4double energy,
                                            const G4String& process)
{
    SecondaryParticleInfo info;
    info.trackID = trackID;
    info.parentID = parentID;
    info.pdgCode = pdgCode;
    info.energy = energy;
    info.creatorProcess = process;

    fSecondariesDownstream.push_back(info);
}

void EventAction::AddKermaEnergy(G4double edep)
{
    fKermaEnergyDeposit += edep;
}

void EventAction::AddKermaFluence(G4double energy, G4double chordLength)
{
    // ═══════════════════════════════════════════════════════════════
    // CALCUL DE LA DOSE PAR FLUENCE (Méthode 2)
    // Utilise les coefficients μ_en/ρ tabulés pour l'EAU (NIST XCOM)
    // Formule : E_déposé = E × L × μ_en/ρ × ρ_eau (identique à Méthode 1bis)
    // ═══════════════════════════════════════════════════════════════

    // Table des coefficients μ_en/ρ pour l'EAU (cm²/g) - NIST XCOM
    static const G4double energyTable[] = {
        10., 15., 20., 30., 40., 50., 60., 80., 100.,
        150., 200., 300., 400., 500., 600., 800., 1000., 1500., 2000.
    };
    static const G4double muEnRhoTable[] = {
        4.944, 1.374, 0.5503, 0.1557, 0.0688, 0.0420, 0.0318, 0.0257, 0.0253,
        0.0278, 0.0299, 0.0319, 0.0328, 0.0330, 0.0329, 0.0321, 0.0311, 0.0283, 0.0260
    };
    static const G4int nPoints = 19;

    G4double E_keV = energy / keV;
    G4double L_cm = chordLength / cm;
    G4double mu_en_rho = 0.;

    // Interpolation log-log pour trouver μ_en/ρ à l'énergie donnée
    if (E_keV <= energyTable[0]) {
        mu_en_rho = muEnRhoTable[0];
    } else if (E_keV >= energyTable[nPoints-1]) {
        mu_en_rho = muEnRhoTable[nPoints-1];
    } else {
        for (G4int i = 0; i < nPoints - 1; i++) {
            if (E_keV >= energyTable[i] && E_keV < energyTable[i+1]) {
                G4double logE = std::log(E_keV);
                G4double logE1 = std::log(energyTable[i]);
                G4double logE2 = std::log(energyTable[i+1]);
                G4double logMu1 = std::log(muEnRhoTable[i]);
                G4double logMu2 = std::log(muEnRhoTable[i+1]);

                G4double logMu = logMu1 + (logMu2 - logMu1) * (logE - logE1) / (logE2 - logE1);
                mu_en_rho = std::exp(logMu);
                break;
            }
        }
    }

    // Densité de l'eau
    G4double rho_water = 1.0;  // g/cm³

    // Contribution à la dose : E × L × μ_en/ρ × ρ (converti en unités Geant4)
    fKermaFluence += E_keV * L_cm * mu_en_rho * rho_water * keV;
}

void EventAction::AddKermaEnergyForced(G4double energy, G4double chordLength)
{
    // ═══════════════════════════════════════════════════════════════
    // CALCUL DE LA DOSE PAR FORÇAGE D'INTERACTION (Méthode 1bis)
    // Pour chaque gamma traversant le détecteur, on calcule l'énergie
    // qu'il aurait déposée en moyenne : E_dep = E × L × μ_en × ρ_eau
    // ═══════════════════════════════════════════════════════════════

    // Table des coefficients μ_en/ρ pour l'EAU (cm²/g) - NIST XCOM
    static const G4double energyTable[] = {
        10., 15., 20., 30., 40., 50., 60., 80., 100.,
        150., 200., 300., 400., 500., 600., 800., 1000., 1500., 2000.
    };
    static const G4double muEnRhoTable[] = {
        4.944, 1.374, 0.5503, 0.1557, 0.0688, 0.0420, 0.0318, 0.0257, 0.0253,
        0.0278, 0.0299, 0.0319, 0.0328, 0.0330, 0.0329, 0.0321, 0.0311, 0.0283, 0.0260
    };
    static const G4int nPoints = 19;

    G4double E_keV = energy / keV;
    G4double L_cm = chordLength / cm;
    G4double mu_en_rho = 0.;

    // Interpolation log-log
    if (E_keV <= energyTable[0]) {
        mu_en_rho = muEnRhoTable[0];
    } else if (E_keV >= energyTable[nPoints-1]) {
        mu_en_rho = muEnRhoTable[nPoints-1];
    } else {
        for (G4int i = 0; i < nPoints - 1; i++) {
            if (E_keV >= energyTable[i] && E_keV < energyTable[i+1]) {
                G4double logE = std::log(E_keV);
                G4double logE1 = std::log(energyTable[i]);
                G4double logE2 = std::log(energyTable[i+1]);
                G4double logMu1 = std::log(muEnRhoTable[i]);
                G4double logMu2 = std::log(muEnRhoTable[i+1]);

                G4double logMu = logMu1 + (logMu2 - logMu1) * (logE - logE1) / (logE2 - logE1);
                mu_en_rho = std::exp(logMu);
                break;
            }
        }
    }

    // Densité de l'eau
    G4double rho_water = 1.0;  // g/cm³

    // Énergie forcée déposée : E × L × μ_en/ρ × ρ
    G4double forcedDeposit_keV = E_keV * L_cm * mu_en_rho * rho_water;

    fKermaForcedDeposit += forcedDeposit_keV * keV;
}

G4int EventAction::GetNumberTransmitted() const
{
    G4int count = 0;
    for (const auto& gamma : fPrimaryGammas) {
        if (gamma.transmitted) count++;
    }
    return count;
}

G4int EventAction::GetNumberAbsorbed() const
{
    G4int count = 0;
    for (const auto& gamma : fPrimaryGammas) {
        if (gamma.detectedUpstream && !gamma.detectedDownstream) count++;
    }
    return count;
}

G4bool EventAction::IsPrimaryTrack(G4int trackID) const
{
    return fTrackIDtoIndex.find(trackID) != fTrackIDtoIndex.end();
}
