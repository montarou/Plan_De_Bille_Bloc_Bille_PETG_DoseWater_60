#include "RunAction.hh"

#include "G4Run.hh"
#include "G4RunManager.hh"
#include "G4AnalysisManager.hh"
#include "G4SystemOfUnits.hh"
#include <cmath>

RunAction::RunAction()
: G4UserRunAction(),
  fKermaTotalEnergy(0.),
  fKermaTotalEnergy2(0.),
  fKermaMass(0.),
  fKermaRadius(2.0*cm),
  fKermaPosition(20.0*cm),
  fKermaEventCount(0),
  fKermaFluenceTotal(0.),
  fKermaFluenceCount(0),
  fKermaForcedTotal(0.),
  fKermaForcedCount(0),
  fActivity4pi(44000.0),
  fConeAngle(60.0*deg),
  fSourcePosZ(2.0*cm),            // Position Z de la source
  fDetectorPosZ(20.0*cm),         // Position Z du détecteur
  fMeanGammasPerDecay(1.924),     // Nombre moyen de gammas par désintégration
  fTotalPrimariesGenerated(0),
  fTotalEventsWithZeroGamma(0),
  fTotalTransmitted(0),
  fTotalAbsorbed(0),
  fOutputFileName("dose_water_output")
{

    // ═══════════════════════════════════════════════════════════════
    // CONFIGURATION DE G4AnalysisManager
    // ═══════════════════════════════════════════════════════════════
    auto analysisManager = G4AnalysisManager::Instance();
    analysisManager->SetDefaultFileType("root");
    analysisManager->SetVerboseLevel(1);
    analysisManager->SetNtupleMerging(true);

    // ═══════════════════════════════════════════════════════════════
    // CALCUL DU FACTEUR DE CORRECTION GÉOMÉTRIQUE
    // ═══════════════════════════════════════════════════════════════
    G4double f_cone = (1.0 - std::cos(fConeAngle)) / 2.0;  // = 0.25 pour 60°

    G4cout << "\n╔════════════════════════════════════════════════════════════╗" << G4endl;
    G4cout << "║  RunAction initialized with " << analysisManager->GetType() << G4endl;
    G4cout << "║  Output file: " << fOutputFileName << ".root" << G4endl;
    G4cout << "║  Detector material: WATER (tissue equivalent)" << G4endl;
    G4cout << "╠════════════════════════════════════════════════════════════╣" << G4endl;
    G4cout << "║  *** NORMALISATION ***                                     ║" << G4endl;
    G4cout << "║  Temps simulé : t = N_events / A_4π                        ║" << G4endl;
    G4cout << "║  Facteur correction : f_cone = " << f_cone << " (cône " << fConeAngle/deg << "°)" << G4endl;
    G4cout << "║  Dose corrigée = Dose_brute × f_cone                       ║" << G4endl;
    G4cout << "╚════════════════════════════════════════════════════════════╝\n" << G4endl;
}

RunAction::~RunAction()
{}

void RunAction::BeginOfRunAction(const G4Run* run)
{
    G4cout << "### Run " << run->GetRunID() << " start." << G4endl;

    // ═══════════════════════════════════════════════════════════════
    // RESET DES COMPTEURS
    // ═══════════════════════════════════════════════════════════════
    fKermaTotalEnergy = 0.;
    fKermaTotalEnergy2 = 0.;
    fKermaEventCount = 0;
    fKermaFluenceTotal = 0.;
    fKermaFluenceCount = 0;
    fKermaForcedTotal = 0.;
    fKermaForcedCount = 0;
    fTotalPrimariesGenerated = 0;
    fTotalEventsWithZeroGamma = 0;
    fTotalTransmitted = 0;
    fTotalAbsorbed = 0;

    // ═══════════════════════════════════════════════════════════════
    // CALCUL DE LA MASSE DU DÉTECTEUR D'EAU
    // ═══════════════════════════════════════════════════════════════
    G4double kerma_volume = (4.0/3.0) * M_PI * std::pow(fKermaRadius, 3);
    G4double water_density = 1.0 * g/cm3;  // Eau
    fKermaMass = kerma_volume * water_density;

    // ═══════════════════════════════════════════════════════════════
    // CRÉATION DES HISTOGRAMMES ET NTUPLES
    // ═══════════════════════════════════════════════════════════════
    auto analysisManager = G4AnalysisManager::Instance();

    analysisManager->OpenFile(fOutputFileName);

    // Histogramme 0 : Nombre de gammas primaires par événement
    analysisManager->CreateH1("nGammasPerEvent",
                              "Number of primary gammas per event;N_{#gamma};Counts",
                              15, -0.5, 14.5);

    // Histogramme 1 : Spectre des énergies générées
    analysisManager->CreateH1("energySpectrum",
                              "Energy spectrum of generated gammas;E (keV);Counts",
                              1500, 0., 1500.);

    // Histogramme 2 : Énergie totale par événement
    analysisManager->CreateH1("totalEnergyPerEvent",
                              "Total primary energy per event;E_{tot} (keV);Counts",
                              500, 0., 5000.);

    // Histogramme 3 : Nombre de gammas transmis par événement
    analysisManager->CreateH1("nTransmittedPerEvent",
                              "Number of transmitted gammas per event;N_{trans};Counts",
                              15, -0.5, 14.5);

    // Histogramme 4 : Nombre de gammas absorbés par événement
    analysisManager->CreateH1("nAbsorbedPerEvent",
                              "Number of absorbed gammas per event;N_{abs};Counts",
                              15, -0.5, 14.5);

    // Histogramme 5 : Dépôt d'énergie dose par événement
    analysisManager->CreateH1("dosePerEvent",
                              "Dose energy deposit per event;E_{dose} (keV);Counts",
                              200, 0., 100.);

    // Histogramme 2D : Corrélation nombre de gammas vs énergie totale
    analysisManager->CreateH2("nGammas_vs_totalEnergy",
                              "Number of gammas vs Total energy;N_{#gamma};E_{tot} (keV)",
                              15, -0.5, 14.5,
                              100, 0., 5000.);

    // Ntuple 0 : Informations par événement
    analysisManager->CreateNtuple("EventData", "Event-level data");
    analysisManager->CreateNtupleIColumn("eventID");
    analysisManager->CreateNtupleIColumn("nPrimaries");
    analysisManager->CreateNtupleDColumn("totalEnergy");
    analysisManager->CreateNtupleIColumn("nTransmitted");
    analysisManager->CreateNtupleIColumn("nAbsorbed");
    analysisManager->CreateNtupleIColumn("nScattered");
    analysisManager->CreateNtupleIColumn("nSecondaries");
    analysisManager->CreateNtupleDColumn("doseDeposit");
    analysisManager->FinishNtuple();

    // Ntuple 1 : Informations par gamma primaire
    analysisManager->CreateNtuple("GammaData", "Primary gamma data");
    analysisManager->CreateNtupleIColumn("eventID");
    analysisManager->CreateNtupleIColumn("gammaIndex");
    analysisManager->CreateNtupleDColumn("energyInitial");
    analysisManager->CreateNtupleDColumn("energyUpstream");
    analysisManager->CreateNtupleDColumn("energyDownstream");
    analysisManager->CreateNtupleDColumn("theta");
    analysisManager->CreateNtupleDColumn("phi");
    analysisManager->CreateNtupleIColumn("detectedUpstream");
    analysisManager->CreateNtupleIColumn("detectedDownstream");
    analysisManager->CreateNtupleIColumn("transmitted");
    analysisManager->FinishNtuple();

    G4cout << "\n╔════════════════════════════════════════════════════════════╗" << G4endl;
    G4cout << "║  HISTOGRAMMES ET NTUPLES CRÉÉS                             ║" << G4endl;
    G4cout << "╠════════════════════════════════════════════════════════════╣" << G4endl;
    G4cout << "║  H1[0] : nGammasPerEvent      - Distribution N_gamma       ║" << G4endl;
    G4cout << "║  H1[1] : energySpectrum       - Spectre des énergies       ║" << G4endl;
    G4cout << "║  H1[2] : totalEnergyPerEvent  - Énergie totale/event       ║" << G4endl;
    G4cout << "║  H1[3] : nTransmittedPerEvent - N transmis/event           ║" << G4endl;
    G4cout << "║  H1[4] : nAbsorbedPerEvent    - N absorbés/event           ║" << G4endl;
    G4cout << "║  H1[5] : dosePerEvent         - Dépôt dose/event           ║" << G4endl;
    G4cout << "║  H2[0] : nGammas_vs_totalEnergy - Corrélation 2D           ║" << G4endl;
    G4cout << "╠════════════════════════════════════════════════════════════╣" << G4endl;
    G4cout << "║  Ntuple[0] : EventData  - 1 entrée par événement           ║" << G4endl;
    G4cout << "║  Ntuple[1] : GammaData  - 1 entrée par gamma primaire      ║" << G4endl;
    G4cout << "╚════════════════════════════════════════════════════════════╝\n" << G4endl;

}

void RunAction::RecordEventStatistics(G4int nPrimaries,
                                      const std::vector<G4double>& primaryEnergies,
                                      G4int nTransmitted,
                                      G4int nAbsorbed,
                                      G4double doseDeposit)
{
    auto analysisManager = G4AnalysisManager::Instance();

    // H0 : Nombre de gammas par événement
    analysisManager->FillH1(0, nPrimaries);

    // H1 : Spectre des énergies
    G4double totalEnergy = 0.;
    for (const auto& energy : primaryEnergies) {
        analysisManager->FillH1(1, energy/keV);
        totalEnergy += energy;
    }

    // H2 : Énergie totale par événement
    analysisManager->FillH1(2, totalEnergy/keV);

    // H3 : Nombre de transmis par événement
    analysisManager->FillH1(3, nTransmitted);

    // H4 : Nombre d'absorbés par événement
    analysisManager->FillH1(4, nAbsorbed);

    // H5 : Dépôt dose par événement
    analysisManager->FillH1(5, doseDeposit/keV);

    // H2D : Corrélation
    analysisManager->FillH2(0, nPrimaries, totalEnergy/keV);

    // Mise à jour des compteurs globaux
    fTotalPrimariesGenerated += nPrimaries;
    if (nPrimaries == 0) fTotalEventsWithZeroGamma++;
    fTotalTransmitted += nTransmitted;
    fTotalAbsorbed += nAbsorbed;
}

void RunAction::AddKermaEnergy(G4double edep)
{
    fKermaTotalEnergy += edep;
    fKermaTotalEnergy2 += edep * edep;
    if (edep > 0) fKermaEventCount++;
}

void RunAction::AddKermaFluence(G4double fluence)
{
    fKermaFluenceTotal += fluence;
    fKermaFluenceCount++;
}

void RunAction::AddKermaEnergyForced(G4double forcedDeposit)
{
    fKermaForcedTotal += forcedDeposit;
    fKermaForcedCount++;
}

void RunAction::EndOfRunAction(const G4Run* run)
{
    G4int nofEvents = run->GetNumberOfEvent();
    if (nofEvents == 0) {
        G4cout << "\n### Run ended: No events processed.\n" << G4endl;
        return;
    }

    // ═══════════════════════════════════════════════════════════════
    // FERMETURE DU FICHIER D'ANALYSE
    // ═══════════════════════════════════════════════════════════════
    auto analysisManager = G4AnalysisManager::Instance();
    analysisManager->Write();
    analysisManager->CloseFile();

    // ═══════════════════════════════════════════════════════════════
    // CALCUL DES STATISTIQUES FINALES
    // ═══════════════════════════════════════════════════════════════
    G4double meanGammasPerEvent = (G4double)fTotalPrimariesGenerated / nofEvents;
    G4double fractionZeroGamma = (G4double)fTotalEventsWithZeroGamma / nofEvents * 100.;

    G4double transmissionRate = 0.;
    G4double absorptionRate = 0.;
    if (fTotalPrimariesGenerated > 0) {
        transmissionRate = (G4double)fTotalTransmitted / fTotalPrimariesGenerated * 100.;
        absorptionRate = (G4double)fTotalAbsorbed / fTotalPrimariesGenerated * 100.;
    }

    // Calcul de la dose (méthode 1)
    G4double meanDosePerEvent = fKermaTotalEnergy / nofEvents;
    G4double doseGy = (fKermaTotalEnergy / fKermaMass) / gray;

    // Variance et erreur
    G4double variance = 0.;
    if (nofEvents > 1) {
        variance = (fKermaTotalEnergy2 / nofEvents) -
                   std::pow(fKermaTotalEnergy / nofEvents, 2);
        variance = std::max(variance, 0.);
    }
    G4double stdDev = std::sqrt(variance);
    G4double stdError = stdDev / std::sqrt(nofEvents);

    // ═══════════════════════════════════════════════════════════════
    // NORMALISATION PAR DÉSINTÉGRATIONS (4π) + CORRECTION GÉOMÉTRIQUE
    // ═══════════════════════════════════════════════════════════════
    // 
    // Principe :
    // 1. Temps simulé = N_événements / A_4π (indépendant du cône)
    // 2. Facteur de correction f_cone = Ω_cône / 4π
    //    - Pour un cône de 60° : f_cone = (1 - cos(60°))/2 = 0.25
    // 3. Dose réelle = Dose_brute × f_cone
    //
    // Cette méthode est équivalente à normaliser par A_cône directement,
    // mais permet de comparer facilement différentes configurations.
    // ═══════════════════════════════════════════════════════════════

    // Temps simulé basé sur l'activité 4π
    G4double simulatedTime_4pi_s = (G4double)nofEvents / fActivity4pi;
    G4double simulatedTime_4pi_h = simulatedTime_4pi_s / 3600.0;

    // Facteur de correction géométrique (fraction d'angle solide du cône)
    G4double f_cone = (1.0 - std::cos(fConeAngle)) / 2.0;  // = 0.25 pour 60°

    // Paramètres géométriques du détecteur
    G4double distance = fDetectorPosZ - fSourcePosZ;
    G4double detectorTheta = std::atan(fKermaRadius / distance);
    G4double detectorSolidAngleFrac = (1.0 - std::cos(detectorTheta)) / 2.0;
    G4double activityToDetector = fActivity4pi * detectorSolidAngleFrac;
    G4double detectorArea_cm2 = M_PI * std::pow(fKermaRadius/cm, 2);

    // Nombre de gammas attendus dans le détecteur (vérification géométrique)
    G4double expectedGammasInDetector = (G4double)fTotalPrimariesGenerated * 
                                        (detectorSolidAngleFrac / f_cone);

    // ═══════════════════════════════════════════════════════════════
    // MÉTHODE 1 : DÉPÔT D'ÉNERGIE (Monte Carlo)
    // ═══════════════════════════════════════════════════════════════

    G4double doseTotal_Gy = (fKermaTotalEnergy / keV) * 1.602e-16 / (fKermaMass / kg);

    // Débit brut (sans correction)
    G4double doseRate_brut_GyPerS = doseTotal_Gy / simulatedTime_4pi_s;
    G4double doseRate_brut_nGyPerH = doseRate_brut_GyPerS * 3600.0 * 1.0e9;

    // Débit corrigé (× f_cone)
    G4double doseRate_corrige_nGyPerH = doseRate_brut_nGyPerH * f_cone;

    G4double convergence = (meanDosePerEvent > 0) ? stdError / meanDosePerEvent * 100. : 0.;
    G4double doseRateError_corrige_nGyPerH = doseRate_corrige_nGyPerH * convergence / 100.;

    // ═══════════════════════════════════════════════════════════════
    // MÉTHODE 1bis : FORÇAGE D'INTERACTION
    // ═══════════════════════════════════════════════════════════════

    G4double doseForced_Gy = (fKermaForcedTotal / fKermaMass) / gray;

    // Débit brut
    G4double doseRateForced_brut_GyPerS = doseForced_Gy / simulatedTime_4pi_s;
    G4double doseRateForced_brut_nGyPerH = doseRateForced_brut_GyPerS * 3600.0 * 1.0e9;

    // Débit corrigé
    G4double doseRateForced_corrige_nGyPerH = doseRateForced_brut_nGyPerH * f_cone;

    G4double convergenceForced = (fKermaForcedCount > 0) ? 
                                  100.0 / std::sqrt((G4double)fKermaForcedCount) : 0.;
    G4double doseRateForcedError_corrige_nGyPerH = doseRateForced_corrige_nGyPerH * convergenceForced / 100.;

    // ═══════════════════════════════════════════════════════════════
    // MÉTHODE 2 : FLUENCE × μ_en/ρ
    // ═══════════════════════════════════════════════════════════════

    G4double doseFluence_Gy = (fKermaFluenceTotal / fKermaMass) / gray;

    // Débit brut
    G4double doseRateFluence_brut_GyPerS = doseFluence_Gy / simulatedTime_4pi_s;
    G4double doseRateFluence_brut_nGyPerH = doseRateFluence_brut_GyPerS * 3600.0 * 1.0e9;

    // Débit corrigé
    G4double doseRateFluence_corrige_nGyPerH = doseRateFluence_brut_nGyPerH * f_cone;

    G4double convergenceFluence = (fKermaFluenceCount > 0) ?
                                   100.0 / std::sqrt((G4double)fKermaFluenceCount) : 0.;
    G4double doseRateFluenceError_corrige_nGyPerH = doseRateFluence_corrige_nGyPerH * convergenceFluence / 100.;

    // ═══════════════════════════════════════════════════════════════
    // MÉTRIQUES DE COMPARAISON (indépendantes de la normalisation)
    // ═══════════════════════════════════════════════════════════════

    G4double energyPerGamma_MC = (fTotalPrimariesGenerated > 0) ?
        (fKermaTotalEnergy/keV) / fTotalPrimariesGenerated : 0.;
    
    G4double energyPerGamma_Forced = (fTotalPrimariesGenerated > 0) ?
        (fKermaForcedTotal/keV) / fTotalPrimariesGenerated : 0.;

    G4double energyPerGamma_Fluence = (fTotalPrimariesGenerated > 0) ?
        (fKermaFluenceTotal/keV) / fTotalPrimariesGenerated : 0.;

    // ═══════════════════════════════════════════════════════════════
    // VALEUR THÉORIQUE DE RÉFÉRENCE (Eu-152, air, 18 cm)
    // ═══════════════════════════════════════════════════════════════
    G4double doseRate_theorique = 174.8;  // nGy/h (calculé analytiquement)

    // Écarts par rapport à la théorie
    G4double ecart_MC = (doseRate_theorique > 0) ? 
        100.0 * std::abs(doseRate_corrige_nGyPerH - doseRate_theorique) / doseRate_theorique : 0.;
    G4double ecart_Forced = (doseRate_theorique > 0) ? 
        100.0 * std::abs(doseRateForced_corrige_nGyPerH - doseRate_theorique) / doseRate_theorique : 0.;

    // ═══════════════════════════════════════════════════════════════
    // AFFICHAGE DU RÉSUMÉ
    // ═══════════════════════════════════════════════════════════════
    G4cout << "\n";
    G4cout << "╔═══════════════════════════════════════════════════════════════════════════╗\n";
    G4cout << "║                         RUN SUMMARY                                       ║\n";
    G4cout << "║          NORMALISATION 4π + CORRECTION GÉOMÉTRIQUE f_cone                ║\n";
    G4cout << "╠═══════════════════════════════════════════════════════════════════════════╣\n";
    G4cout << "║  Number of events processed: " << nofEvents << G4endl;
    G4cout << "╠═══════════════════════════════════════════════════════════════════════════╣\n";
    G4cout << "║  PRIMARY GAMMA GENERATION STATISTICS:                                     ║\n";
    G4cout << "║    Total gammas generated     : " << fTotalPrimariesGenerated << G4endl;
    G4cout << "║    Mean gammas per event      : " << meanGammasPerEvent << G4endl;
    G4cout << "║    Expected (theory)          : 1.924" << G4endl;
    G4cout << "║    Events with 0 gamma        : " << fTotalEventsWithZeroGamma
           << " (" << fractionZeroGamma << "%)" << G4endl;
    G4cout << "║    Expected P(N=0)            : ~11.7%" << G4endl;
    G4cout << "╠═══════════════════════════════════════════════════════════════════════════╣\n";
    G4cout << "║  TRANSMISSION STATISTICS:                                                 ║\n";
    G4cout << "║    Gammas transmitted         : " << fTotalTransmitted
           << " (" << transmissionRate << "%)" << G4endl;
    G4cout << "║    Gammas absorbed            : " << fTotalAbsorbed
           << " (" << absorptionRate << "%)" << G4endl;
    G4cout << "╠═══════════════════════════════════════════════════════════════════════════╣\n";
    G4cout << "║  SOURCE PARAMETERS:                                                       ║\n";
    G4cout << "║    Activity (4π)              : " << fActivity4pi << " Bq" << G4endl;
    G4cout << "║    Cone half-angle            : " << fConeAngle/deg << " deg" << G4endl;
    G4cout << "║    Cone solid angle fraction  : " << f_cone*100 << " %" << G4endl;
    G4cout << "╠═══════════════════════════════════════════════════════════════════════════╣\n";
    G4cout << "║  *** NORMALISATION ***                                                    ║\n";
    G4cout << "║    Temps simulé (4π)          : " << simulatedTime_4pi_s << " s ("
           << simulatedTime_4pi_h << " h)" << G4endl;
    G4cout << "║    Formule                    : t = N_events / A_4π" << G4endl;
    G4cout << "╠═══════════════════════════════════════════════════════════════════════════╣\n";
    G4cout << "║  *** CORRECTION GÉOMÉTRIQUE ***                                           ║\n";
    G4cout << "║    Facteur f_cone             : " << f_cone << " (Ω_cône/4π)" << G4endl;
    G4cout << "║    Formule                    : Dose_réelle = Dose_brute × f_cone" << G4endl;
    G4cout << "╠═══════════════════════════════════════════════════════════════════════════╣\n";
    G4cout << "║  DETECTOR GEOMETRY:                                                       ║\n";
    G4cout << "║    Position (z)               : " << fDetectorPosZ/cm << " cm" << G4endl;
    G4cout << "║    Distance source-detecteur  : " << distance/cm << " cm" << G4endl;
    G4cout << "║    Detector radius            : " << fKermaRadius/cm << " cm" << G4endl;
    G4cout << "║    Detector area (πr²)        : " << detectorArea_cm2 << " cm²" << G4endl;
    G4cout << "║    Detector mass              : " << fKermaMass/g << " g" << G4endl;
    G4cout << "║    Detector half-angle        : " << detectorTheta/deg << " deg" << G4endl;
    G4cout << "║    Detector solid angle frac  : " << detectorSolidAngleFrac*100 << " %" << G4endl;
    G4cout << "╠═══════════════════════════════════════════════════════════════════════════╣\n";
    G4cout << "║  GAMMAS DANS LE DÉTECTEUR:                                                ║\n";
    G4cout << "║    Gammas entrants (observé)  : " << fKermaForcedCount << G4endl;
    G4cout << "║    Gammas attendus (géom.)    : " << (G4int)expectedGammasInDetector << G4endl;
    G4cout << "║    Fraction du cône           : " << 100.0*fKermaForcedCount/fTotalPrimariesGenerated << " %" << G4endl;
    G4cout << "╠═══════════════════════════════════════════════════════════════════════════╣\n";
    G4cout << "║  *** MÉTHODE 1 : DÉPÔT D'ÉNERGIE (Monte Carlo) ***                        ║\n";
    G4cout << "║    Total energy deposited     : " << fKermaTotalEnergy/keV << " keV" << G4endl;
    G4cout << "║    Events with deposit        : " << fKermaEventCount << G4endl;
    G4cout << "║    Energy per gamma generated : " << energyPerGamma_MC << " keV/gamma" << G4endl;
    G4cout << "║    Débit BRUT                 : " << doseRate_brut_nGyPerH << " nGy/h" << G4endl;
    G4cout << "║  ┌─────────────────────────────────────────────────────────────────────┐  ║\n";
    G4cout << "║  │ Débit CORRIGÉ (× " << f_cone << ")     : " << doseRate_corrige_nGyPerH 
           << " ± " << doseRateError_corrige_nGyPerH << " nGy/h" << G4endl;
    G4cout << "║  └─────────────────────────────────────────────────────────────────────┘  ║\n";
    G4cout << "╠═══════════════════════════════════════════════════════════════════════════╣\n";
    G4cout << "║  *** MÉTHODE 1bis : FORÇAGE D'INTERACTION ***                             ║\n";
    G4cout << "║    Gammas avec forçage        : " << fKermaForcedCount << G4endl;
    G4cout << "║    Énergie forcée totale      : " << fKermaForcedTotal/keV << " keV" << G4endl;
    G4cout << "║    Energy per gamma generated : " << energyPerGamma_Forced << " keV/gamma" << G4endl;
    G4cout << "║    Débit BRUT                 : " << doseRateForced_brut_nGyPerH << " nGy/h" << G4endl;
    G4cout << "║  ┌─────────────────────────────────────────────────────────────────────┐  ║\n";
    G4cout << "║  │ Débit CORRIGÉ (× " << f_cone << ")     : " << doseRateForced_corrige_nGyPerH 
           << " ± " << doseRateForcedError_corrige_nGyPerH << " nGy/h" << G4endl;
    G4cout << "║  └─────────────────────────────────────────────────────────────────────┘  ║\n";
    G4cout << "╠═══════════════════════════════════════════════════════════════════════════╣\n";
    G4cout << "║  *** MÉTHODE 2 : FLUENCE × μ_en/ρ ***                                     ║\n";
    G4cout << "║    Gammas traversant          : " << fKermaFluenceCount << G4endl;
    G4cout << "║    Σ(E × L × μ_en)            : " << fKermaFluenceTotal/keV << " keV" << G4endl;
    G4cout << "║    Energy per gamma generated : " << energyPerGamma_Fluence << " keV/gamma" << G4endl;
    G4cout << "║    Débit BRUT                 : " << doseRateFluence_brut_nGyPerH << " nGy/h" << G4endl;
    G4cout << "║  ┌─────────────────────────────────────────────────────────────────────┐  ║\n";
    G4cout << "║  │ Débit CORRIGÉ (× " << f_cone << ")     : " << doseRateFluence_corrige_nGyPerH 
           << " ± " << doseRateFluenceError_corrige_nGyPerH << " nGy/h" << G4endl;
    G4cout << "║  └─────────────────────────────────────────────────────────────────────┘  ║\n";
    G4cout << "╠═══════════════════════════════════════════════════════════════════════════╣\n";
    G4cout << "║  *** COMPARAISON AVEC VALEUR THÉORIQUE ***                                ║\n";
    G4cout << "║    Valeur théorique (Eu-152)  : " << doseRate_theorique << " nGy/h" << G4endl;
    G4cout << "║    Écart Méthode 1            : " << ecart_MC << " %" << G4endl;
    G4cout << "║    Écart Méthode 1bis         : " << ecart_Forced << " %" << G4endl;
    if (ecart_Forced < 5.0) {
        G4cout << "║    → EXCELLENT ! Écart < 5%                                              ║\n";
    } else if (ecart_Forced < 10.0) {
        G4cout << "║    → BON. Écart < 10%                                                    ║\n";
    } else {
        G4cout << "║    → ATTENTION : Écart > 10%                                             ║\n";
    }
    G4cout << "╠═══════════════════════════════════════════════════════════════════════════╣\n";
    G4cout << "║  *** MÉTRIQUES DE COMPARAISON (indépendantes du temps) ***                ║\n";
    G4cout << "║    E_déposée/gamma (MC)       : " << energyPerGamma_MC << " keV/gamma" << G4endl;
    G4cout << "║    E_forcée/gamma             : " << energyPerGamma_Forced << " keV/gamma" << G4endl;
    G4cout << "║    E_fluence/gamma            : " << energyPerGamma_Fluence << " keV/gamma" << G4endl;
    G4cout << "╠═══════════════════════════════════════════════════════════════════════════╣\n";
    G4cout << "║  OUTPUT FILE: " << fOutputFileName << ".root" << G4endl;
    G4cout << "╚═══════════════════════════════════════════════════════════════════════════╝\n";
    G4cout << G4endl;
}
