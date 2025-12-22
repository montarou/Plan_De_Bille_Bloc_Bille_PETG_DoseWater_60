#include "PrimaryGeneratorAction.hh"

#include "G4LogicalVolumeStore.hh"
#include "G4LogicalVolume.hh"
#include "G4Box.hh"
#include "G4RunManager.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4SystemOfUnits.hh"
#include "Randomize.hh"
#include <cmath>

PrimaryGeneratorAction::PrimaryGeneratorAction()
: G4VUserPrimaryGeneratorAction(),
  fParticleGun(nullptr),
  fLastEventGammaCount(0),
  fConeAngle(60.0 * deg),
  fSourcePosition(0.*mm, 0.*mm, 2.*cm)
{
    G4int n_particle = 1;
    fParticleGun = new G4ParticleGun(n_particle);
    
    // ═══════════════════════════════════════════════════════════════
    // DÉFINITION DU TYPE DE PARTICULE (GAMMA)
    // ═══════════════════════════════════════════════════════════════
    G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
    G4ParticleDefinition* particle = particleTable->FindParticle("gamma");
    fParticleGun->SetParticleDefinition(particle);

    // ═══════════════════════════════════════════════════════════════
    // SPECTRE GAMMA DE L'EUROPIUM-152
    // ═══════════════════════════════════════════════════════════════
    // Raies gamma de l'Europium-152 : énergies (keV) et intensités relatives (%)
    // INTENSITÉ = probabilité d'émission PAR DÉSINTÉGRATION
    // Somme > 100% car plusieurs gammas peuvent être émis par désintégration

    fGammaEnergies = {
        40.12, 39.52, 121.78, 244.70, 344.28, 411.12,
        443.96, 778.90, 867.38, 964.08, 1112.07, 1408.01
    };

    fGammaIntensities = {
        37.7, 20.8, 28.5, 7.6, 26.5, 2.2,
        2.8, 12.9, 4.2, 14.6, 13.6, 21.0
    };

    // ═══════════════════════════════════════════════════════════════
    // CONVERSION EN PROBABILITÉS (fraction entre 0 et 1)
    // ═══════════════════════════════════════════════════════════════
    fGammaProbabilities.clear();
    G4double totalIntensity = 0.0;
    for (size_t i = 0; i < fGammaIntensities.size(); ++i) {
        G4double prob = fGammaIntensities[i] / 100.0;
        fGammaProbabilities.push_back(prob);
        totalIntensity += fGammaIntensities[i];
    }

    // ═══════════════════════════════════════════════════════════════
    // AFFICHAGE DU SPECTRE
    // ═══════════════════════════════════════════════════════════════
    G4cout << "\n╔═══════════════════════════════════════════===════════╗\n";
    G4cout << "║        SPECTRE GAMMA EUROPIUM-152 (Eu-152)             ║\n";
    G4cout << "╠══════════════════════════════════════════════===═==════╣\n";
    G4cout << "║  Énergie (keV)  │  Intensité (%)  │Prob.d'émission     ║\n";
    G4cout << "╠═════════════════╪═════════════════╪════════════======══╣\n";

    for (size_t i = 0; i < fGammaEnergies.size(); ++i) {
        char buffer[120];
        sprintf(buffer, "║    %7.2f      │      %5.1f      │         %6.4f            ║",
                fGammaEnergies[i], fGammaIntensities[i], fGammaProbabilities[i]);
        G4cout << buffer << G4endl;
    }

    G4cout << "╠════════════════════════════════════════════════════════════════=============╣" << G4endl;
    G4cout << "║  Intensité totale : " << totalIntensity << " %                              ║" << G4endl;
    G4cout << "║  Nombre moyen de gammas/désintégration : " << totalIntensity/100.0 << "     ║" << G4endl;
    G4cout << "╠══════════════════════════════════════════════════════════════=============══╣" << G4endl;
    G4cout << "║  MÉTHODE : Pour chaque raie, tirage indépendant avec sa                     ║" << G4endl;
    G4cout << "║            probabilité propre. Permet 0, 1, 2, ... gammas                   ║" << G4endl;
    G4cout << "║            par désintégration (événement).                                  ║" << G4endl;
    G4cout << "╚═════════════════════════════════════════════════════════════=============═══╝\n" << G4endl;

    // ═══════════════════════════════════════════════════════════════
    // POSITION DE LA SOURCE
    // ═══════════════════════════════════════════════════════════════
    // Position de départ : décalée pour viser le centre d'une bille centrale
    // x=1mm, y=1mm aligne avec une des 4 billes centrales du plan

    fParticleGun->SetParticlePosition(fSourcePosition);
}

PrimaryGeneratorAction::~PrimaryGeneratorAction()
{
    delete fParticleGun;
}

void PrimaryGeneratorAction::GenerateDirectionInCone(G4double coneAngle,
                                                     G4double& theta,
                                                     G4double& phi,
                                                     G4ThreeVector& direction)
{
    // ═══════════════════════════════════════════════════════════════
    // GÉNÉRATION D'UNE DIRECTION ALÉATOIRE DANS UN CÔNE
    // ═══════════════════════════════════════════════════════════════
    // Distribution uniforme en angle solide

    G4double cosTheta = 1.0 - G4UniformRand() * (1.0 - std::cos(coneAngle));
    theta = std::acos(cosTheta);
    phi = G4UniformRand() * 2.0 * M_PI;

    G4double sinTheta = std::sin(theta);
    G4double u_x = sinTheta * std::cos(phi);
    G4double u_y = sinTheta * std::sin(phi);
    G4double u_z = cosTheta;

    direction = G4ThreeVector(u_x, u_y, u_z);
}

void PrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{

    // ═══════════════════════════════════════════════════════════════
    // SIMULATION RÉALISTE D'UNE DÉSINTÉGRATION Eu-152
    // ═══════════════════════════════════════════════════════════════
    //
    // Pour chaque raie gamma/X du spectre :
    //   - Tirer un nombre aléatoire r ∈ [0, 1]
    //   - Si r < probabilité de la raie → émettre un gamma de cette énergie
    //
    // Cela permet d'avoir 0, 1, 2, 3... gammas par désintégration,
    // avec une moyenne de 1.924 gammas/désintégration.
    //
    // Chaque gamma généré crée un PrimaryVertex dans G4Event.
    // Ces informations sont récupérées dans EventAction::BeginOfEventAction()
    //
    // ═══════════════════════════════════════════════════════════════

    fLastEventGammaCount = 0;

    // Parcourir toutes les raies
    for (size_t i = 0; i < fGammaEnergies.size(); ++i) {

        // Tirage aléatoire pour cette raie
        G4double random = G4UniformRand();

        // Si le tirage est inférieur à la probabilité, on émet ce gamma
        if (random < fGammaProbabilities[i]) {

            // Énergie de cette raie
            G4double energy = fGammaEnergies[i] * keV;
            fParticleGun->SetParticleEnergy(energy);

            // Direction aléatoire dans le cône (indépendante pour chaque gamma)
            G4double theta, phi;
            G4ThreeVector direction;
            GenerateDirectionInCone(fConeAngle, theta, phi, direction);
            fParticleGun->SetParticleMomentumDirection(direction);

            // ═══════════════════════════════════════════════════════
            // GÉNÉRATION DU VERTEX PRIMAIRE
            // ═══════════════════════════════════════════════════════
            // Cette information est automatiquement stockée dans G4Event
            // et sera accessible via anEvent->GetPrimaryVertex(i)
            // Le trackID sera assigné par Geant4 : 1, 2, 3, ...
            fParticleGun->GeneratePrimaryVertex(anEvent);
            fLastEventGammaCount++;
        }
    }

    // ═══════════════════════════════════════════════════════════════
    // AFFICHAGE DE DIAGNOSTIC
    // ═══════════════════════════════════════════════════════════════

    G4int eventID = anEvent->GetEventID();
    if (eventID < 10 || eventID % 10000 == 0) {
        G4cout << "PrimaryGenerator | Event " << eventID
        << " : " << fLastEventGammaCount << " gamma(s) generated";

        if (fLastEventGammaCount > 0) {
            // Lister les énergies générées
            G4cout << " [";
            for (G4int v = 0; v < anEvent->GetNumberOfPrimaryVertex(); ++v) {
                G4PrimaryVertex* vertex = anEvent->GetPrimaryVertex(v);
                if (vertex && vertex->GetPrimary()) {
                    if (v > 0) G4cout << ", ";
                    G4cout << vertex->GetPrimary()->GetKineticEnergy()/keV << " keV";
                }
            }
            G4cout << "]";
        }
        G4cout << G4endl;
    }
}
