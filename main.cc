//
// ********************************************************************
// * Simulation Geant4 - Dose dans l'eau                             *
// * Source Eu-152, détecteur sphérique d'eau à z = 20 cm            *
// * MODE SÉQUENTIEL (mono-thread)                                   *
// ********************************************************************
//

#include "DetectorConstruction.hh"
#include "ActionInitialization.hh"
#include "PhysicsList.hh"

#include "G4RunManager.hh"
#include "G4UImanager.hh"
#include "G4VisExecutive.hh"
#include "G4UIExecutive.hh"

#include <fstream>
#include <iostream>

int main(int argc, char** argv)
{
    // ═══════════════════════════════════════════════════════════════
    // REDIRECTION DE G4cout/G4cerr VERS FICHIER output.log
    // ═══════════════════════════════════════════════════════════════
    std::ofstream logFile("output.log");
    std::streambuf* coutBuffer = std::cout.rdbuf();
    std::streambuf* cerrBuffer = std::cerr.rdbuf();
    std::cout.rdbuf(logFile.rdbuf());
    std::cerr.rdbuf(logFile.rdbuf());
    
    G4cout << "\n";
    G4cout << "╔═══════════════════════════════════════════════════════════════╗\n";
    G4cout << "║         SIMULATION GEANT4 - DOSE DANS L'EAU                   ║\n";
    G4cout << "║                                                               ║\n";
    G4cout << "║  Source     : Eu-152 (44 kBq)                                 ║\n";
    G4cout << "║  Détecteur  : Sphère d'EAU (r = 2.0 cm) à z = 20 cm           ║\n";
    G4cout << "║  Objectif   : Mesure de la dose absorbée (tissu mou)          ║\n";
    G4cout << "╚═══════════════════════════════════════════════════════════════╝\n";
    G4cout << G4endl;

    // Détection du mode interactif
    G4UIExecutive* ui = nullptr;
    if (argc == 1) {
        ui = new G4UIExecutive(argc, argv);
    }

    // ═══════════════════════════════════════════════════════════════
    // Création du run manager en mode SÉQUENTIEL (mono-thread)
    // ═══════════════════════════════════════════════════════════════
    G4RunManager* runManager = new G4RunManager;
    
    // Géométrie obligatoire
    runManager->SetUserInitialization(new DetectorConstruction());
    
    // Liste de physique personnalisée (FTFP_BERT + Livermore + StepLimiter)
    runManager->SetUserInitialization(new PhysicsList());
    
    // Initialisation des actions utilisateur
    runManager->SetUserInitialization(new ActionInitialization());

    // Initialisation du manager de visualisation
    G4VisManager* visManager = new G4VisExecutive;
    visManager->Initialize();

    // Pointeur vers le UI manager
    G4UImanager* UImanager = G4UImanager::GetUIpointer();

    if (!ui) {
        // Mode batch
        G4String command = "/control/execute ";
        G4String fileName = argv[1];
        UImanager->ApplyCommand(command + fileName);
    }
    else {
        // Mode interactif
        UImanager->ApplyCommand("/control/execute init_vis.mac");
        ui->SessionStart();
        delete ui;
    }

    // Nettoyage
    delete visManager;
    delete runManager;
    
    // ═══════════════════════════════════════════════════════════════
    // RESTAURER LES BUFFERS DE SORTIE ORIGINAUX
    // ═══════════════════════════════════════════════════════════════
    std::cout.rdbuf(coutBuffer);
    std::cerr.rdbuf(cerrBuffer);
    logFile.close();
    
    return 0;
}
