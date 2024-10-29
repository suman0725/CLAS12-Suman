#include <cstdlib>
#include <iostream>
#include <chrono>
#include <TFile.h>
#include <TTree.h>
#include <TApplication.h>
#include <TROOT.h>
#include <TDatabasePDG.h>
#include <TLorentzVector.h>
#include <TH1.h>
#include <TH2.h>
#include <TChain.h>
#include <TCanvas.h>
#include <TBenchmark.h>
#include "clas12reader.h"
#include "HipoChain.h"
#include <cmath>
#include <filesystem>


namespace fs = std::filesystem;
using namespace clas12;

void SetLorentzVector(TLorentzVector &p4, clas12::region_part_ptr rp) {
    p4.SetXYZM(rp->par()->getPx(), rp->par()->getPy(), rp->par()->getPz(), p4.M());
}

// Function to calculate Pt^2 for the produced hadron
double CalculatePt2(TLorentzVector P_h, TLorentzVector q) {
    // Get the angle between P_h and q
    double theta_a = P_h.Vect().Angle(q.Vect());

    // Calculate Pt^2 using |P_h|^2 * sin^2(theta_a)
    double Pt2 = P_h.P() * P_h.P() * std::sin(theta_a) * std::sin(theta_a);
    
    return Pt2;
}




// Function to compute φh
double CalculatePhih(TLorentzVector q, TLorentzVector p, TLorentzVector p_h) {
    // Extract 3-vectors from the Lorentz vectors
    TVector3 q3 = q.Vect();   // q vector
    TVector3 p3 = p.Vect();   // another particle's 
    TVector3 p_h3 = p_h.Vect(); // hadron's vector

    // Compute the cross products
    TVector3 cross_qp = q3.Cross(p3);      // (q x p)
    TVector3 cross_qph = q3.Cross(p_h3);   // (q x p_h)

    // Compute the dot product of the cross products
    double dotProduct = cross_qp.Dot(cross_qph);

    // Compute the magnitudes of the cross products
    double mag_qp = cross_qp.Mag();
    double mag_qph = cross_qph.Mag();

    // Calculate cos(φh)
    double cos_phi_h = dotProduct / (mag_qp * mag_qph);

    // Return φh by applying the inverse cosine (acos)
    return acos(cos_phi_h) * TMath::RadToDeg(); // result in radians
}


void addHipoFiles(clas12root::HipoChain& chain, const fs::path& baseDir) {
    std::vector<fs::path> files;
    for (const auto& entry : fs::recursive_directory_iterator(baseDir)) {
        if (entry.path().extension() == ".hipo") {
            files.push_back(entry.path());
        }
    }
    std::sort(files.begin(), files.end());
    for (const auto& file : files) {
        chain.Add(file.string());
        std::cout << "Added file: " << file << std::endl;
    }
}



void old1() {
    auto db = TDatabasePDG::Instance();
    double beamEnergy = 10.532; // GeV
    TLorentzVector beam(0, 0, beamEnergy, beamEnergy); // Assuming the beam is along the z-axis
    TLorentzVector target(0, 0, 0, db->GetParticle(2212)->Mass()); // Proton target
    TLorentzVector el, pip;

    // Histograms for kinematic variables
    auto* hnu = new TH1F("h_nu", "LD2_nu;nu(GeV);Counts", 200, 1, 11);
    auto* hQ2 = new TH1F("Q2", "Q^{2} distribution;Q^{2} (GeV)^{2};Counts", 200, 0, 8);
    auto* hxB = new TH1F("xB", "Bjorken x distribution;xB;Counts", 100, 0, 1);
    auto* hW2 = new TH1F("W2", "W^{2} distribution;W^{2} (GeV)^{2};Counts", 200, 0, 20);
    auto* hPt2 = new TH1F("Pt2", "Pt^2 distribution;Pt^2;Counts", 100, 0, 3);
    auto* hz = new TH1F("z", "z distribution;z;Counts", 100, 0, 1); // New histogram for z
     auto* hPhi_h = new TH1F("Phi_h", "Phi_h distribution;Phi_h (Degree);Counts", 100, 0, 180); // New histogram for z

    clas12root::HipoChain chain;
    //string baseDirectory = "/lustre24/expphy/cache/hallb/scratch/rg-d/production/Bspot/v5dstLD2/dst/recon/018309";
    //addHipoFiles(chain, baseDirectory);
    chain.Add("/lustre24/expphy/cache/hallb/scratch/rg-d/production/Bspot/v5dstLD2/dst/recon/018309/rec_clas_018309.evio.00035-00039.hipo");
    chain.db()->turnOffQADB();

    auto& c12 = chain.C12ref();
    while (chain.Next()) {
        auto particles = c12->getDetParticles();
        TLorentzVector el_temp, pip_temp;

        bool foundElectron = false;
        bool foundPion = false;

        double nu = 0; // Declare nu outside of the loop
        TLorentzVector q; // Declare q outside of the loop

        for (auto& p : particles) {
            int pid = p->par()->getPid();
            int status = p->par()->getStatus();
            int chi2Pid = p->par()->getChi2Pid();

            // Electron selection
            if (pid == 11 && status < 0 && chi2Pid < 5 && chi2Pid > -5) {
                SetLorentzVector(el_temp, p);
                foundElectron = true;
                nu = beam.Energy() - el_temp.Energy(); // Calculate nu
                q = beam - el_temp; // Calculate q after finding the electron
                double Q2 = -q.Mag2();
                double W2 = target.M() * target.M() - Q2 + 2 * target.M() * nu;
                double xB = Q2 / (2 * target.M() * nu);

                hnu->Fill(nu);
                hQ2->Fill(Q2);
                hxB->Fill(xB);
                hW2->Fill(W2);
            }

            // Pion selection
            if (pid == 211 && chi2Pid > -10 && chi2Pid < 10) {
                SetLorentzVector(pip_temp, p);
                foundPion = true;
            }
        }

        // Only calculate and fill histograms if both electron and pion are found
        if (foundElectron && foundPion) {
            // Calculate Pt^2 using the selected pion and q

            
            double Pt2 = CalculatePt2(pip_temp, q); 
            hPt2->Fill(Pt2);
            std::cout << "Calculated Pt^2: " << Pt2 << std::endl;

            // Calculate z using the energy of the pion and nu
            double z = pip_temp.E() / nu; // z = E_h / nu
            hz->Fill(z); // Fill the z histogram

            // Calculate Phi_h using the beam, target, pip_temp
            double Phi_ppi =  CalculatePhih(q,beam,pip_temp);
            hPhi_h->Fill(Phi_ppi);
        }
    }

    // Drawing histograms and saving as PDF
    TCanvas* can1 = new TCanvas();
    hnu->Draw();
    can1->Print("Kinematic_Variables_LD2_018309.pdf(", "pdf");

    TCanvas* can2 = new TCanvas();
    hQ2->Draw();
    can2->Print("Kinematic_Variables_LD2_018309.pdf", "pdf");

    TCanvas* can3 = new TCanvas();
    hxB->Draw();
    can3->Print("Kinematic_Variables_LD2_018309.pdf", "pdf");

    TCanvas* can4 = new TCanvas();
    hW2->Draw();
    can4->Print("Kinematic_Variables_LD2_018309.pdf", "pdf");

    // Old Pt2 histogram (log scale and x-axis range)
    TCanvas* can5 = new TCanvas();
    can5->SetLogy();            // Set y-axis to log scale
    hPt2->GetXaxis()->SetRangeUser(0, 3);  // Set x-axis range
    hPt2->Draw();
    can5->Print("Kinematic_Variables_LD2_018309.pdf", "pdf");

    TCanvas* can6 = new TCanvas();
    hPhi_h->Draw();
    can6->Print("Kinematic_Variables_LD2_018309.pdf", "pdf");

    // New z histogram (log scale and x-axis range)
    TCanvas* can7 = new TCanvas();
    can7->SetLogy();            // Set y-axis to log scale
    hz->GetXaxis()->SetRangeUser(0, 1); // Set x-axis range
    hz->Draw();
    can7->Print("Kinematic_Variables_LD2_018309.pdf)", "pdf");
}
