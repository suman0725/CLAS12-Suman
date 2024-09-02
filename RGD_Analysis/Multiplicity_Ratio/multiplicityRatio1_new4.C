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
#include <TGraph.h>

// Function to set the Lorentz vector for a given particle
void SetLorentzVector(TLorentzVector &p4, clas12::region_part_ptr rp) {
    p4.SetXYZM(rp->par()->getPx(), rp->par()->getPy(), rp->par()->getPz(), p4.M());
}

void multiplicityRatio1_new4() {
    auto db = TDatabasePDG::Instance();
    const double m_e = db->GetParticle(11)->Mass(); // Mass of electron
    const double m_pion = db->GetParticle(211)->Mass(); // Mass of pion

    int nBins = 7;
    double nuMin = 1.5;
    double nuMax = 8.5;

    // Histograms for LD and Nuclear targets
    auto* h_LD_electrons = new TH1F("h_LD_electrons", "LD Electrons;nu;Counts", nBins, nuMin, nuMax);
    auto* h_LD_hadrons = new TH1F("h_LD_hadrons", "LD Hadrons;nu;Counts", nBins, nuMin, nuMax);
    auto* h_Nuc_electrons = new TH1F("h_Nuc_electrons", "Nuclear Target Electrons;nu;Counts", nBins, nuMin, nuMax);
    auto* h_Nuc_hadrons = new TH1F("h_Nuc_hadrons", "Nuclear Target Hadrons;nu;Counts", nBins, nuMin, nuMax);

    // Histogram for nu 
    auto* h_nu = new TH1F("h_nu", "LD2_nu", 200, 1, 10);
    auto* h_Nuc_nu = new TH1F("h_Nuc_nu", "Carbon_nu", 200, 1, 10);
    // Process LD Target Files
    clas12root::HipoChain chainLD;
    chainLD.Add("/lustre19/expphy/cache/hallb/scratch/rg-d/production/prod/v4ob_aideLD2/dst/recon/018528/rec_clas_018528.evio.00310-00314.hipo");
    chainLD.db()->turnOffQADB();
    auto& c12LD = chainLD.C12ref();

    TLorentzVector beam(0, 0, 10.6, 10.6); // Beam Lorentz vector
    TLorentzVector target(0, 0, 0, db->GetParticle(2212)->Mass()); // Target proton

    while (chainLD.Next()) {
        auto particles = c12LD->getDetParticles();
        TLorentzVector el; // Electron Lorentz vector

        double max_E = 0.0;
        double nu = 0.0;
        bool eventPassesCuts = false;

        // Process all particles
        for (auto& p : particles) {
            int pid = p->par()->getPid();

            if (pid == 11) { // Electron
                SetLorentzVector(el, p);
                double E = el.E();
                if (E > max_E) {
                    max_E = E;
                    nu = beam.E() - el.E();
		    h_nu->Fill(nu); 
                    // Kinematic cuts for electron
                    TLorentzVector q = beam - el;
                    double Q2 = -q.Mag2();
                    double W = sqrt(target.M() * target.M() + 2 * target.M() * nu - Q2);
                    double y = nu / beam.E();

                    if (Q2 > 1 && y > 0.25 && y < 0.85 && W > 2) {
                        eventPassesCuts = true; // Mark event as passing the cuts

                        // Fill electron histogram
                        h_LD_electrons->Fill(nu);
                    }
                }
            } else if (pid == 211) { // Positive pion (hadron)
                // Only fill hadron histograms if the event passes the electron cuts
                if (eventPassesCuts) {
                    TLorentzVector hadron;
                    SetLorentzVector(hadron, p);

                    // Fill histogram for hadrons
                    h_LD_hadrons->Fill(nu);
                }
            }
        }
    }

    // Process Nuclear Target Files
    clas12root::HipoChain chainNuc;
    chainNuc.Add("/lustre19/expphy/cache/hallb/scratch/rg-d/production/prod/v4ob_aideCxC/dst/recon/018451/rec_clas_018451.evio.00310-00314.hipo");
    chainNuc.db()->turnOffQADB();
    auto& c12Nuc = chainNuc.C12ref();

    while (chainNuc.Next()) {
        auto particles = c12Nuc->getDetParticles();
        TLorentzVector el; // Electron Lorentz vector

        double max_E_c = 0.0;
        double nu_c = 0.0;
        bool eventPassesCuts = false;

        // Process all particles
        for (auto& p : particles) {
            int pid = p->par()->getPid();

            if (pid == 11) { // Electron
                SetLorentzVector(el, p);
                double E = el.E();
                if (E > max_E_c) {
                    max_E_c = E;
                    nu_c = beam.E() - el.E();
	  	    h_Nuc_nu->Fill(nu_c); 
                    // Kinematic cuts for electron
                    TLorentzVector q = beam - el;
                    double Q2 = -q.Mag2();
                    double W = sqrt(target.M() * target.M() + 2 * target.M() * nu_c - Q2);
                    double y = nu_c / beam.E();

                    if (Q2 > 1 && y > 0.25 && y < 0.85 && W > 2) {
                        eventPassesCuts = true; // Mark event as passing the cuts

                        // Fill electron histogram
                        h_Nuc_electrons->Fill(nu_c);
                    }
                }
            } else if (pid == 211) { // Positive pion (hadron)
                // Only fill hadron histograms if the event passes the electron cuts
                if (eventPassesCuts) {
                    TLorentzVector hadron;
                    SetLorentzVector(hadron, p);

                    // Fill histogram for hadrons
                    h_Nuc_hadrons->Fill(nu_c);
                }
            }
        }
    }

    // Apply double ratio calculation
    TGraph* graph = new TGraph();
    int pointIndex = 0;

    for (int bin = 1; bin <= nBins; ++bin) {
        double binCenter = h_LD_electrons->GetBinCenter(bin);
        double count_LD_electrons = h_LD_electrons->GetBinContent(bin);
        double count_LD_hadrons = h_LD_hadrons->GetBinContent(bin);
        double count_Nuc_electrons = h_Nuc_electrons->GetBinContent(bin);
        double count_Nuc_hadrons = h_Nuc_hadrons->GetBinContent(bin);

        double ratio = 0.0;
        if (count_LD_electrons > 0 && count_Nuc_electrons > 0) {
            ratio = (count_Nuc_hadrons / count_Nuc_electrons) / (count_LD_hadrons / count_LD_electrons);
        }

        graph->SetPoint(pointIndex++, binCenter, ratio);
    }

    // Save histograms and scatter plot
    TFile* outFile = new TFile("multiplicity_ratio_with_cuts.root", "RECREATE");
    h_LD_electrons->Write();
    h_LD_hadrons->Write();
    h_Nuc_electrons->Write();
    h_Nuc_hadrons->Write();
    graph->Write();

    // Optionally display scatter plot in a canvas
    TCanvas* canvas = new TCanvas("canvas", "Multiplicity Ratio vs nu", 800, 600);
    graph->SetTitle("Multiplicity Ratio vs nu;nu;R");
    graph->SetMarkerStyle(kFullCircle);
    graph->SetMarkerColor(kBlue);
    graph->Draw("AP");
    canvas->SaveAs("multiplicity_ratio_vs_nu_with_cuts.png");

    TCanvas* can1 = new TCanvas();
    h_nu->Draw(); 

    TCanvas* can2 = new TCanvas();
    h_Nuc_nu->Draw(); 



}

