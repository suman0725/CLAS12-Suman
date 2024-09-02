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
#include <TGraph.h>
#include "clas12reader.h"
#include "HipoChain.h"

using namespace clas12;
using namespace std;

void multiplicityRatio2_new() {
    auto db = TDatabasePDG::Instance();
    // Define the number of bins and the range for z
    int nBins = 7;
    double zMin = 0.1;
    double zMax = 0.8;
    const double m_e = db->GetParticle(11)->Mass(); // Mass of electron
    const double m_pion = db->GetParticle(211)->Mass(); // Mass of pion

    // Histograms for LD and Nuclear targets
    auto* h_LD_electrons = new TH1F("h_LD_electrons", "LD Electrons;z;Counts", nBins, zMin, zMax);
    auto* h_LD_hadrons = new TH1F("h_LD_hadrons", "LD Hadrons;z;Counts", nBins, zMin, zMax);
    auto* h_Nuc_electrons = new TH1F("h_Nuc_electrons", "Nuclear Target Electrons;z;Counts", nBins, zMin, zMax);
    auto* h_Nuc_hadrons = new TH1F("h_Nuc_hadrons", "Nuclear Target Hadrons;z;Counts", nBins, zMin, zMax);

    // Process LD Target Files
    clas12root::HipoChain chainLD;
    chainLD.Add("/lustre19/expphy/cache/hallb/scratch/rg-d/production/prod/v4ob_aideLD2/dst/recon/018528/rec_clas_018528.evio.00310-00314.hipo");
    chainLD.db()->turnOffQADB();
    auto& c12LD = chainLD.C12ref();

    while (chainLD.Next()) {
        auto particles = c12LD->getDetParticles();
        double nu = 0.0;

        for (auto& p : particles) {
            int pid = p->par()->getPid();

            if (pid == 11) { // Electron
                auto px = p->par()->getPx();
                auto py = p->par()->getPy();
                auto pz = p->par()->getPz();
                double E = sqrt(px * px + py * py + pz * pz + m_e * m_e);
                nu = 10.532 - E;
                h_LD_electrons->Fill(nu);
            } else if (pid == 211) { // Hadron (e.g., pions)
                auto px = p->par()->getPx();
                auto py = p->par()->getPy();
                auto pz = p->par()->getPz();
                double E_pion = sqrt(px * px + py * py + pz * pz + m_pion * m_pion);
                double z = E_pion / nu;
                h_LD_hadrons->Fill(z);
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
        double nu_c = 0.0;

        for (auto& p : particles) {
            int pid = p->par()->getPid();

            if (pid == 11) { // Electron
                auto px = p->par()->getPx();
                auto py = p->par()->getPy();
                auto pz = p->par()->getPz();
                double E = sqrt(px * px + py * py + pz * pz + m_e * m_e);
                nu_c = 10.532 - E;
                h_Nuc_electrons->Fill(nu_c);
            } else if (abs(pid) == 211) { // Hadron (e.g., pions)
                auto px = p->par()->getPx();
                auto py = p->par()->getPy();
                auto pz = p->par()->getPz();
                double E_pion = sqrt(px * px + py * py + pz * pz + m_pion * m_pion);
                double z = E_pion / nu_c;
                h_Nuc_hadrons->Fill(z);
            }
        }
    }

    // Apply double ratio calculation
    TGraph* graph = new TGraph();
    int pointIndex = 0;

    for (int bin = 1; bin <= nBins; ++bin) {
        double binCenter = h_LD_hadrons->GetBinCenter(bin);
        double count_LD_hadrons = h_LD_hadrons->GetBinContent(bin);
        double count_LD_electrons = h_LD_electrons->GetBinContent(bin);
        double count_Nuc_hadrons = h_Nuc_hadrons->GetBinContent(bin);
        double count_Nuc_electrons = h_Nuc_electrons->GetBinContent(bin);

        double ratio = 0.0;
        if (count_LD_electrons > 0 && count_Nuc_electrons > 0) {
            ratio = (count_Nuc_hadrons / count_Nuc_electrons) / (count_LD_hadrons / count_LD_electrons);
        }

        graph->SetPoint(pointIndex++, binCenter, ratio);
    }

    // Save histograms and scatter plot
    TFile* outFile = new TFile("multiplicity_ratio_z.root", "RECREATE");
    h_LD_electrons->Write();
    h_LD_hadrons->Write();
    h_Nuc_electrons->Write();
    h_Nuc_hadrons->Write();
    graph->Write();
    outFile->Close();

    // Optionally display scatter plot in a canvas
    TCanvas* canvas = new TCanvas("canvas", "Multiplicity Ratio vs z", 800, 600);
    graph->SetTitle("Multiplicity Ratio vs z;z;R");
    graph->SetMarkerStyle(kFullCircle);
    graph->SetMarkerColor(kBlue);
    graph->Draw("AP");
    canvas->SaveAs("multiplicity_ratio_vs_z.png");
}
	
