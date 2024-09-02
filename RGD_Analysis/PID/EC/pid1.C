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
#include <TLine.h>

using namespace clas12;
using namespace std;
int getSector(double phi) {
    // Assuming phi is already in the expected range of [-180, 180]
    if (phi >= -30 && phi < 30) return 1;
    if (phi >= 30 && phi < 90) return 2;
    if (phi >= 90 && phi < 150) return 3;
    if ((phi >= 150 && phi <= 180) || (phi >= -180 && phi < -150)) return 4;
    if (phi >= -150 && phi < -90) return 5;
    if (phi >= -90 && phi < -30) return 6;

    return -1; // Default if not in any sector
}

void pid1() {     

// Process LD Target Files
    clas12root::HipoChain chainLD;
    chainLD.Add("/lustre24/expphy/cache/hallb/scratch/rg-d/production/prod/v4ob_aideLD2/dst/recon/018528/rec_clas_018528.evio.00310-00314.hipo");
    chainLD.db()->turnOffQADB();
    auto& c12LD = chainLD.C12ref();

    int status; 
    int sector; 
    
     TH1F* h_test = new TH1F("h_test", "Test Histogram", 100, 0, 0.2);
     TH1F* h_SimpleMomentum = new TH1F("h_SimpleMomentum", "Momentum", 100, 0, 10);

    int event_count = 0;
    std::vector<TH1F*> h_nPhe_negative(6);
    std::vector<TH2F*> h_SamplingFraction_vs_Momentum(6);
    // Create 2D histograms for ECin + ECout vs Epcal for each sector
    std::vector<TH2F*> h_Ecin_Ecout_vs_Epcal(6);
    for (int i = 0; i < 6; ++i) {
        h_Ecin_Ecout_vs_Epcal[i] = new TH2F(Form("h_Ecin_Ecout_vs_Epcal_sector%d", i+1), 
                                            Form("ECin + ECout vs Epcal for Sector %d", i+1), 
                                            200, 0, 1.8,  // X-axis range for ECin + ECout
                                            200, 0, 2);  // Y-axis range for Epcal
    }
    
    for (int i=0; i<6; ++i){
     h_SamplingFraction_vs_Momentum[i] = new TH2F(Form("h_SamplingFraction_vs_Momentum_sector%d", i+1),
                                                     Form("Sampling Fraction vs Momentum for Sector %d", i+1),
                                                     200, 0, 13,  // X-axis range for Momentum
                                                     200, 0, 0.3); // Y-axis range for Sampling Fraction
    }

    while (chainLD.Next()) {        
    auto particles = c12LD->getDetParticles();
        for (auto& p : particles) { 
            // Check if the particle is negatively charged
            if (p->par()->getPid() == 11) { 
		int status = p->par()->getStatus(); 

		//Apply forward detector status cut only 
		if (status > -4000 && status <= 2000) {
		double phi = p->getPhi()* TMath::RadToDeg();

		double ECin = p->cal(ECIN)->getEnergy();
		double ECout = p->cal(ECOUT)->getEnergy();
		double Epcal = p->cal(PCAL)->getEnergy();
		double momentum = p->par()->getP(); 
		double samplingFraction = (ECin + ECout) / momentum; // Calculate sampling fraction   
		sector = getSector(phi);


		h_test->Fill(samplingFraction);

		 h_SimpleMomentum->Fill(momentum); 
		
//		cout << "Phi: " << phi << " Sector: " << sector << endl;

			if (sector > 0 && sector <= 6) {
				// Fill 2D histogram for ECin + ECout vs Epcal
	                        h_Ecin_Ecout_vs_Epcal[sector - 1]->Fill(Epcal,ECin + ECout);
                               
 			        h_SamplingFraction_vs_Momentum[sector - 1]->Fill(momentum,samplingFraction);

	                }
		}
            }
        }

    }
    std::cout << "Histogram entry count for sector " << sector << ": " 
          << h_Ecin_Ecout_vs_Epcal[sector - 1]->GetEntries() << std::endl;
     // Create a canvas to display the histograms
    TCanvas* canvas = new TCanvas("canvas", "ECin + ECout vs Epcal", 1200, 800);
    canvas->Divide(3, 2); // Divide the canvas into a 3x2 grid
    for (int i = 0; i < 6; ++i) {
        canvas->cd(i+1);  // Move to the i-th pad
        h_Ecin_Ecout_vs_Epcal[i]->Draw("COLZ");

     }
     canvas->Update(); 
    // Save the canvas as a PNG file
    canvas->SaveAs("Ecin_Ecout_vs_Epcal_electrons_LD2.png");
    // Create the second canvas for Sampling Fraction vs Momentum
    TCanvas* canvas2 = new TCanvas("canvas2", "Sampling Fraction vs Momentum", 1200, 800);
    canvas2->Divide(3, 2); // Divide canvas into a 3x2 grid
    for (int i = 0; i < 6; ++i) {
        canvas2->cd(i+1);
        h_SamplingFraction_vs_Momentum[i]->Draw("COLZ"); // Draw Sampling Fraction vs Momentum
    }
    canvas2->Update(); 
    // Save the second canvas as a PNG file
    canvas2->SaveAs("SamplingFraction_vs_Momentum_electron_LD2.png");
    TCanvas* c1 = new TCanvas("c1", "Sampling Fraction", 800, 600);
    c1->Divide(2,1);
    c1->cd(1);
    // Draw the 1D histogram on the canvas
    h_test->Draw();

    // Display the canvas
    c1->cd(2); 
    h_SimpleMomentum->Draw(); 
    c1->Update();

}
