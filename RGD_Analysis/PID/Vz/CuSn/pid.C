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

void pid() {     

// Process LD Target Files
    clas12root::HipoChain chainLD;
    chainLD.Add("/lustre24/expphy/cache/hallb/scratch/rg-d/production/prod/v4ob_aideCuSn/dst/recon/018624/rec_clas_018624.evio.00310-00314.hipo");
    chainLD.db()->turnOffQADB();
    auto& c12LD = chainLD.C12ref();

    int status; 
    int sector; 
    
    // Define histograms for vertex coordinates vx, vy, vz
TH1F* h_vx[6];
TH1F* h_vy[6];
TH1F* h_vz[6];

// Initialize histograms for each sector
for (int sector = 1; sector <= 6; ++sector) {
    h_vx[sector - 1] = new TH1F(Form("h_vx_sector%d", sector), Form("Vertex vx for sector %d", sector), 100, -10, 10);
    h_vy[sector - 1] = new TH1F(Form("h_vy_sector%d", sector), Form("Vertex vy for sector %d", sector), 100, -10, 10);
    h_vz[sector - 1] = new TH1F(Form("h_vz_sector%d", sector), Form("Vertex vz for sector %d", sector), 100, -20, 20);
}

    while (chainLD.Next()) {        
        auto particles = c12LD->getDetParticles(); 
        for (auto& p : particles) { 
            // Check if the particle is an electron
            if (p->par()->getPid() == 11) { 
		int status = p->par()->getStatus(); 

		//Apply forward detector status cut only 
		if (status > -4000 && status <= 2000) {
		double phi = p->getPhi()* TMath::RadToDeg();
		sector = getSector(phi);
		double vx = p->par()->getVx(); 
		double vy = p->par()->getVy();
		double vz = p->par()->getVz();
		
//		cout << "Phi: " << phi << " Sector: " << sector << endl;

			if (sector > 0 && sector <= 6) {
				h_vx[sector - 1]->Fill(vx);
			        h_vy[sector - 1]->Fill(vy);
			        h_vz[sector - 1]->Fill(vz);	
	                }
		}
            }
        }

    }

	// Create and save canvas for vx
TCanvas* canvas_vx = new TCanvas("canvas_vx", "Vertex vx", 1200, 800);
canvas_vx->Divide(3, 2); // Divide the canvas into a 3x2 grid
for (int i = 0; i < 6; ++i) {
    canvas_vx->cd(i + 1); // Move to the i-th pad
    h_vx[i]->Draw();
}
canvas_vx->Update();
canvas_vx->SaveAs("Vertex_vx.png");

// Create and save canvas for vy
TCanvas* canvas_vy = new TCanvas("canvas_vy", "Vertex vy", 1200, 800);
canvas_vy->Divide(3, 2); // Divide the canvas into a 3x2 grid
for (int i = 0; i < 6; ++i) {
    canvas_vy->cd(i + 1); // Move to the i-th pad
    h_vy[i]->Draw();
}
canvas_vy->Update();
canvas_vy->SaveAs("Vertex_vy.png");

// Create and save canvas for vz
TCanvas* canvas_vz = new TCanvas("canvas_vz", "Vertex vz", 1200, 800);
canvas_vz->Divide(3, 2); // Divide the canvas into a 3x2 grid
for (int i = 0; i < 6; ++i) {
    canvas_vz->cd(i + 1); // Move to the i-th pad
    h_vz[i]->Draw();
}
canvas_vz->Update();
canvas_vz->SaveAs("Vertex_vz.png");
}
