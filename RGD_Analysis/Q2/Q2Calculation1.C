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
#include <TChain.h>
#include <TCanvas.h>
#include <TBenchmark.h>
#include "clas12reader.h"
#include "HipoChain.h"

using namespace clas12;

void Q2Calculation1() {        

    auto db = TDatabasePDG::Instance();
    TLorentzVector target(0, 0, 0, db->GetParticle(2212)->Mass());
    TLorentzVector beam(0, 0, 10.6, 10.6); // Beam four-momentum
    TLorentzVector pip;
    const double m_pip = db->GetParticle(211)->Mass(); 
    const double m_e = db->GetParticle(11)->Mass(); // Mass of electron
    
    auto* hQ2 = new TH1F("Q^2", "Q^2", 200, 0, 6); // Histogram for Q² values
    auto* hPionPt2 = new TH1F("Pt^2", "Pt^2", 200,0,3);  
    auto* hxB = new TH1F("xB", "xB", 200, 0, 1);
    auto* hz = new TH1F("z", "z", 200, 0, 2); 
    auto* hW = new TH1F("W", "W", 200, 0, 5);
     auto* hnu = new TH1F("nu", "nu", 200, 0, 10);
   auto* hy = new TH1F("y", "y", 200, 0, 1);

    clas12root::HipoChain chain;
    chain.Add("/lustre19/expphy/cache/hallb/scratch/rg-d/production/prod/v4ob_aideLD2/dst/recon/018528/rec_clas_018528.evio.00310-00314.hipo");
    chain.db()->turnOffQADB();
    auto& c12 = chain.C12ref();
    
    while (chain.Next()) {   // Loop over events
        TLorentzVector el; // Reset for the new event 
	double maxEnergy = -1.0;
	double theta; 
	double Q2_value; 
	double nu;
	double pT2_pion; 

        
        auto electrons = c12->getByID(11); // Get electrons (ID 11)
        for (auto& e : electrons) {      // Loop through each electron
            auto px = e->par()->getPx(); 
            auto py = e->par()->getPy(); 
            auto pz = e->par()->getPz();
	     theta = atan2(sqrt(px*px+py*py),pz);
            auto E_s = sqrt(px * px + py * py + pz * pz + m_e * m_e); // Calculate energy
	   	    if (E_s>maxEnergy){             
                    maxEnergy = E_s; 
	    	    el.SetPxPyPzE(px,py,pz,E_s);
		    }
	}
	    if (maxEnergy > -1.0) {
	      nu = beam.E() - el.E();
	     double y = nu / beam.E();    
     	     Q2_value = 4 * beam.E() * el.E() * sin(theta / 2) * sin(theta / 2);
	     double xB = Q2_value / (2 * target.M() * nu);
	     double W = sqrt(target.M()*target.M() + 2*target.M()*nu - Q2_value); // Calculate W^2
	     hnu->Fill(nu); 
	     hW->Fill(W); 
	      hy->Fill(y); 
	     hQ2->Fill(Q2_value);
	     hxB->Fill(xB); 
	   }
	auto ppions = c12->getByID(211); // Get pions (ID 211)
	for (auto& ppion : ppions) { // Loop through each pion
    	auto px_pion = ppion->par()->getPx(); 
    	auto py_pion = ppion->par()->getPy(); 
	auto pz_pion = ppion->par()->getPz(); 
    	pT2_pion = px_pion * px_pion + py_pion * py_pion;
	double energy_pion = sqrt(px_pion * px_pion + py_pion * py_pion + pz_pion * pz_pion + m_pip * m_pip);
	double z = energy_pion / nu;
    	hPionPt2->Fill(pT2_pion);
	hz->Fill(z); 
	
	}



    }
    TCanvas* can = new TCanvas();
    hQ2->Draw(); // Draw the histogram

    TCanvas* can1 = new TCanvas(); 
    hPionPt2->Draw(); 
    
   TCanvas* can2 = new TCanvas(); 
   hxB->Draw(); 

  TCanvas* can3 = new TCanvas(); 
  hz->Draw(); 

  TCanvas* can4 = new TCanvas(); 
  hW->Draw(); 
  
   TCanvas* can5 = new TCanvas(); 
  hnu->Draw(); 
 

  TCanvas* can6 = new TCanvas(); 
  hy->Draw(); 


}

