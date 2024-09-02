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
void SetLorentzVector(TLorentzVector &p4,clas12::region_part_ptr rp)
{
p4.SetXYZM(rp->par()->getPx(), rp->par()->getPy(), rp->par()->getPz(), p4.M());} 

void Q2Calculation (){        

	auto db=TDatabasePDG::Instance();
	TLorentzVector beam(0,0,10.6,10.6);
  	TLorentzVector el(0,0,0,0);
  	TLorentzVector target(0,0,0,db->GetParticle(2212)->Mass());

  	double Q2; 
  //TLorentzVector Q2; 
  TLorentzVector q; 
  const double m_e = db->GetParticle(11)->Mass();
  
  auto* hQ2=new TH1F("Q^2","Q^2",200,0,15);
  int counter=0;
  
  clas12root::HipoChain chain;
  chain.Add("/lustre19/expphy/cache/hallb/scratch/rg-d/production/prod/v4ob_aideLD2/dst/recon/018528/rec_clas_018528.evio.00310-00314.hipo");
  chain.db()->turnOffQADB();
  auto& c12=chain.C12ref();
  while (chain.Next()) {   

  	auto electrons = c12->getByID(11);
	double maxEnergy = -1.0;
	clas12::region_part_ptr maxEnergyElectron = nullptr;

	for (auto& electron : electrons) {
                auto  px = electron->par()->getPx();
                auto py = electron->par()->getPy();
                auto pz = electron->par()->getPz();
		
		double energy = sqrt(px * px + py * py + pz * pz + m_e * m_e);

		if (energy > maxEnergy) {
                    maxEnergy = energy;
	  	    maxEnergyElectron = electron;
		}
	}

	if (maxEnergyElectron) {
                SetLorentzVector(el, maxEnergyElectron); // Use the max energy electron
                q = beam - el;
                Q2 = -q.Mag2();
                hQ2->Fill(Q2);
         }	  
       

   counter++;
 }
 TCanvas* can=new TCanvas();
 hQ2->Draw(); 
}
