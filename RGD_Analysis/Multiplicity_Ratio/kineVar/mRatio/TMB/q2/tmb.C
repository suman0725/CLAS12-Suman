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
#include <TGraphErrors.h>
#include <TLatex.h>
#include <vector>


namespace fs = std::filesystem;
using namespace clas12;

void SetLorentzVector(TLorentzVector &p4, clas12::region_part_ptr rp) {
    p4.SetXYZM(rp->par()->getPx(), rp->par()->getPy(), rp->par()->getPz(), p4.M());
}

double CalculatePt2(TLorentzVector P_h, TLorentzVector q) {
    TVector3 q3 = q.Vect(); 
    TVector3 p_h3 = P_h.Vect(); 
    TVector3 cross_qph = q3.Cross(p_h3);
    double magnitudeCrossProductSquared = cross_qph.Mag2(); // |q × p_h|^2
    double magnitudeQSquared = q3.Mag2();                     // |q|^2
    return magnitudeCrossProductSquared / magnitudeQSquared; // Pt^2
}




/* void addHipoFiles(clas12root::HipoChain& chain, const fs::path& baseDir) {
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
} */

void tmb() {
    auto db = TDatabasePDG::Instance();
    double beamEnergy = 10.532; // GeV
    TLorentzVector beam(0, 0, beamEnergy, beamEnergy); // Beam vector
    TLorentzVector target(0, 0, 0, db->GetParticle(2212)->Mass()); // Proton target
    TLorentzVector el, pip;

    // Histograms for Q^2 for pions
    /* auto* h_LD2_Q2_piP = new TH1F("h_LD2_Q2_pion", "LD2 Q^{2} (Pions);Q^{2} (GeV)^{2};Counts", 200, 0, 8.5);
    auto* h_C_Q2_piP = new TH1F("h_C_Q2_pion", "Carbon Q^{2} (Pions);Q^{2} (GeV)^{2};Counts", 200, 0, 8.5);
    */
    // Binning for Q^2
    int nBins = 8;
    double q2Min = 0;
    double q2Max = 8;
    double binWidth = (q2Max - q2Min) / nBins;

    std::vector<double> pT2_D(nBins, 0.0); // Store <Pt^2> for LD2
    std::vector<double> pT2_A(nBins, 0.0); // Store <Pt^2> for Carbon
    std::vector<double> pT4_D(nBins, 0.0); // Store <Pt^4> for LD2
    std::vector<double> pT4_A(nBins, 0.0); // Store <Pt^4> for Carbon
    std::vector<int> counts_D(nBins, 0);   // Count events for LD2
    std::vector<int> counts_A(nBins, 0);   // Count events for Carbon
    int totalEntries = 0;
    int totalEntries_C = 0;
    std::vector<double> errors_C(nBins, 0.0);
    std::vector<double> errors_D(nBins, 0.0);
    std::vector<double> errors(nBins, 0.0);  // Standard error for each bin
    std::vector<double> average_pT2(nBins, 0.0);
    std::vector<double> average_pT2_C(nBins, 0.0);
    std::vector<double> deltapT2(nBins, 0.0);
    std::vector<double> q2Values(nBins);
    std::vector<double> average_pT4(nBins, 0.0);   // For LD2
    std::vector<double> average_pT4_C(nBins, 0.0); // For Carbon
    std::vector<double> variance(nBins, 0.0);      // Variance for LD2
    std::vector<double> variance_C(nBins, 0.0);    // Variance for Carbon



    // Hipo files
    clas12root::HipoChain chainLD2, chainC;
    chainLD2.Add("/lustre24/expphy/cache/hallb/scratch/rg-d/production/Bspot/v5dstLD2/dst/recon/018309/rec_clas_018309.evio.00035-00039.hipo");
    chainC.Add("/lustre24/expphy/cache/hallb/scratch/rg-d/production/Bspot/v5dstCxC/dst/recon/018339/rec_clas_018339.evio.00015-00019.hipo");
    chainLD2.db()->turnOffQADB();
    chainC.db()->turnOffQADB();

    // Process LD2 target
    auto& c12 = chainLD2.C12ref();
    while (chainLD2.Next()) {
        auto particles = c12->getDetParticles();
        TLorentzVector el_temp, pip_temp;
        bool foundElectron = false, foundPion = false;

        for (auto& p : particles) {
            int pid = p->par()->getPid();
            int status = p->par()->getStatus();
            int chi2Pid = p->par()->getChi2Pid();

            // Check for electron
            if (pid == 11 && status < 0 && chi2Pid < 5 && chi2Pid > -5) {
                SetLorentzVector(el_temp, p);
                foundElectron = true;
            }

            // Check for pion
            if (pid == 211 && chi2Pid > -10 && chi2Pid < 10) {
                SetLorentzVector(pip_temp, p);
                foundPion = true;
            }
        }

        // Calculate nu, q, and Q^2 if an electron is found
        if (foundElectron) {
            double nu = beam.Energy() - el_temp.Energy(); 
            TLorentzVector q = beam - el_temp; 
            double Q2 = -q.Mag2();

            // If both electron and pion are found, calculate Pt^2
            if (foundPion) {
                int binIndex = static_cast<int>(Q2 / binWidth);
                if (binIndex >= 0 && binIndex < nBins) {
                    double pT2_value = CalculatePt2(pip_temp, q);
                    pT2_D[binIndex] += pT2_value;
                    pT4_D[binIndex] += pT2_value*pT2_value; 
                    counts_D[binIndex]++;
                    totalEntries ++;
                    //h_LD2_Q2_piP->Fill(Q2);

                
                }
            }  
        }
    }

                for (int i = 0; i < nBins; ++i) {
                    if (counts_D[i] > 0) {

                        average_pT2[i] = pT2_D[i] / totalEntries; 
                        average_pT4[i] = pT4_D[i] / totalEntries;               
                        variance[i] = (average_pT4[i] - average_pT2[i] * average_pT2[i]) / totalEntries ;
                        errors_D[i] = (variance[i] > 0) ? std::sqrt(variance[i] / counts_D[i]) : 0.0;
                    }
                        
                }

    // Process Carbon target
    auto& c12_C = chainC.C12ref();
    while (chainC.Next()) {
        auto particles = c12_C->getDetParticles();
        TLorentzVector el_temp_C, pip_temp_C;
        bool foundElectron_C = false, foundPion_C = false;

        for (auto& p : particles) {
            int pid = p->par()->getPid();
            int status = p->par()->getStatus();
            int chi2Pid = p->par()->getChi2Pid();

            // Check for electron
            if (pid == 11 && status < 0 && chi2Pid < 5 && chi2Pid > -5) {
                SetLorentzVector(el_temp_C, p);
                foundElectron_C = true;
            }

            // Check for pion
            if (pid == 211 && chi2Pid > -10 && chi2Pid < 10) {
                SetLorentzVector(pip_temp_C, p);
                foundPion_C = true;
            }
        }

        // Calculate nu, q, and Q^2 if an electron is found
        if (foundElectron_C) {
            double nu_C = beam.Energy() - el_temp_C.Energy(); 
            TLorentzVector q_C = beam - el_temp_C; 
            double Q2_C = -q_C.Mag2();

            // If both electron and pion are found, calculate Pt^2
            if (foundPion_C) {
                int binIndex_C = static_cast<int>(Q2_C / binWidth);
                if (binIndex_C >= 0 && binIndex_C < nBins) {
                    double pT2_value_C = CalculatePt2(pip_temp_C, q_C);
                    pT2_A[binIndex_C] += pT2_value_C;
                    pT4_A[binIndex_C] += pT2_value_C*pT2_value_C; 
                    counts_A[binIndex_C]++;
                    totalEntries_C++;
                    //h_C_Q2_piP->Fill(Q2_C);
                    
                    
                }
            }
        }
    }


                for (int i = 0; i < nBins; ++i) {
                    if (counts_A[i] > 0) {

                        average_pT2_C[i] = pT2_A[i] / totalEntries_C; 
                        average_pT4_C[i] = pT4_A[i] / totalEntries_C ;  
                        variance_C[i] = (average_pT4_C[i] - average_pT2_C[i] * average_pT2_C[i]) / totalEntries_C ;
                        errors_C[i] = (variance_C[i] > 0) ? std::sqrt(variance_C[i] / counts_A[i]) : 0.0;
                    }
                        
                }
                
    TGraphErrors* graph = new TGraphErrors(nBins);
    for (int i = 0; i < nBins; ++i){
        deltapT2[i] = average_pT2_C[i] - average_pT2[i]; 
        errors[i] = sqrt(errors_C[i]*errors_C[i] + errors_D[i]*errors_D[i]) ;
        q2Values[i] = q2Min + (i + 0.5) * binWidth;
        graph->SetPoint(i, q2Values[i], deltapT2[i]);  // Set Q² value and deltapT2
        graph->SetPointError(i, 0, errors[i]);  // Set the error for each point (0 for x-error)
        std::cout << "Bin " << i << ": deltapT2 = " << deltapT2[i] << ", Error = " << errors[i] << std::endl;
    }

    
    graph->SetTitle("Scatter Plot of #Delta p_T^2 vs Q^{2}");
    graph->GetXaxis()->SetTitle("Q^{2} (GeV)^{2}");
    graph->GetYaxis()->SetTitle("#Delta p_T^2");
    // Customize the marker style and size
    graph->SetMarkerStyle(20);  // Solid circle marker
    graph->SetMarkerSize(1.2);  // Adjust marker size as needed
    graph->SetMarkerColor(kBlue);  // Set marker color (black)

    // Optionally customize the line for error bars
    graph->SetLineWidth(1);  // Set the width of the error bars
    graph->SetLineColor(kBlack);  // Set the color of the error bars (black)

    auto* c1 = new TCanvas("c1", "Scatter Plot", 800, 600);
    graph->Draw("AP");  // Draw with axes and points




    
}
