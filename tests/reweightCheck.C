#include <iostream>
#include <TFile.h>
#include <TTree.h>
#include <TMath.h>
#include <TVector3.h>
#include <vector>
#include "../common.h"
#include <TSystem.h>
#include <TLorentzVector.h>
#include <TH1.h>
#include "../higgs_drv.cc"
//#include "../mylibrary.h

using namespace std;

void reweightCheck(){
  Bool_t debugInfo = false;
  Bool_t debugEventData = false;


  int libpdf_check =  gSystem->Load("../libpdf.so");
  int compiledrv_check = gSystem->CompileMacro("../higgs_drv.cc","0", "higgs_drv");
  
  TFile *fin = TFile::Open("test_inputs/MG5_h0p_converted.root");
  TTree *tin = (TTree*) fin->Get("Delphes");
 
   Float_t         elec_mass;
   Float_t         elec_pt;
   Float_t         elec_px;
   Float_t         elec_py;
   Float_t         elec_pz;
   Float_t         elec_eta;
   Float_t         elec_phi;
   Float_t         elec_e;
   Int_t           elec_pid;
   Float_t         ve_mass;
   Float_t         ve_pt;
   Float_t         ve_px;
   Float_t         ve_py;
   Float_t         ve_pz;
   Float_t         ve_eta;
   Float_t         ve_phi;
   Float_t         ve_e;
   Int_t           ve_pid;
   Float_t         mu_mass;
   Float_t         mu_pt;
   Float_t         mu_px;
   Float_t         mu_py;
   Float_t         mu_pz;
   Float_t         mu_eta;
   Float_t         mu_phi;
   Float_t         mu_e;
   Int_t           mu_pid;
   Float_t         vm_mass;
   Float_t         vm_pt;
   Float_t         vm_px;
   Float_t         vm_py;
   Float_t         vm_pz;
   Float_t         vm_eta;
   Float_t         vm_phi;
   Float_t         vm_e;
   Int_t           vm_pid;
   Float_t         g1_mass;
   Float_t         g1_pt;
   Float_t         g1_px;
   Float_t         g1_py;
   Float_t         g1_pz;
   Float_t         g1_eta;
   Float_t         g1_phi;
   Float_t         g1_e;
   Int_t           g1_pid;
   Float_t         g2_mass;
   Float_t         g2_pt;
   Float_t         g2_px;
   Float_t         g2_py;
   Float_t         g2_pz;
   Float_t         g2_eta;
   Float_t         g2_phi;
   Float_t         g2_e;
   Int_t           g2_pid;
   Float_t         h_mass;
   Float_t         h_pt;
   Float_t         h_px;
   Float_t         h_py;
   Float_t         h_pz;
   Float_t         h_eta;
   Float_t         h_phi;
   Float_t         h_e;
   Int_t           h_pid;

   // List of branches
   TBranch        *b_elec_mass;   //!
   TBranch        *b_elec_pt;   //!
   TBranch        *b_elec_px;   //!
   TBranch        *b_elec_py;   //!
   TBranch        *b_elec_pz;   //!
   TBranch        *b_elec_eta;   //!
   TBranch        *b_elec_phi;   //!
   TBranch        *b_elec_e;   //!
   TBranch        *b_elec_pid;   //!
   TBranch        *b_ve_mass;   //!
   TBranch        *b_ve_pt;   //!
   TBranch        *b_ve_px;   //!
   TBranch        *b_ve_py;   //!
   TBranch        *b_ve_pz;   //!
   TBranch        *b_ve_eta;   //!
   TBranch        *b_ve_phi;   //!
   TBranch        *b_ve_e;   //!
   TBranch        *b_ve_pid;   //!
   TBranch        *b_mu_mass;   //!
   TBranch        *b_mu_pt;   //!
   TBranch        *b_mu_px;   //!
   TBranch        *b_mu_py;   //!
   TBranch        *b_mu_pz;   //!
   TBranch        *b_mu_eta;   //!
   TBranch        *b_mu_phi;   //!
   TBranch        *b_mu_e;   //!
   TBranch        *b_mu_pid;   //!
   TBranch        *b_vm_mass;   //!
   TBranch        *b_vm_pt;   //!
   TBranch        *b_vm_px;   //!
   TBranch        *b_vm_py;   //!
   TBranch        *b_vm_pz;   //!
   TBranch        *b_vm_eta;   //!
   TBranch        *b_vm_phi;   //!
   TBranch        *b_vm_e;   //!
   TBranch        *b_vm_pid;   //!
   TBranch        *b_g1_mass;   //!
   TBranch        *b_g1_pt;   //!
   TBranch        *b_g1_px;   //!
   TBranch        *b_g1_py;   //!
   TBranch        *b_g1_pz;   //!
   TBranch        *b_g1_eta;   //!
   TBranch        *b_g1_phi;   //!
   TBranch        *b_g1_e;   //!
   TBranch        *b_g1_pid;   //!
   TBranch        *b_g2_mass;   //!
   TBranch        *b_g2_pt;   //!
   TBranch        *b_g2_px;   //!
   TBranch        *b_g2_py;   //!
   TBranch        *b_g2_pz;   //!
   TBranch        *b_g2_eta;   //!
   TBranch        *b_g2_phi;   //!
   TBranch        *b_g2_e;   //!
   TBranch        *b_g2_pid;   //!
   TBranch        *b_h_mass;   //!
   TBranch        *b_h_pt;   //!
   TBranch        *b_h_px;   //!
   TBranch        *b_h_py;   //!
   TBranch        *b_h_pz;   //!
   TBranch        *b_h_eta;   //!
   TBranch        *b_h_phi;   //!
   TBranch        *b_h_e;   //!
   TBranch        *b_h_pid;   //!

  tin->SetBranchAddress("elec_mass", &elec_mass, &b_elec_mass);
   tin->SetBranchAddress("elec_pt", &elec_pt, &b_elec_pt);
   tin->SetBranchAddress("elec_px", &elec_px, &b_elec_px);
   tin->SetBranchAddress("elec_py", &elec_py, &b_elec_py);
   tin->SetBranchAddress("elec_pz", &elec_pz, &b_elec_pz);
   tin->SetBranchAddress("elec_eta", &elec_eta, &b_elec_eta);
   tin->SetBranchAddress("elec_phi", &elec_phi, &b_elec_phi);
   tin->SetBranchAddress("elec_e", &elec_e, &b_elec_e);
   tin->SetBranchAddress("elec_pid", &elec_pid, &b_elec_pid);
   tin->SetBranchAddress("ve_mass", &ve_mass, &b_ve_mass);
   tin->SetBranchAddress("ve_pt", &ve_pt, &b_ve_pt);
   tin->SetBranchAddress("ve_px", &ve_px, &b_ve_px);
   tin->SetBranchAddress("ve_py", &ve_py, &b_ve_py);
   tin->SetBranchAddress("ve_pz", &ve_pz, &b_ve_pz);
   tin->SetBranchAddress("ve_eta", &ve_eta, &b_ve_eta);
   tin->SetBranchAddress("ve_phi", &ve_phi, &b_ve_phi);
   tin->SetBranchAddress("ve_e", &ve_e, &b_ve_e);
   tin->SetBranchAddress("ve_pid", &ve_pid, &b_ve_pid);
   tin->SetBranchAddress("mu_mass", &mu_mass, &b_mu_mass);
   tin->SetBranchAddress("mu_pt", &mu_pt, &b_mu_pt);
   tin->SetBranchAddress("mu_px", &mu_px, &b_mu_px);
   tin->SetBranchAddress("mu_py", &mu_py, &b_mu_py);
   tin->SetBranchAddress("mu_pz", &mu_pz, &b_mu_pz);
   tin->SetBranchAddress("mu_eta", &mu_eta, &b_mu_eta);
   tin->SetBranchAddress("mu_phi", &mu_phi, &b_mu_phi);
   tin->SetBranchAddress("mu_e", &mu_e, &b_mu_e);
   tin->SetBranchAddress("mu_pid", &mu_pid, &b_mu_pid);
   tin->SetBranchAddress("vm_mass", &vm_mass, &b_vm_mass);
   tin->SetBranchAddress("vm_pt", &vm_pt, &b_vm_pt);
   tin->SetBranchAddress("vm_px", &vm_px, &b_vm_px);
   tin->SetBranchAddress("vm_py", &vm_py, &b_vm_py);
   tin->SetBranchAddress("vm_pz", &vm_pz, &b_vm_pz);
   tin->SetBranchAddress("vm_eta", &vm_eta, &b_vm_eta);
   tin->SetBranchAddress("vm_phi", &vm_phi, &b_vm_phi);
   tin->SetBranchAddress("vm_e", &vm_e, &b_vm_e);
   tin->SetBranchAddress("vm_pid", &vm_pid, &b_vm_pid);
   tin->SetBranchAddress("g1_mass", &g1_mass, &b_g1_mass);
   tin->SetBranchAddress("g1_pt", &g1_pt, &b_g1_pt);
   tin->SetBranchAddress("g1_px", &g1_px, &b_g1_px);
   tin->SetBranchAddress("g1_py", &g1_py, &b_g1_py);
   tin->SetBranchAddress("g1_pz", &g1_pz, &b_g1_pz);
   tin->SetBranchAddress("g1_eta", &g1_eta, &b_g1_eta);
   tin->SetBranchAddress("g1_phi", &g1_phi, &b_g1_phi);
   tin->SetBranchAddress("g1_e", &g1_e, &b_g1_e);
   tin->SetBranchAddress("g1_pid", &g1_pid, &b_g1_pid);
   tin->SetBranchAddress("g2_mass", &g2_mass, &b_g2_mass);
   tin->SetBranchAddress("g2_pt", &g2_pt, &b_g2_pt);
   tin->SetBranchAddress("g2_px", &g2_px, &b_g2_px);
   tin->SetBranchAddress("g2_py", &g2_py, &b_g2_py);
   tin->SetBranchAddress("g2_pz", &g2_pz, &b_g2_pz);
   tin->SetBranchAddress("g2_eta", &g2_eta, &b_g2_eta);
   tin->SetBranchAddress("g2_phi", &g2_phi, &b_g2_phi);
   tin->SetBranchAddress("g2_e", &g2_e, &b_g2_e);
   tin->SetBranchAddress("g2_pid", &g2_pid, &b_g2_pid);
   tin->SetBranchAddress("h_mass", &h_mass, &b_h_mass);
   tin->SetBranchAddress("h_pt", &h_pt, &b_h_pt);
   tin->SetBranchAddress("h_px", &h_px, &b_h_px);
   tin->SetBranchAddress("h_py", &h_py, &b_h_py);
   tin->SetBranchAddress("h_pz", &h_pz, &b_h_pz);
   tin->SetBranchAddress("h_eta", &h_eta, &b_h_eta);
   tin->SetBranchAddress("h_phi", &h_phi, &b_h_phi);
   tin->SetBranchAddress("h_e", &h_e, &b_h_e);
   tin->SetBranchAddress("h_pid", &h_pid, &b_h_pid);

   afloat_t event_data[24];
   TFile *fout = new TFile("MG5_h0p_NLO_to_h0m_histos.root", "RECREATE");
 
   TH1F *h_mll = new TH1F("mll","mll", 20,0,100);
   TH1F *h_mll_rew03 = new TH1F("mllrew03","mllrew03", 20,0,100);
   TH1F *h_mll_rew0 = new TH1F("mllrew0","mllrew0", 20,0,100);

   TH1F *h_ptll = new TH1F("ptll","ptll", 20,0,120);
   TH1F *h_ptll_rew03 = new TH1F("ptllrew03","ptllrew03", 20,0,120);
   TH1F *h_ptll_rew0 = new TH1F("ptllrew0","ptllrew0", 20,0,120);

   TH1F *h_dphill = new TH1F("dphill","dphill", 16,0,3.2);
   TH1F *h_dphill_rew03 = new TH1F("dphillrew03","dphillrew03", 16,0,3.2);
   TH1F *h_dphill_rew0 = new TH1F("dphillrew0","dphillrew0", 16,0,3.2);

   TH1F *h_dpt = new TH1F("dpt","dpt",          20,0,60);
   TH1F *h_dpt_rew03 = new TH1F("dptrew03","dptrew03",          20,0,60);
   TH1F *h_dpt_rew0 = new TH1F("dptrew0","dptrew0",          20,0,60);


   TH1F *h_ratio = new TH1F("ratio","ratio", 100,0,10);
   TH1F *h_ME_ca1p = new TH1F("ME_ca1p", "ME_ca1p", 15,0, 0.0000015);
   TH1F *h_ME_ca1p_small = new TH1F("ME_ca1p_small", "ME_ca1p_small", 10,0, 0.000001);
   TH1F *h_ME_ca00 = new TH1F("ME_ca00", "ME_ca00", 15,0, 0.0000000015);
   TH1F *h_ME_ca00_small = new TH1F("ME_ca00_small", "ME_ca00_small", 10,0, 0.000000001);
   TH1F *h_ME_ca03p = new TH1F("ME_ca03p", "ME_ca03p", 15,0, 0.0000015);


   TH1F *h_logME_ca1p = new TH1F ("logME_ca1p", "logME_ca1p", 15, -30, 0);
   TH1F *h_logME_ca00 = new TH1F ("logME_ca00", "logME_ca00", 15, -30, 0);

   TH1F *h_ME_ca1p_weig = new TH1F("ME_ca1p_weig", "ME_ca1p_weig", 15,0, 0.0000015);
   TH1F *h_ME_ca00_weig = new TH1F("ME_ca00_weig", "ME_ca00_weig", 15,0, 0.0000000015);

   TH1F *h_ME_ca1p_weig_small = new TH1F("ME_ca1p_weig_small", "ME_ca1p_weig_small", 10,0, 0.000001);
   TH1F *h_ME_ca00_weig_small = new TH1F("ME_ca00_weig_small", "ME_ca00_weig_small", 10,0, 0.000000001);

   TH1F *h_logME_ca1p_weig = new TH1F ("logME_ca1p_weig", "logME_ca1p_weig", 15, -30, 0);
   TH1F *h_logME_ca00_weig = new TH1F ("logME_ca00_weig", "logME_ca00_weig", 15, -30, 0);
   
   h_mll->Sumw2();
   h_mll_rew03->Sumw2();
   h_mll_rew0->Sumw2();

   h_ptll->Sumw2();
   h_ptll_rew03->Sumw2();
   h_ptll_rew0->Sumw2();

   h_dphill->Sumw2();
   h_dphill_rew03->Sumw2();
   h_dphill_rew0->Sumw2();

   h_dpt->Sumw2();
   h_dpt_rew03->Sumw2();
   h_dpt_rew0->Sumw2();

   h_ME_ca1p->Sumw2();
   h_ME_ca00->Sumw2();

   h_ME_ca1p_small->Sumw2();
   h_ME_ca00_small->Sumw2();

   h_logME_ca1p->Sumw2();
   h_logME_ca00->Sumw2();

   h_ME_ca1p_weig->Sumw2();
   h_ME_ca00_weig->Sumw2();

   h_ME_ca1p_weig_small->Sumw2();
   h_ME_ca00_weig_small->Sumw2();

   h_logME_ca1p_weig->Sumw2();
   h_logME_ca00_weig->Sumw2();

   Float_t cthstar, cth1, cth2, phi;
   
   Long64_t nentries = tin->GetEntries();
   
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
     tin->GetEntry(jentry);
     
     if(jentry%10000==0) cout << jentry << endl;
     TLorentzVector lm, lp, vbar, v;
     TLorentzVector ga, gb;
     ga.SetPxPyPzE(g1_px, g1_py, g1_pz, g1_e);
     gb.SetPxPyPzE(g2_px, g2_py, g2_pz, g2_e);
     
     // in epmum electron is always -11 (positive charge) and muon 13 (negative charge)
     if(elec_pid>0){
       lm.SetPxPyPzE(elec_px, elec_py, elec_pz, elec_e);
       vbar.SetPxPyPzE(ve_px, ve_py, ve_pz, ve_e);
       lp.SetPxPyPzE(mu_px, mu_py, mu_pz, mu_e);
       v.SetPxPyPzE(vm_px, vm_py, vm_pz, vm_e);
     }else{
       lm.SetPxPyPzE(mu_px, mu_py, mu_pz, mu_e);
       vbar.SetPxPyPzE(vm_px, vm_py, vm_pz, vm_e);
       lp.SetPxPyPzE(elec_px, elec_py, elec_pz, elec_e);
       v.SetPxPyPzE(ve_px, ve_py, ve_pz, ve_e);
     }
       
       
     
     TLorentzVector particles[6]= {ga, gb, lp, v, lm, vbar};
     // fill the event_data array
     for(int ip=0; ip<6; ++ip){
       event_data[ip*4+0] = particles[ip].E();
       event_data[ip*4+1] = particles[ip].Px();
       event_data[ip*4+2] = particles[ip].Py();
       event_data[ip*4+3] = particles[ip].Pz();
     }//for ip    
     
     
     // some selection cuts
     bool good_event = true;
     if(TMath::Abs(elec_eta)>2.47) good_event=false;
     if(TMath::Abs(elec_eta)<1.52 && TMath::Abs(elec_eta)>1.37) good_event=false;
     if(TMath::Abs(mu_eta)>2.5) good_event=false;
     
     if(elec_pt > mu_pt ){
       if(elec_pt<22.) good_event=false;
       if(mu_pt<15.) good_event=false;
     }
     else{
       if(elec_pt<15.) good_event=false;
       if(mu_pt<22.) good_event=false;
     }

     //if(!good_event){
       //  cout << elec_pt <<"\t" << mu_pt << "\t" << elec_eta <<"\t" << mu_eta << endl;
     //}
     if(!good_event) continue;



     TLorentzVector dilep;
     dilep = lp+lm;
     Float_t mll, ptll, dpt, dphill;
     mll=dilep.M();
     ptll=dilep.Pt();

     h_mll->Fill(mll);
     h_ptll->Fill(ptll);

     dpt=TMath::Abs(elec_pt-mu_pt);
     dphill=lp.DeltaPhi(lm);

     h_dpt->Fill(dpt);
     h_dphill->Fill(dphill);

     afloat_t result[1];
     afloat_t result03[1];
     afloat_t result00[1];
     afloat_t ca   = 1.0;
     afloat_t kSM  = 1.0;
     afloat_t kHWW = 0.0;
     afloat_t kAWW = 0.0;
     afloat_t kHdw = 0.0;
     afloat_t Lambda = 1000.;
     quadme(result, event_data, 1,0,0,0,0,0,0,0,0,0,0, ca, kSM, kHWW, kAWW, kHdw, Lambda, -1);

     ca = 0.3;
     kSM = 1.0;
     kHWW = 2.0;
     kAWW = 2.0;
     kHdw = 1.0;
     Lambda = 125.5;
     quadme(result03, event_data, 1,0,0,0,0,0,0,0,0,0,0, ca, kSM, kHWW, kAWW, kHdw, Lambda, -1);

     ca = 0.0;
     kSM = 1.0;
     kHWW = 0.0;
     kAWW = 1.0;
     kHdw = 0.0;
     Lambda = 1000.;
     quadme(result00, event_data, 1,0,0,0,0,0,0,0,0,0,0, ca, kSM, kHWW, kAWW, kHdw, Lambda, -1);
    

     afloat_t norm = 1.9336535e-09;
     afloat_t reweight1p = (result03[0]/result[0]);
     afloat_t reweight0 = (result00[0]/result[0]);

     //   cout << result00[0] << endl;

     //cout << result[0] << endl;
     h_ME_ca1p->Fill(result[0]);
     h_ME_ca00->Fill(result00[0]);

     h_ME_ca1p_small->Fill(result[0]);
     h_ME_ca00_small->Fill(result00[0]);

     h_mll_rew03->Fill(mll, reweight1p);
     h_mll_rew0->Fill(mll, reweight0);

     h_dphill_rew03->Fill(dphill, reweight1p);
     h_dphill_rew0->Fill(dphill, reweight0);

     h_ptll_rew03->Fill(ptll, reweight1p);
     h_ptll_rew0->Fill(ptll, reweight0);

     h_dpt_rew03->Fill(dpt, reweight1p);
     h_dpt_rew0->Fill(dpt, reweight0);

     h_ME_ca03p->Fill(result03[0]);

     h_logME_ca1p->Fill(TMath::Log(result[0]));
     h_logME_ca00->Fill(TMath::Log(result00[0]));

     h_ME_ca1p_weig->Fill(result[0], reweight1p);
     h_ME_ca00_weig->Fill(result00[0], reweight0);

     h_ME_ca1p_weig_small->Fill(result[0], reweight1p);
     h_ME_ca00_weig_small->Fill(result00[0], reweight0);

     h_logME_ca1p_weig->Fill(TMath::Log(result[0]), reweight1p);
     h_logME_ca00_weig->Fill(TMath::Log(result00[0]), reweight0);



     //h_ratio->Fill(reweight);
     /*
     // Now calculate some ZRF angles:
     TLorentzVector W1 = lp+v;
     TLorentzVector W2 = lm+vbar;
     const TLorentzVector H = W1+W2;
     
     TLorentzVector q;
     q.SetPxPyPzE(0,0, (H.E()+H.Pz())/2. , (H.E()+H.Pz())/2.);
     q.Boost(-(H.BoostVector()));

     const TVector3 parton = q.Vect().Unit();

     W1.Boost(-(H.BoostVector()));
     W2.Boost(-(H.BoostVector()));


     const TVector3 w1 = W1.Vect().Unit();
     const TVector3 w2 = W2.Vect().Unit();

     //costh*
     cthstr = TMath::Cos(parton.Angle(w1));
     
     TLorentzVector nlp(lp), nlm(lm), nv(v), nvbar(vbar);
     nlp.Boost(-(H.BoostVector()));
     nlm.Boost(-(H.BoostVector()));
     nv.Boost(-(H.BoostVector()));
     nvbar.Boost(-(H.BoostVector()));
     
     const TVector3 veclp = nlp.Vect();
     const TVector3 veclm = nlm.Vect();
     const TVector3 vecv = nv.Vect();
     const TVector3 vecvbar = nvbar.Vect();
     const TVector3 vecz(0.,0.,1.);

     //phi phi1
     const TVector3
       */

   }//for nentries


      fout->cd();

   h_ME_ca1p->Write();
   h_ME_ca00->Write();
   h_ME_ca1p_small->Write();
   h_ME_ca00_small->Write();

   h_mll->Write();
   h_ptll->Write();
   h_dphill->Write();
   h_dpt->Write();

   h_mll_rew03->Write();
   h_mll_rew0->Write();
   h_dphill_rew03->Write();
   h_dphill_rew0->Write();
   h_ptll_rew03->Write();
   h_ptll_rew0->Write();
   h_dpt_rew03->Write();
   h_dpt_rew0->Write();
   h_ME_ca03p->Write();
   h_logME_ca1p->Write();
   h_ME_ca1p_weig->Write();
   h_ME_ca00_weig->Write();
   h_ME_ca1p_weig_small->Write();
   h_ME_ca00_weig_small->Write();

   h_logME_ca1p_weig->Write();
   h_logME_ca00_weig->Write();
   // h_ratio->Write();
   

    fout->Close();
  
   //   makeRatioPlot(h_ME_ca1p, h_ME_ca1p_weig, "0^{+}", "Same Reweig 0^{+}", "meorig", "merew", "me", "log(ME)", "Norm. Entries", "Orig/Rew","./");



}// end of script
