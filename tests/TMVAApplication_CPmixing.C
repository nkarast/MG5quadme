#include <cstdlib>
#include <iostream>
#include <map>
#include <string>
#include <math.h>
#include <cmath>

#include "TChain.h"
#include "TFile.h"
#include "TTree.h"
#include "TBranch.h"
#include "TString.h"
#include "TObjString.h"
#include "TSystem.h"
#include "TROOT.h"
#include "TH1.h"

#if not defined(__CINT__) || defined (__MAKECINT__)
#include "TMVA/Tools.h"
#include "TMVA/Reader.h"
#include "TMVA/MethodCuts.h"
#endif

#include "TCanvas.h"
#include "THStack.h"
#include "TLatex.h"
#include "TStyle.h"
#include "TMath.h"

// -- nikos CP mixing
#include <vector>
#include  <TLorentzVector.h>
#include "reweighting/common.h"
#include "reweighting/higgs_drv.cc"

using namespace std;

void TMVAApplication_CPmixing (TString HiggsMass = "130", TString analysis = "BDT3", TString ana_name = "BDT_spin0", bool runEvenOdd = false, 
															 TString ana2_name = "BDT_spin2", TString sub_ana_name = "BDT_spin0_wjets", 
															 TString sys = "Nominal", TString weightsdir = "weights_spin0", TString weightsdir2 = "weights_spin2", 
															 TString sub_weightsdir = "weights_spin0_wjets", TString input_dir = "inputs/") {
    
  std::cout << "HiggsMass = " << HiggsMass << ", sys = " << sys << std::endl;
  std::cout << "analysis = " << ana_name << ", analysis2 = " << ana2_name << std::endl;
  std::cout << "sub analysis = " << sub_ana_name << std::endl;
  std::cout << "weightsdir = " << weightsdir << std::endl;
	std::cout << "weightsdir2 = " << weightsdir2 << std::endl;
    
  float invpb=20280.2; //20693.7;//13232.9;//19876.4;//12980.2 (HCP)
  bool doDebug = false;
  bool doNtuplesForLimit = true; // if this is true, you will get ntuples ready for limit setting, organized in the correct structure!
	bool doMVAPropertiesCuts = true;
	bool useXjetsSubtraining = false;
	bool runSpin2training = true;
  bool doTopDataDrivenNtuples = false; 
	bool doPtllReweighting = true;
  TString trainodd = "trainodd_";
  TString traineven = "traineven_";
  bool dokHWWScan = true;
    
	TMVA::Tools::Instance();
	std::cout << std::endl << "Start!" << std::endl;
    
	TString output_dir = "./";
    
	Int_t sampN = 29;
	TString fnames[sampN];
	fnames[0] = "ggH"+HiggsMass+"_EF";
	fnames[1] = "ttbar_dilep_powheg_AFII";
	fnames[2] = "singletop_powheg_AFII";
	fnames[3] = "WW_gg";
	fnames[4] = "WW_qq_powheg_py6perugia2011c";
	fnames[5] = "WW2jet_sherpa";
	fnames[6] = "Wgamma_alpgen";
	fnames[7] = "Wgammastar_sherpa";
	fnames[8] = "ZZ_powheg";
	fnames[9] = "ZZ_gg2ZZ";
	fnames[10] = "ZZ_sherpa2j";
	fnames[11] = "WZ_powheg";
	fnames[12] = "WZ_sherpa2j";
	fnames[13] = "Zjets";
	fnames[14] = "Zjets_tautau";
	fnames[15] = "DrellYan";
	fnames[16] = "DrellYan_tautau";
	fnames[17] = "Zgamma_sherpa";
	fnames[18] = "Zgamma_tautau_sherpa";
	fnames[19] = "ZjetsEW_sherpa_jj";
	fnames[20] = "ZjetsEW_tautau_sherpa_jj";
	fnames[21] = "ggH125_CPMixing_mg5";
	fnames[22] = "ggH125_CPMixing_mg5";
	fnames[23] = "ggH125_CPMixing_mg5";
	fnames[24] = "ggH125_CPMixing_mg5";
	fnames[25] = "ggH125_CPMixing_mg5";
	fnames[26] = "ggH125_CPMixing_mg5";
	fnames[27] = "Wjets_datadriven";
	fnames[28] = "data_Run1";

  //fnames[20] = "data_Moriond_Blind";
  // fnames[1] = "VBF"+HiggsMass+"_EF";
  // fnames[2] = "WH"+HiggsMass;
  // fnames[3] = "ZH"+HiggsMass;
  // fnames[16] = "ggH125_spin0p";

    
	TString hnames[sampN];
	hnames[0] = "ggf"+HiggsMass;
	hnames[1] = "ttbar";
	hnames[2] = "st";
	hnames[3] = "wwgg";
	hnames[4] = "wwqq";
	hnames[5] = "ww2j";
	hnames[6] = "wg";
	hnames[7] = "wgs";
	hnames[8] = "zz";
	hnames[9] = "zzgg";
	hnames[10] = "zz2j";
	hnames[11] = "wz";
	hnames[12] = "wz2j";
	hnames[13] = "zhighjets";
	hnames[14] = "zhighjetstautau";
	hnames[15] = "dy";
	hnames[16] = "dytautau";
	hnames[17] = "zg";
	hnames[18] = "zgtautau";
	hnames[19] = "zew";
	hnames[20] = "zewtautau";
	if(dokHWWScan) {
		hnames[21] = "ggH125kHWW0dot1p";
		hnames[22] = "ggH125kHWW0dot1m";
		hnames[23] = "ggH125kHWW0dot15p";
		hnames[24] = "ggH125kHWW0dot15m";
		hnames[25] = "ggH125kHWW0dot2p";
		hnames[26] = "ggH125kHWW0dot2m";
	} else {
		hnames[21] = "ggH125mixca08m";
		hnames[22] = "ggH125mixca085m";
		hnames[23] = "ggH125mixca09m";
		hnames[24] = "ggH125mixca095";
		hnames[25] = "ggH125mixca098";
		hnames[26] = "ggH125mixca099";
	}
	hnames[27] = "wjets";
	hnames[28] = "data";
	//hnames[19] = "dataBlind";
	// hnames[1] = "vbf"+HiggsMass;
	// hnames[2] = "wh"+HiggsMass;
	// hnames[3] = "zh"+HiggsMass;
	// hnames[16] = "ggH125spin0p";

    
	if(sys!="Nominal") sampN = 27;
	if(sys=="EL_SysFakeUp_MVA" || sys=="EL_SysFakeDown_MVA" || sys=="MU_SysFakeUp_MVA" || sys=="MU_SysFakeDown_MVA") sampN = 28;
	output_dir = output_dir+"H"+HiggsMass;
	if(doNtuplesForLimit){
	  if(sys=="Nominal") output_dir += "/Nominal/Normal/";
	  else if(sys=="ElecResolutionUp") output_dir += "/ATLAS_EL_RES/ElecResolutionUp/";
	  else if(sys=="ElecResolutionDown") output_dir += "/ATLAS_EL_RES/ElecResolutionDown/";
	  else if(sys=="ElecScaleUp") output_dir += "/ATLAS_EL_ESCALE/ElecScaleUp/";
	  else if(sys=="ElecScaleDown") output_dir += "/ATLAS_EL_ESCALE/ElecScaleDown/";
	  else if(sys=="IDLOW") output_dir += "/ATLAS_MU_ID_RES/IDLOW/";
	  else if(sys=="IDUP") output_dir += "/ATLAS_MU_ID_RES/IDUP/";
	  else if(sys=="JERUp") output_dir += "/ATLAS_JER/JERUp/";
	  else if(sys=="JESDown") output_dir += "/ATLAS_JES/JESDown/";
	  else if(sys=="JESUp") output_dir += "/ATLAS_JES/JESUp/";
	  else if(sys=="MSLOW") output_dir += "/ATLAS_MU_MS_RES/MSLOW/";
	  else if(sys=="MSUP") output_dir += "/ATLAS_MU_MS_RES/MSUP/";
	  else if(sys=="MUONSCALELOW") output_dir += "/ATLAS_MU_ESCALE/MUONSCALELOW/";
	  else if(sys=="MUONSCALEUP") output_dir += "/ATLAS_MU_ESCALE/MUONSCALEUP/";
	  else if(sys=="ResoSoftTermsDown_ptHard") output_dir += "/ATLAS_MET_RESOSOFT/ResoSoftTermsDown_ptHard/";
	  else if(sys=="ResoSoftTermsUp_ptHard") output_dir += "/ATLAS_MET_RESOSOFT/ResoSoftTermsUp_ptHard/";
	  else if(sys=="ScaleSoftTermsDown_ptHard") output_dir += "/ATLAS_MET_SCALESOFT/ScaleSoftTermsDown_ptHard/";
	  else if(sys=="ScaleSoftTermsUp_ptHard") output_dir += "/ATLAS_MET_SCALESOFT/ScaleSoftTermsUp_ptHard/";
		else if(sys=="EL_SysFakeUp_MVA") output_dir += "/FakeRate_EL_HWW/SysFakeUp/";
		else if(sys=="EL_SysFakeDown_MVA") output_dir += "/FakeRate_EL_HWW/SysFakeDown/";
		else if(sys=="MU_SysFakeUp_MVA") output_dir += "/FakeRate_MU_HWW/SysFakeUp/";
		else if(sys=="MU_SysFakeDown_MVA") output_dir += "/FakeRate_MU_HWW/SysFakeDown/";
	  else if(sys=="ResoSoftTrackMetDown_ptHard") output_dir += "/ATLAS_TRACKMET_RESOSOFT/ResoSoftTrackMetDown_ptHard/";
	  else if(sys=="ResoSoftTrackMetUp_ptHard") output_dir += "/ATLAS_TRACKMET_RESOSOFT/ResoSoftTrackMetUp_ptHard/";
	  else if(sys=="ScaleSoftTrackMetDown_ptHard") output_dir += "/ATLAS_TRACKMET_SCALESOFT/ScaleSoftTrackMetDown_ptHard/";
	  else if(sys=="ScaleSoftTrackMetUp_ptHard") output_dir += "/ATLAS_TRACKMET_SCALESOFT/ScaleSoftTrackMetUp_ptHard/";
	  else if(sys=="FlavJESDown") output_dir += "/ATLAS_JES_FLAV/FlavJESDown/";
	  else if(sys=="FlavJESUp") output_dir += "/ATLAS_JES_FLAV/FlavJESUp/";
	  else if(sys=="bJESDown") output_dir += "/ATLAS_JES_BJET/bJESDown/";
	  else if(sys=="bJESUp") output_dir += "/ATLAS_JES_BJET/bJESUp/";
	  else if(sys=="FWDJESDown") output_dir += "/ATLAS_JES_FWD/FWDJESDown/";
	  else if(sys=="FWDJESUp") output_dir += "/ATLAS_JES_FWD/FWDJESUp/";
	  else if(sys=="NPVJESDown") output_dir += "/ATLAS_JES_NPV/NPVJESDown/";
	  else if(sys=="NPVJESUp") output_dir += "/ATLAS_JES_NPV/NPVJESUp/";
	  else if(sys=="MuJESDown") output_dir += "/ATLAS_JES_MU/MuJESDown/";
	  else if(sys=="MuJESUp") output_dir += "/ATLAS_JES_MU/MuJESUp/";
	  else if(sys=="FlavRespJESDown") output_dir += "/ATLAS_JES_FlavResp/FlavRespJESDown/";
	  else if(sys=="FlavRespJESUp") output_dir += "/ATLAS_JES_FlavResp/FlavRespJESUp/";
	  else if(sys=="PileRhoJESDown") output_dir += "/ATLAS_JES_2012_PileRho_HWW/PileRhoJESDown/";
	  else if(sys=="PileRhoJESUp") output_dir += "/ATLAS_JES_2012_PileRho_HWW/PileRhoJESUp/";
	  else if(sys=="HighPtJESDown") output_dir += "/ATLAS_JES_HighPt/HighPtJESDown/";
	  else if(sys=="HighPtJESUp") output_dir += "/ATLAS_JES_HighPt/HighPtJESUp/";
	  else if(sys=="FlavCompJESDown") output_dir += "/ATLAS_JES_FlavComp_HWW_other/FlavCompJESDown/";
	  else if(sys=="FlavCompJESUp") output_dir += "/ATLAS_JES_FlavComp_HWW_other/FlavCompJESUp/";
	  else if(sys=="PilePtJESDown") output_dir += "/ATLAS_JES_2012_PilePt/PilePtJESDown/";
	  else if(sys=="PilePtJESUp") output_dir += "/ATLAS_JES_2012_PilePt/PilePtJESUp/";
	  else if(sys=="Eta_ModellingJESDown") output_dir += "/ATLAS_JES_Eta_Modelling/Eta_ModellingJESDown/";
	  else if(sys=="Eta_ModellingJESUp") output_dir += "/ATLAS_JES_Eta_Modelling/Eta_ModellingJESUp/";
	  else if(sys=="NonClosure_AFIIJESDown") output_dir += "/ATLAS_JES_NonClosure_AFII/NonClosure_AFIIJESDown/";
	  else if(sys=="NonClosure_AFIIJESUp") output_dir += "/ATLAS_JES_NonClosure_AFII/NonClosure_AFIIJESUp/";
	  else if(sys=="Eta_StatMethodJESDown") output_dir += "/ATLAS_JES_2012_Eta_StatMethod/Eta_StatMethodJESDown/";
	  else if(sys=="Eta_StatMethodJESUp") output_dir += "/ATLAS_JES_2012_Eta_StatMethod/Eta_StatMethodJESUp/";
	  else if(sys=="NP_Detector1JESDown") output_dir += "/ATLAS_JES_2012_Detector1/NP_Detector1JESDown/";
	  else if(sys=="NP_Detector1JESUp") output_dir += "/ATLAS_JES_2012_Detector1/NP_Detector1JESUp/";
	  else if(sys=="NP_Modelling1JESDown") output_dir += "/ATLAS_JES_2012_Modelling1/NP_Modelling1JESDown/";
	  else if(sys=="NP_Modelling1JESUp") output_dir += "/ATLAS_JES_2012_Modelling1/NP_Modelling1JESUp/";
	  else if(sys=="NPVJESDown") output_dir += "/ATLAS_JES_NPV/NPVJESDown/";
	  else if(sys=="NPVJESUp") output_dir += "/ATLAS_JES_NPV/NPVJESUp/";
		else if(sys=="SoftTrackResoPara") output_dir += "/ATLAS_TRACKMET_RESOPARASOFT/SoftTrackResoPara/";
		else if(sys=="SoftTrackResoPerp") output_dir += "/ATLAS_TRACKMET_RESOPERPSOFT/SoftTrackResoPerp/";
		else if(sys=="SoftTrackScaleDown") output_dir += "/ATLAS_TRACKMET_SCALESOFT/SoftTrackScaleDown/";
		else if(sys=="SoftTrackScaleUp") output_dir += "/ATLAS_TRACKMET_SCALESOFT/SoftTrackScaleUp/";

	} else output_dir = output_dir+"/"+sys+"/";
	gSystem->Exec("mkdir -p "+output_dir);
    
	TMVA::Reader *reader_emu_0jet = new TMVA::Reader("!Color:!Silent");
	TMVA::Reader *reader_emu_1jet = new TMVA::Reader("!Color:!Silent");
	TMVA::Reader *reader_evenodd_emu_0jet(0), *reader_evenodd_emu_1jet(0);
	if(runEvenOdd) {
		reader_evenodd_emu_0jet = new TMVA::Reader("!Color:!Silent");
		reader_evenodd_emu_1jet = new TMVA::Reader("!Color:!Silent");
	}

	TMVA::Reader *reader_train2_emu_0jet(0), *reader_train2_emu_1jet(0);
  TMVA::Reader *reader_evenodd_train2_emu_0jet(0), *reader_evenodd_train2_emu_1jet(0);
	if(runSpin2training) {
		reader_train2_emu_0jet = new TMVA::Reader("!Color:!Silent");
		reader_train2_emu_1jet = new TMVA::Reader("!Color:!Silent");
		if(runEvenOdd) {
			reader_evenodd_train2_emu_0jet = new TMVA::Reader("!Color:!Silent");
			reader_evenodd_train2_emu_1jet = new TMVA::Reader("!Color:!Silent");
		}
	}

	TMVA::Reader *reader_subtrain_emu_0jet(0), *reader_subtrain_emu_1jet(0);
  TMVA::Reader *reader_evenodd_subtrain_emu_0jet(0), *reader_evenodd_subtrain_emu_1jet(0);
	if(useXjetsSubtraining) {
		reader_subtrain_emu_0jet = new TMVA::Reader("!Color:!Silent");
		reader_subtrain_emu_1jet = new TMVA::Reader("!Color:!Silent");
		if(runEvenOdd) {
			reader_evenodd_subtrain_emu_0jet = new TMVA::Reader("!Color:!Silent");
			reader_evenodd_subtrain_emu_1jet = new TMVA::Reader("!Color:!Silent");
		}
	}

	Bool_t overlapWZ, passMueORsimple;
	Int_t intoverlapWZ, intpassMueORsimple, isZll, notZgamma;
	Double_t  MVAEventWeight, WgStarEventWeight, higgsPtEventWeight; 
	Float_t MT, MT_TrackHWW_Clj, DPhill, Ptll, Mll, MET, METRel, MET_phi, jetPt0, jetEta0, jetPhi0, DEtall, MinDPhi_TrackHWW_Clj;
	Float_t METRel_TrackHWW_Clj, MET_TrackHWW_Clj, MET_phi_TrackHWW_Clj, MET_x_TrackHWW_Clj, MET_y_TrackHWW_Clj;
	Float_t lepPt0, lepPt1, lepEta0, lepEta1, lepPhi0, lepPhi1, Mtt, Mtt_TrackHWW_Clj, lepID0, lepID1, HPt, Efun, Ptll_truth;
	Float_t BDPhill, BDPsill, BPLeadLep, BPSubLeadLep, BELeadNeu, BESubLeadNeu, Esum, DPt, leadLepMT, leadLepPt, subleadLepPt, MTtruth; 
	Float_t HPttruth, higgsPt;

	Int_t m_jet_n, m_el_n, m_mu_n, isLowPtCand, nJets_Pt20_MV1_85, isBlinded;
	UInt_t RunNumber, EventNumber, mc_channel_number;
	Int_t wjets_RunNumber, wjets_EventNumber;
	Double_t PDFWeight[2], pileupEventWeight_080, pileupEventWeight, pileupEventWeight_090;
	Double_t lepTrigSFEventWeight, lepTrigSFEventWeightUp, lepTrigSFEventWeightDown, lepSF0Error, lepSF1Error;
	Double_t MV120_85_EventWeight, MV120_85_BJetWeight, MV120_85_BJetWeightUp, MV120_85_BJetWeightDown;
	Double_t MV120_85_CTJetWeight, MV120_85_CTJetWeightUp, MV120_85_CTJetWeightDown;
	Double_t MV120_85_MisTagWeight, MV120_85_MisTagWeightUp, MV120_85_MisTagWeightDown;
	Double_t lepSF0ErrorIso, lepSF1ErrorIso, lepSF0EventWeight, lepSF1EventWeight;
	//	Double_t pt_weight = 1.;
	//for wjets/qcd systematics
	Double_t SysFakeWeight = 1; 	Int_t isFake0 = 0; Int_t isFake1 = 0;
	Double_t SysFakeStat = 1; Double_t SysFakeFlav = 1; Double_t SysFakeOther = 1;
        
	Int_t                mc_n;
	std::vector<float>   *mc_pt;
	std::vector<float>   *mc_m;
	std::vector<float>   *mc_eta;
	std::vector<float>   *mc_phi;
	std::vector<int>     *mc_status;
	std::vector<int>     *mc_pdgId;
	
	std::vector<int> *m_mcevt_pdf_id1;
	std::vector<int> *m_mcevt_pdf_id2;
	std::vector<double> *m_mcevt_pdf_x1;
	std::vector<double> *m_mcevt_pdf_x2;
	std::vector<double> *m_mcevt_pdf_scale;
    
	mc_pt = 0; mc_m = 0; mc_eta = 0; mc_phi = 0; mc_status = 0; mc_pdgId = 0;
	m_mcevt_pdf_id1=0; m_mcevt_pdf_id2=0; m_mcevt_pdf_x1=0; m_mcevt_pdf_x2=0; m_mcevt_pdf_scale=0;
                            
	if(ana_name.Contains("mtmllptllmtt")) {
		std::cout << "ana_name: mtmllptllmtt" << std::endl;
		reader_emu_0jet->AddVariable("MT_TrackHWW_Clj", &MT_TrackHWW_Clj);
		reader_emu_0jet->AddVariable("Mll", &Mll);
		reader_emu_0jet->AddVariable("Ptll", &Ptll);
		reader_emu_0jet->AddVariable("DPhill", &DPhill);

		reader_emu_1jet->AddVariable("MT_TrackHWW_Clj", &MT_TrackHWW_Clj);
		reader_emu_1jet->AddVariable("Mll", &Mll);
		reader_emu_1jet->AddVariable("Ptll", &Ptll);
		reader_emu_1jet->AddVariable("max(-1,Mtt_TrackHWW_Clj)", &Mtt_TrackHWW_Clj);

		if(runEvenOdd) {
			reader_evenodd_emu_0jet->AddVariable("MT_TrackHWW_Clj", &MT_TrackHWW_Clj);
			reader_evenodd_emu_0jet->AddVariable("Mll", &Mll);
			reader_evenodd_emu_0jet->AddVariable("Ptll", &Ptll);
			reader_evenodd_emu_0jet->AddVariable("DPhill", &DPhill);

			reader_evenodd_emu_1jet->AddVariable("MT_TrackHWW_Clj", &MT_TrackHWW_Clj);
			reader_evenodd_emu_1jet->AddVariable("Mll", &Mll);
			reader_evenodd_emu_1jet->AddVariable("Ptll", &Ptll);
			reader_evenodd_emu_1jet->AddVariable("max(-1,Mtt_TrackHWW_Clj)", &Mtt_TrackHWW_Clj);
		}
	}

	else if(ana_name.Contains("CPboosted")) {
		std::cout << "ana_name: CPboosted" << std::endl;
		reader_emu_0jet->AddVariable("Esum", &Esum);
		reader_emu_0jet->AddVariable("Mll", &Mll);
		reader_emu_0jet->AddVariable("abs(DPt)", &DPt);
		reader_emu_0jet->AddVariable("DPhill", &DPhill);

		reader_emu_1jet->AddVariable("Esum", &Esum);
		reader_emu_1jet->AddVariable("Mll", &Mll);
		reader_emu_1jet->AddVariable("abs(DPt)", &DPt);
		reader_emu_1jet->AddVariable("DPhill", &DPhill);

		if(runEvenOdd) {
			reader_evenodd_emu_0jet->AddVariable("Esum", &Esum);
			reader_evenodd_emu_0jet->AddVariable("Mll", &Mll);
			reader_evenodd_emu_0jet->AddVariable("abs(DPt)", &DPt);
			reader_evenodd_emu_0jet->AddVariable("DPhill", &DPhill);

			reader_evenodd_emu_1jet->AddVariable("Esum", &Esum);
			reader_evenodd_emu_1jet->AddVariable("Mll", &Mll);
			reader_evenodd_emu_1jet->AddVariable("abs(DPt)", &DPt);
			reader_evenodd_emu_1jet->AddVariable("DPhill", &DPhill);
		}
	}

	else if(ana_name.Contains("CPEfun")) {
		std::cout << "ana_name: CPEfun" << std::endl;
		reader_emu_0jet->AddVariable("lepPt0>lepPt1?(lepPt0-0.5*lepPt1+0.5*MET_TrackHWW_Clj):(lepPt1-0.5*lepPt0+0.5*MET_TrackHWW_Clj)", &Efun);
		reader_emu_0jet->AddVariable("Mll", &Mll);
		reader_emu_0jet->AddVariable("abs(DPt)", &DPt);
		reader_emu_0jet->AddVariable("DPhill", &DPhill);

		reader_emu_1jet->AddVariable("lepPt0>lepPt1?(lepPt0-0.5*lepPt1+0.5*MET_TrackHWW_Clj):(lepPt1-0.5*lepPt0+0.5*MET_TrackHWW_Clj)", &Efun);
		reader_emu_1jet->AddVariable("Mll", &Mll);
		reader_emu_1jet->AddVariable("abs(DPt)", &DPt);
		reader_emu_1jet->AddVariable("DPhill", &DPhill);

		if(runEvenOdd) {
			reader_evenodd_emu_0jet->AddVariable("lepPt0>lepPt1?(lepPt0-0.5*lepPt1+0.5*MET_TrackHWW_Clj):(lepPt1-0.5*lepPt0+0.5*MET_TrackHWW_Clj)", &Efun);
			reader_evenodd_emu_0jet->AddVariable("Mll", &Mll);
			reader_evenodd_emu_0jet->AddVariable("abs(DPt)", &DPt);
			reader_evenodd_emu_0jet->AddVariable("DPhill", &DPhill);

			reader_evenodd_emu_1jet->AddVariable("lepPt0>lepPt1?(lepPt0-0.5*lepPt1+0.5*MET_TrackHWW_Clj):(lepPt1-0.5*lepPt0+0.5*MET_TrackHWW_Clj)", &Efun);
			reader_evenodd_emu_1jet->AddVariable("Mll", &Mll);
			reader_evenodd_emu_1jet->AddVariable("abs(DPt)", &DPt);
			reader_evenodd_emu_1jet->AddVariable("DPhill", &DPhill);
		}
	}

	else if(ana_name.Contains("4var")) {
		std::cout << "ana_name: 4var" << std::endl;
		reader_emu_0jet->AddVariable("MT_TrackHWW_Clj", &MT_TrackHWW_Clj);
		reader_emu_0jet->AddVariable("Mll", &Mll);
		reader_emu_0jet->AddVariable("Ptll", &Ptll);
		reader_emu_0jet->AddVariable("DPhill", &DPhill);

		reader_emu_1jet->AddVariable("MT_TrackHWW_Clj", &MT_TrackHWW_Clj);
		reader_emu_1jet->AddVariable("Mll", &Mll);
		reader_emu_1jet->AddVariable("DPhill", &DPhill);
		reader_emu_1jet->AddVariable("Ptll", &Ptll);


		if(runEvenOdd) {
			reader_evenodd_emu_0jet->AddVariable("MT_TrackHWW_Clj", &MT_TrackHWW_Clj);
			reader_evenodd_emu_0jet->AddVariable("Mll", &Mll);
			reader_evenodd_emu_0jet->AddVariable("Ptll", &Ptll);
			reader_evenodd_emu_0jet->AddVariable("DPhill", &DPhill);

			reader_evenodd_emu_1jet->AddVariable("MT_TrackHWW_Clj", &MT_TrackHWW_Clj);
			reader_evenodd_emu_1jet->AddVariable("Mll", &Mll);
			reader_evenodd_emu_1jet->AddVariable("DPhill", &DPhill);
			reader_evenodd_emu_1jet->AddVariable("Ptll", &Ptll);

		}
	}	else if(ana_name.Contains("kHWWScan")) {
		std::cout << "ana_name: kHWWScan" << std::endl;
		reader_emu_0jet->AddVariable("MET_TrackHWW_Clj", &MET_TrackHWW_Clj);
		reader_emu_0jet->AddVariable("Mll", &Mll);
		reader_emu_0jet->AddVariable("Ptll", &Ptll);
		reader_emu_0jet->AddVariable("DPhill", &DPhill);

		reader_emu_1jet->AddVariable("MET_TrackHWW_Clj", &MET_TrackHWW_Clj);
		reader_emu_1jet->AddVariable("Mll", &Mll);
		reader_emu_1jet->AddVariable("Ptll", &Ptll);
		reader_emu_1jet->AddVariable("DPhill", &DPhill);



		if(runEvenOdd) {
			reader_evenodd_emu_0jet->AddVariable("MET_TrackHWW_Clj", &MET_TrackHWW_Clj);
			reader_evenodd_emu_0jet->AddVariable("Mll", &Mll);
			reader_evenodd_emu_0jet->AddVariable("Ptll", &Ptll);
			reader_evenodd_emu_0jet->AddVariable("DPhill", &DPhill);

			reader_evenodd_emu_1jet->AddVariable("MET_TrackHWW_Clj", &MET_TrackHWW_Clj);
			reader_evenodd_emu_1jet->AddVariable("Mll", &Mll);
			reader_evenodd_emu_1jet->AddVariable("Ptll", &Ptll);
			reader_evenodd_emu_1jet->AddVariable("DPhill", &DPhill);


		}
	} else if(ana_name.Contains("sigVsBgr")) {
		std::cout << "ana_name: sigVsBgr" << std::endl;
		reader_emu_0jet->AddVariable("Esum", &Esum);
		reader_emu_0jet->AddVariable("Mll", &Mll);
		reader_emu_0jet->AddVariable("MT_TrackHWW_Clj", &MT_TrackHWW_Clj);
		reader_emu_0jet->AddVariable("DPhill", &DPhill);

		reader_emu_1jet->AddVariable("Esum", &Esum);
		reader_emu_1jet->AddVariable("Mll", &Mll);
		reader_emu_1jet->AddVariable("MT_TrackHWW_Clj", &MT_TrackHWW_Clj);
		reader_emu_1jet->AddVariable("DPhill", &DPhill);

		if(runEvenOdd) {
			reader_evenodd_emu_0jet->AddVariable("Esum", &Esum);
			reader_evenodd_emu_0jet->AddVariable("Mll", &Mll);
			reader_evenodd_emu_0jet->AddVariable("MT_TrackHWW_Clj", &MT_TrackHWW_Clj);
			reader_evenodd_emu_0jet->AddVariable("DPhill", &DPhill);

			reader_evenodd_emu_1jet->AddVariable("Esum", &Esum);
			reader_evenodd_emu_1jet->AddVariable("Mll", &Mll);
			reader_evenodd_emu_1jet->AddVariable("MT_TrackHWW_Clj", &MT_TrackHWW_Clj);
			reader_evenodd_emu_1jet->AddVariable("DPhill", &DPhill);


		}
	}

	if(runSpin2training){
		if(ana2_name.Contains("mtmllptllmtt")) {
		std::cout << "ana2_name: mtmllptllmtt" << std::endl;
			reader_train2_emu_0jet->AddVariable("MT_TrackHWW_Clj", &MT_TrackHWW_Clj);
			reader_train2_emu_0jet->AddVariable("Mll", &Mll);
			reader_train2_emu_0jet->AddVariable("Ptll", &Ptll);
			reader_train2_emu_0jet->AddVariable("DPhill", &DPhill);

			reader_train2_emu_1jet->AddVariable("MT_TrackHWW_Clj", &MT_TrackHWW_Clj);
			reader_train2_emu_1jet->AddVariable("Mll", &Mll);
			reader_train2_emu_1jet->AddVariable("Ptll", &Ptll);
			reader_train2_emu_1jet->AddVariable("max(-1,Mtt_TrackHWW_Clj)", &Mtt_TrackHWW_Clj);

			if(runEvenOdd) {
				reader_evenodd_train2_emu_0jet->AddVariable("MT_TrackHWW_Clj", &MT_TrackHWW_Clj);
				reader_evenodd_train2_emu_0jet->AddVariable("Mll", &Mll);
				reader_evenodd_train2_emu_0jet->AddVariable("Ptll", &Ptll);
				reader_evenodd_train2_emu_0jet->AddVariable("DPhill", &DPhill);

				reader_evenodd_train2_emu_1jet->AddVariable("MT_TrackHWW_Clj", &MT_TrackHWW_Clj);
				reader_evenodd_train2_emu_1jet->AddVariable("Mll", &Mll);
				reader_evenodd_train2_emu_1jet->AddVariable("Ptll", &Ptll);
				reader_evenodd_train2_emu_1jet->AddVariable("max(-1,Mtt_TrackHWW_Clj)", &Mtt_TrackHWW_Clj);
			}
		}
		else if(ana2_name.Contains("CPboosted")) {
			std::cout << "ana2_name: CPboosted" << std::endl;
			reader_train2_emu_0jet->AddVariable("Esum", &Esum);
			reader_train2_emu_0jet->AddVariable("Mll", &Mll);
			reader_train2_emu_0jet->AddVariable("abs(DPt)", &DPt);
			reader_train2_emu_0jet->AddVariable("DPhill", &DPhill);

			reader_train2_emu_1jet->AddVariable("Esum", &Esum);
			reader_train2_emu_1jet->AddVariable("Mll", &Mll);
			reader_train2_emu_1jet->AddVariable("abs(DPt)", &DPt);
			reader_train2_emu_1jet->AddVariable("DPhill", &DPhill);

			if(runEvenOdd) {
				reader_evenodd_train2_emu_0jet->AddVariable("Esum", &Esum);
				reader_evenodd_train2_emu_0jet->AddVariable("Mll", &Mll);
				reader_evenodd_train2_emu_0jet->AddVariable("abs(DPt)", &DPt);
				reader_evenodd_train2_emu_0jet->AddVariable("DPhill", &DPhill);

				reader_evenodd_train2_emu_1jet->AddVariable("Esum", &Esum);
				reader_evenodd_train2_emu_1jet->AddVariable("Mll", &Mll);
				reader_evenodd_train2_emu_1jet->AddVariable("abs(DPt)", &DPt);
				reader_evenodd_train2_emu_1jet->AddVariable("DPhill", &DPhill);
			}
		}//CPboosted

		else if(ana2_name.Contains("CPEfun")) {
		std::cout << "ana2_name: CPEfun" << std::endl;
			reader_train2_emu_0jet->AddVariable("lepPt0>lepPt1?(lepPt0-0.5*lepPt1+0.5*MET_TrackHWW_Clj):(lepPt1-0.5*lepPt0+0.5*MET_TrackHWW_Clj)", &Efun);
			reader_train2_emu_0jet->AddVariable("Mll", &Mll);
			reader_train2_emu_0jet->AddVariable("abs(DPt)", &DPt);
			reader_train2_emu_0jet->AddVariable("DPhill", &DPhill);

			reader_train2_emu_1jet->AddVariable("lepPt0>lepPt1?(lepPt0-0.5*lepPt1+0.5*MET_TrackHWW_Clj):(lepPt1-0.5*lepPt0+0.5*MET_TrackHWW_Clj)", &Efun);
			reader_train2_emu_1jet->AddVariable("Mll", &Mll);
			reader_train2_emu_1jet->AddVariable("abs(DPt)", &DPt);
			reader_train2_emu_1jet->AddVariable("DPhill", &DPhill);

			if(runEvenOdd) {
				reader_evenodd_train2_emu_0jet->AddVariable("lepPt0>lepPt1?(lepPt0-0.5*lepPt1+0.5*MET_TrackHWW_Clj):(lepPt1-0.5*lepPt0+0.5*MET_TrackHWW_Clj)", &Efun);
				reader_evenodd_train2_emu_0jet->AddVariable("Mll", &Mll);
				reader_evenodd_train2_emu_0jet->AddVariable("abs(DPt)", &DPt);
				reader_evenodd_train2_emu_0jet->AddVariable("DPhill", &DPhill);

				reader_evenodd_train2_emu_1jet->AddVariable("lepPt0>lepPt1?(lepPt0-0.5*lepPt1+0.5*MET_TrackHWW_Clj):(lepPt1-0.5*lepPt0+0.5*MET_TrackHWW_Clj)", &Efun);
				reader_evenodd_train2_emu_1jet->AddVariable("Mll", &Mll);
				reader_evenodd_train2_emu_1jet->AddVariable("abs(DPt)", &DPt);
				reader_evenodd_train2_emu_1jet->AddVariable("DPhill", &DPhill);
			}
		}//CPEfun

		else if(ana2_name.Contains("4var")) {
			std::cout << "ana2_name: 4var" << std::endl;
			reader_train2_emu_0jet->AddVariable("MT_TrackHWW_Clj", &MT_TrackHWW_Clj);
			reader_train2_emu_0jet->AddVariable("Mll", &Mll);
			reader_train2_emu_0jet->AddVariable("Ptll", &Ptll);
			reader_train2_emu_0jet->AddVariable("DPhill", &DPhill);

			reader_train2_emu_1jet->AddVariable("MT_TrackHWW_Clj", &MT_TrackHWW_Clj);
			reader_train2_emu_1jet->AddVariable("Mll", &Mll);
			reader_train2_emu_1jet->AddVariable("DPhill", &DPhill);
			reader_train2_emu_1jet->AddVariable("Ptll", &Ptll);


			if(runEvenOdd) {
				reader_evenodd_train2_emu_0jet->AddVariable("MT_TrackHWW_Clj", &MT_TrackHWW_Clj);
				reader_evenodd_train2_emu_0jet->AddVariable("Mll", &Mll);
				reader_evenodd_train2_emu_0jet->AddVariable("Ptll", &Ptll);
				reader_evenodd_train2_emu_0jet->AddVariable("DPhill", &DPhill);

				reader_evenodd_train2_emu_1jet->AddVariable("MT_TrackHWW_Clj", &MT_TrackHWW_Clj);
				reader_evenodd_train2_emu_1jet->AddVariable("Mll", &Mll);
				reader_evenodd_train2_emu_1jet->AddVariable("DPhill", &DPhill);
				reader_evenodd_train2_emu_1jet->AddVariable("Ptll", &Ptll);


			}	
		} else if(ana2_name.Contains("kHWWScan")) {
			std::cout << "ana2_name: kHWWScan" << std::endl;
			reader_train2_emu_0jet->AddVariable("MET_TrackHWW_Clj", &MET_TrackHWW_Clj);
			reader_train2_emu_0jet->AddVariable("Mll", &Mll);
			reader_train2_emu_0jet->AddVariable("Ptll", &Ptll);
			reader_train2_emu_0jet->AddVariable("DPhill", &DPhill);

			reader_train2_emu_1jet->AddVariable("MET_TrackHWW_Clj", &MET_TrackHWW_Clj);
			reader_train2_emu_1jet->AddVariable("Mll", &Mll);
			reader_train2_emu_1jet->AddVariable("Ptll", &Ptll);
			reader_train2_emu_1jet->AddVariable("DPhill", &DPhill);



			if(runEvenOdd) {
				reader_evenodd_train2_emu_0jet->AddVariable("MET_TrackHWW_Clj", &MET_TrackHWW_Clj);
				reader_evenodd_train2_emu_0jet->AddVariable("Mll", &Mll);
				reader_evenodd_train2_emu_0jet->AddVariable("Ptll", &Ptll);
				reader_evenodd_train2_emu_0jet->AddVariable("DPhill", &DPhill);

				reader_evenodd_train2_emu_1jet->AddVariable("MET_TrackHWW_Clj", &MET_TrackHWW_Clj);
				reader_evenodd_train2_emu_1jet->AddVariable("Mll", &Mll);
				reader_evenodd_train2_emu_1jet->AddVariable("Ptll", &Ptll);
				reader_evenodd_train2_emu_1jet->AddVariable("DPhill", &DPhill);


			}
		} else if(ana2_name.Contains("sigVsBgr")) {
			std::cout << "ana2_name: sigVsBrg" << std::endl;
			reader_train2_emu_0jet->AddVariable("Esum", &Esum);
			reader_train2_emu_0jet->AddVariable("Mll", &Mll);
			reader_train2_emu_0jet->AddVariable("MT_TrackHWW_Clj", &MT_TrackHWW_Clj);
			reader_train2_emu_0jet->AddVariable("DPhill", &DPhill);

			reader_train2_emu_1jet->AddVariable("Esum", &Esum);
			reader_train2_emu_1jet->AddVariable("Mll", &Mll);
			reader_train2_emu_1jet->AddVariable("MT_TrackHWW_Clj", &MT_TrackHWW_Clj);
			reader_train2_emu_1jet->AddVariable("DPhill", &DPhill);

			if(runEvenOdd) {
				reader_evenodd_train2_emu_0jet->AddVariable("Esum", &Esum);
				reader_evenodd_train2_emu_0jet->AddVariable("Mll", &Mll);
				reader_evenodd_train2_emu_0jet->AddVariable("MT_TrackHWW_Clj", &MT_TrackHWW_Clj);
				reader_evenodd_train2_emu_0jet->AddVariable("DPhill", &DPhill);

				reader_evenodd_train2_emu_1jet->AddVariable("Esum", &Esum);
				reader_evenodd_train2_emu_1jet->AddVariable("Mll", &Mll);
				reader_evenodd_train2_emu_1jet->AddVariable("MT_TrackHWW_Clj", &MT_TrackHWW_Clj);
				reader_evenodd_train2_emu_1jet->AddVariable("DPhill", &DPhill);


			}
		}

	}//if runSpin2training

	if (useXjetsSubtraining) {
		if(sub_ana_name.Contains("wjets")) {
			std::cout << "sub_ana_name: wjets" << std::endl;
			reader_subtrain_emu_0jet->AddVariable("MinDPhi_TrackHWW_Clj", &MinDPhi_TrackHWW_Clj);
			reader_subtrain_emu_0jet->AddVariable("DEtall", &DEtall);
			reader_subtrain_emu_0jet->AddVariable("leadLepPt", &leadLepPt);
			reader_subtrain_emu_0jet->AddVariable("MT_TrackHWW_Clj", &MT_TrackHWW_Clj);
			reader_subtrain_emu_0jet->AddVariable("subleadLepPt", &subleadLepPt);
			reader_subtrain_emu_0jet->AddVariable("leadLepMT", &leadLepMT);
			reader_subtrain_emu_0jet->AddVariable("MET_TrackHWW_Clj", &MET_TrackHWW_Clj);
			reader_subtrain_emu_0jet->AddVariable("Mll", &Mll);

			reader_subtrain_emu_1jet->AddVariable("MinDPhi_TrackHWW_Clj", &MinDPhi_TrackHWW_Clj);
			reader_subtrain_emu_1jet->AddVariable("DEtall", &DEtall);
			reader_subtrain_emu_1jet->AddVariable("leadLepPt", &leadLepPt);
			reader_subtrain_emu_1jet->AddVariable("MT_TrackHWW_Clj", &MT_TrackHWW_Clj);
			reader_subtrain_emu_1jet->AddVariable("subleadLepPt", &subleadLepPt);
			reader_subtrain_emu_1jet->AddVariable("leadLepMT", &leadLepMT);
			reader_subtrain_emu_1jet->AddVariable("MET_TrackHWW_Clj", &MET_TrackHWW_Clj);
			reader_subtrain_emu_1jet->AddVariable("Mll", &Mll);

			if(runEvenOdd) {
				reader_evenodd_subtrain_emu_0jet->AddVariable("MT_TrackHWW_Clj", &MT_TrackHWW_Clj);
				reader_evenodd_subtrain_emu_0jet->AddVariable("DEtall", &DEtall);
				reader_evenodd_subtrain_emu_0jet->AddVariable("MinDPhi_TrackHWW_Clj", &MinDPhi_TrackHWW_Clj);
				reader_evenodd_subtrain_emu_0jet->AddVariable("leadLepPt", &leadLepPt);
				reader_evenodd_subtrain_emu_0jet->AddVariable("subleadLepPt", &subleadLepPt);
				reader_evenodd_subtrain_emu_0jet->AddVariable("leadLepMT", &leadLepMT);
				reader_evenodd_subtrain_emu_0jet->AddVariable("MET_TrackHWW_Clj", &MET_TrackHWW_Clj);
				reader_evenodd_subtrain_emu_0jet->AddVariable("Mll", &Mll);

				reader_evenodd_subtrain_emu_1jet->AddVariable("MT_TrackHWW_Clj", &MT_TrackHWW_Clj);
				reader_evenodd_subtrain_emu_1jet->AddVariable("DEtall", &DEtall);
				reader_evenodd_subtrain_emu_1jet->AddVariable("MinDPhi_TrackHWW_Clj", &MinDPhi_TrackHWW_Clj);
				reader_evenodd_subtrain_emu_1jet->AddVariable("leadLepPt", &leadLepPt);
				reader_evenodd_subtrain_emu_1jet->AddVariable("subleadLepPt", &subleadLepPt);
				reader_evenodd_subtrain_emu_1jet->AddVariable("leadLepMT", &leadLepMT);
				reader_evenodd_subtrain_emu_1jet->AddVariable("MET_TrackHWW_Clj", &MET_TrackHWW_Clj);
				reader_evenodd_subtrain_emu_1jet->AddVariable("Mll", &Mll);
			}
		}
	}// if useXjetsSubtraining
	
	//BookMVA
	reader_emu_0jet->BookMVA(analysis, weightsdir+traineven+"H"+HiggsMass+"_0jet_emu/TMVAClassification_"+ana_name+"_"+traineven+"H"+HiggsMass+"_0jet_emu_"+analysis+".weights.xml");
	reader_emu_1jet->BookMVA(analysis, weightsdir+traineven+"H"+HiggsMass+"_1jet_emu/TMVAClassification_"+ana_name+"_"+traineven+"H"+HiggsMass+"_1jet_emu_"+analysis+".weights.xml");
	if(runEvenOdd){
		reader_evenodd_emu_0jet->BookMVA(analysis, weightsdir+trainodd+"H"+HiggsMass+"_0jet_emu/TMVAClassification_"+ana_name+"_"+trainodd+"H"+HiggsMass+"_0jet_emu_"+analysis+".weights.xml");
		reader_evenodd_emu_1jet->BookMVA(analysis, weightsdir+trainodd+"H"+HiggsMass+"_1jet_emu/TMVAClassification_"+ana_name+"_"+trainodd+"H"+HiggsMass+"_1jet_emu_"+analysis+".weights.xml");
	}

	if(runSpin2training){
		reader_train2_emu_0jet->BookMVA(analysis, weightsdir2+traineven+"H"+HiggsMass+"_0jet_emu/TMVAClassification_"+ana2_name+"_"+traineven+"H"+HiggsMass+"_0jet_emu_"+analysis+".weights.xml");
		reader_train2_emu_1jet->BookMVA(analysis, weightsdir2+traineven+"H"+HiggsMass+"_1jet_emu/TMVAClassification_"+ana2_name+"_"+traineven+"H"+HiggsMass+"_1jet_emu_"+analysis+".weights.xml");
		if(runEvenOdd){
			reader_evenodd_train2_emu_0jet->BookMVA(analysis, weightsdir2+trainodd+"H"+HiggsMass+"_0jet_emu/TMVAClassification_"+ana2_name+"_"+trainodd+"H"+HiggsMass+"_0jet_emu_"+analysis+".weights.xml");
			reader_evenodd_train2_emu_1jet->BookMVA(analysis, weightsdir2+trainodd+"H"+HiggsMass+"_1jet_emu/TMVAClassification_"+ana2_name+"_"+trainodd+"H"+HiggsMass+"_1jet_emu_"+analysis+".weights.xml");
		}
	}

	if(useXjetsSubtraining){
		reader_subtrain_emu_0jet->BookMVA("WjetsBDT_8var_cfg_242_trainOnEven", sub_weightsdir+traineven+"H"+HiggsMass+"_0jet_emu/TMVAClassification_"+sub_ana_name+"_"+traineven+"H"+HiggsMass+"_0jet_emu.weights.xml");
		reader_subtrain_emu_1jet->BookMVA("WjetsBDT_8var_cfg_242_trainOnEven", sub_weightsdir+traineven+"H"+HiggsMass+"_1jet_emu/TMVAClassification_"+sub_ana_name+"_"+traineven+"H"+HiggsMass+"_1jet_emu.weights.xml");
		if(runEvenOdd){
			reader_evenodd_subtrain_emu_0jet->BookMVA("WjetsBDT_8var_cfg_242_trainOnOdd", sub_weightsdir+trainodd+"H"+HiggsMass+"_0jet_emu/TMVAClassification_"+sub_ana_name+"_"+trainodd+"H"+HiggsMass+"_0jet_emu.weights.xml");
			reader_evenodd_subtrain_emu_1jet->BookMVA("WjetsBDT_8var_cfg_242_trainOnOdd", sub_weightsdir+trainodd+"H"+HiggsMass+"_1jet_emu/TMVAClassification_"+sub_ana_name+"_"+trainodd+"H"+HiggsMass+"_1jet_emu.weights.xml");
		}
	}


	Int_t samp_i = 21; sampN=27;
	Bool_t runoverallsamples = false;
	if(doTopDataDrivenNtuples) samp_i = 3;
	if(sys=="EL_SysFakeUp_MVA" || sys=="EL_SysFakeDown_MVA" || sys=="MU_SysFakeUp_MVA" || sys=="MU_SysFakeDown_MVA") samp_i = 27;
	if(sys=="NonClosure_AFIIJESDown" || sys=="NonClosure_AFIIJESUp") {samp_i = 1; sampN = 27;}
	if(samp_i==0 && sampN>20) runoverallsamples = true;

	for(; samp_i < sampN; samp_i++){
	  
		//get input files
		TString fname = input_dir + fnames[samp_i] + "_"+ sys +".root";
		if ( TString(fnames[samp_i]).Contains("data")) fname = input_dir + fnames[samp_i] + ".root";

		std::cout << "##### sample = " << fname << std::endl;
		TFile *f_input = TFile::Open(fname);
		if(!f_input) continue;
		TString treename = "";
		if( TString(fnames[samp_i]).Contains("Wjets_datadriven") ) treename="WjetsCR"; //WjetsCR_dil
		else if ( TString(fnames[samp_i]).Contains("QCD_datadriven") ) treename="QCDCR";
		else  treename="HWWTree";
		TTree *t_input = (TTree*)f_input->Get(treename);
	  
		if(TString(fnames[samp_i]).Contains("Zjets") || TString(fnames[samp_i]).Contains("DrellYan")) {
			t_input->SetBranchAddress("PDFWeight", PDFWeight);
			t_input->SetBranchAddress("Ptll_truth", &Ptll_truth);
		}
		if(!(TString(fnames[samp_i]).Contains("Wjets_datadriven") || TString(fnames[samp_i]).Contains("QCD_datadriven")) ) {

			t_input->SetBranchAddress("pileupEventWeight", &pileupEventWeight);
			t_input->SetBranchAddress("pileupEventWeight_090", &pileupEventWeight_090);
			t_input->SetBranchAddress("pileupEventWeight_080", &pileupEventWeight_080);
			t_input->SetBranchAddress("lepTrigSFEventWeight", &lepTrigSFEventWeight);
			t_input->SetBranchAddress("MV120_85_EventWeight", &MV120_85_EventWeight);
			t_input->SetBranchAddress("MV120_85_CTJetWeight", &MV120_85_CTJetWeight);
			t_input->SetBranchAddress("MV120_85_BJetWeight", &MV120_85_BJetWeight);
			t_input->SetBranchAddress("MV120_85_MisTagWeight", &MV120_85_MisTagWeight);
			t_input->SetBranchAddress("lepSF0EventWeight", &lepSF0EventWeight);
			t_input->SetBranchAddress("lepSF1EventWeight", &lepSF1EventWeight);
			t_input->SetBranchAddress("RunNumber", &RunNumber);
			t_input->SetBranchAddress("EventNumber", &EventNumber);
			t_input->SetBranchAddress("mc_channel_number", &mc_channel_number);
			t_input->SetBranchAddress("passMueORsimple", &passMueORsimple);
			t_input->SetBranchAddress("higgsPtEventWeight", &higgsPtEventWeight);
			if(sys=="Nominal"){
				t_input->SetBranchAddress("lepTrigSFEventWeightUp", &lepTrigSFEventWeightUp);
				t_input->SetBranchAddress("lepTrigSFEventWeightDown", &lepTrigSFEventWeightDown);
				t_input->SetBranchAddress("MV120_85_CTJetWeightUp", &MV120_85_CTJetWeightUp);
				t_input->SetBranchAddress("MV120_85_BJetWeightUp", &MV120_85_BJetWeightUp);
				t_input->SetBranchAddress("MV120_85_MisTagWeightUp", &MV120_85_MisTagWeightUp);
				t_input->SetBranchAddress("MV120_85_CTJetWeightDown", &MV120_85_CTJetWeightDown);
				t_input->SetBranchAddress("MV120_85_BJetWeightDown", &MV120_85_BJetWeightDown);
				t_input->SetBranchAddress("MV120_85_MisTagWeightDown", &MV120_85_MisTagWeightDown);
				t_input->SetBranchAddress("lepSF0Error", &lepSF0Error);
				t_input->SetBranchAddress("lepSF1Error", &lepSF1Error);
				t_input->SetBranchAddress("lepSF0ErrorIso", &lepSF0ErrorIso);
				t_input->SetBranchAddress("lepSF1ErrorIso", &lepSF1ErrorIso);
			}
		} else {
			t_input->SetBranchAddress("RunNumber", &wjets_RunNumber);
			t_input->SetBranchAddress("EventNumber", &wjets_EventNumber);
		} 
		t_input->SetBranchAddress("isLowPtCand", &isLowPtCand);
		t_input->SetBranchAddress("isZll", &isZll);
		t_input->SetBranchAddress("notZgamma", &notZgamma);
		t_input->SetBranchAddress("Mll", &Mll);
		t_input->SetBranchAddress("overlapWZ", &overlapWZ);
		t_input->SetBranchAddress("WgStarEventWeight", &WgStarEventWeight);
		t_input->SetBranchAddress("MT", &MT);
		t_input->SetBranchAddress("MT_TrackHWW_Clj", &MT_TrackHWW_Clj);
		t_input->SetBranchAddress("Ptll", &Ptll);
		t_input->SetBranchAddress("DPhill", &DPhill);
		t_input->SetBranchAddress("METRel", &METRel);
		t_input->SetBranchAddress("METRel_TrackHWW_Clj", &METRel_TrackHWW_Clj);
		t_input->SetBranchAddress("MET_TrackHWW_Clj", &MET_TrackHWW_Clj);
		t_input->SetBranchAddress("MET_phi_TrackHWW_Clj", &MET_phi_TrackHWW_Clj);
		t_input->SetBranchAddress("MET_x_TrackHWW_Clj", &MET_x_TrackHWW_Clj);
		t_input->SetBranchAddress("MET_y_TrackHWW_Clj", &MET_y_TrackHWW_Clj);
		t_input->SetBranchAddress("MET", &MET);
		t_input->SetBranchAddress("MET_phi", &MET_phi);
		t_input->SetBranchAddress("lepPt0", &lepPt0);
		t_input->SetBranchAddress("lepPt1", &lepPt1);
		t_input->SetBranchAddress("lepEta0", &lepEta0);
		t_input->SetBranchAddress("lepEta1", &lepEta1);
		t_input->SetBranchAddress("lepPhi0", &lepPhi0);
		t_input->SetBranchAddress("lepPhi1", &lepPhi1);
		t_input->SetBranchAddress("lepID0", &lepID0);
		t_input->SetBranchAddress("lepID1", &lepID1);
		t_input->SetBranchAddress("jetPt0", &jetPt0);
		t_input->SetBranchAddress("jetEta0", &jetEta0);
		t_input->SetBranchAddress("jetPhi0", &jetPhi0);
		t_input->SetBranchAddress("nJets_Pt20_MV1_85", &nJets_Pt20_MV1_85);
		t_input->SetBranchAddress("Mtt", &Mtt);
		t_input->SetBranchAddress("Mtt_TrackHWW_Clj", &Mtt_TrackHWW_Clj);
		t_input->SetBranchAddress("m_el_n", &m_el_n);
		t_input->SetBranchAddress("m_mu_n", &m_mu_n);
		t_input->SetBranchAddress("m_jet_n", &m_jet_n);
		t_input->SetBranchAddress("MVAEventWeight", &MVAEventWeight);
		t_input->SetBranchAddress("isBlinded", &isBlinded);
		t_input->SetBranchAddress("DPt", &DPt);
		t_input->SetBranchAddress("HPt", &HPt);
		t_input->SetBranchAddress("MinDPhi_TrackHWW_Clj", &MinDPhi_TrackHWW_Clj);
		t_input->SetBranchAddress("DEtall", &DEtall);

		//boosted vars
		t_input->SetBranchAddress("Esum",&Esum);
		t_input->SetBranchAddress("BDPhill",&BDPhill);
		t_input->SetBranchAddress("BDPsill",&BDPsill);
		t_input->SetBranchAddress("BPLeadLep",&BPLeadLep);
		t_input->SetBranchAddress("BPSubLeadLep",&BPSubLeadLep);
		t_input->SetBranchAddress("BELeadNeu",&BELeadNeu);
		t_input->SetBranchAddress("BESubLeadNeu",&BESubLeadNeu);
		//	t_input->SetBranchAddress("higgs_pt_weight", &pt_weight);

		if( !(TString(fnames[samp_i]).Contains("data") ) ) {
			t_input->SetBranchAddress("m_mcevt_pdf_id1", &m_mcevt_pdf_id1);
			t_input->SetBranchAddress("m_mcevt_pdf_id2", &m_mcevt_pdf_id2);
			t_input->SetBranchAddress("m_mcevt_pdf_x1", &m_mcevt_pdf_x1);
			t_input->SetBranchAddress("m_mcevt_pdf_x2", &m_mcevt_pdf_x2);
			t_input->SetBranchAddress("m_mcevt_pdf_scale", &m_mcevt_pdf_scale);
			t_input->SetBranchAddress("MTtruth", &MTtruth);
		}

		if(TString(fnames[samp_i]).Contains("Wjets_datadriven") || TString(fnames[samp_i]).Contains("QCD_datadriven") ){
			t_input->SetBranchAddress("SysFakeWeightCorr", &SysFakeWeight);
			t_input->SetBranchAddress("SysFakeStatCorr", &SysFakeStat);
			t_input->SetBranchAddress("SysFakeFlavCorr", &SysFakeFlav);
			t_input->SetBranchAddress("SysFakeOtherCorr", &SysFakeOther);
			t_input->SetBranchAddress("isFake0", &isFake0);
			t_input->SetBranchAddress("isFake1", &isFake1);
			t_input->SetBranchAddress("passMueORsimple", &intpassMueORsimple);
		}

		if( (TString(fnames[samp_i]).Contains("CPMixing_mg5") ) ) {
			t_input->SetBranchAddress("mc_n", &mc_n);
			t_input->SetBranchAddress("mc_pt", &mc_pt);
			t_input->SetBranchAddress("mc_m", &mc_m);
			t_input->SetBranchAddress("mc_eta", &mc_eta);
			t_input->SetBranchAddress("mc_phi", &mc_phi);
			t_input->SetBranchAddress("mc_status", &mc_status);
			t_input->SetBranchAddress("mc_pdgId", &mc_pdgId);
		}

    if(sys=="Nominal" && TString(fnames[samp_i]).Contains("ggH125_EF")) t_input->SetBranchAddress("higgsPt",&higgsPt);   

		TString outputfilename = output_dir+hnames[samp_i]+".root";
		if(!doNtuplesForLimit) outputfilename = output_dir+fnames[samp_i]+".root";
		TFile *target;
		target  = new TFile( outputfilename, "RECREATE");
		std::cout << "output file = " << outputfilename << std::endl;
        
		TTree *newtree = new TTree("MVATree","mva output tree");
                
		Double_t nMVAEventWeight, eventweight, mva_weight, npileupEventWeight, npileupEventWeight_090, npileupEventWeight_080;
		double mcevt_pdf_x1, mcevt_pdf_x2, mcevt_pdf_scale;
		int mcevt_pdf_id1, mcevt_pdf_id2;
		Double_t nlepTrigSFEventWeight, nlepTrigSFEventWeightUp, nlepTrigSFEventWeightDown, nlepSF0Error, nlepSF1Error, nlepSF0ErrorIso, nlepSF1ErrorIso;
		Double_t nlepSF0EventWeight, nlepSF1EventWeight, nMV120_85_EventWeight, nMV120_85_CTJetWeight, nMV120_85_MisTagWeight,  nMV120_85_BJetWeight;
		Double_t nMV120_85_CTJetWeightUp, nMV120_85_CTJetWeightDown, nMV120_85_MisTagWeightUp, nMV120_85_MisTagWeightDown, nMV120_85_BJetWeightUp, nMV120_85_BJetWeightDown;

		Float_t nMT, nMT_TrackHWW_Clj, nDPhill, nPtll, nMll, nMETRel, nMET, nMET_phi, nDEtall, nMinDPhi_TrackHWW_Clj;
		Float_t nMETRel_TrackHWW_Clj, nMET_TrackHWW_Clj, nMET_phi_TrackHWW_Clj, nMW0_TrackHWW_Clj, nMW1_TrackHWW_Clj;
		Float_t nMtt, nMtt_TrackHWW_Clj, njetPt0, njetEta0, njetPhi0, nHPt, nMTtruth;
		Float_t nlepID0, nlepID1, nlepPt0, nlepPt1, nlepEta0, nlepEta1, nlepPhi0, nlepPhi1;
		Float_t nBDPhill, nBDPsill, nBPLeadLep, nBPSubLeadLep, nBELeadNeu, nBESubLeadNeu, nDPt, nEsum, nhiggsPt;
		Int_t njet_n, nel_n, nmu_n, nBJets, runnumber, nEventNumber, nmc_channel_number, nisBlinded;
		Double_t xjets_mva_weight = -999.;
		Double_t spin2_mva_weight = -999.;
		Double_t spin2_xjets_mva_weight = -999.;
		//for wjets systematics
		Double_t nSysFakeWeight; Int_t nisFake0; Int_t nisFake1;
		Double_t nSysFakeStat; 		Double_t nSysFakeFlav; 		Double_t nSysFakeOther;

		Double_t ca00_ME, ca01p_ME, ca02p_ME, ca03p_ME, ca04p_ME, ca05p_ME, source_ME;
		ca00_ME = 1.0; ca01p_ME = 1.0; ca02p_ME = 1.0; ca03p_ME = 1.0; ca04p_ME = 1.0; ca05p_ME = 1.0; source_ME = 1.0;
		Double_t kHWW1p_ME, kHWW2p_ME, kHWW3p_ME, kHWW4p_ME, kHWW5p_ME, kHWW6p_ME;
		kHWW1p_ME = 1.0; kHWW2p_ME = 1.0; kHWW3p_ME = 1.0; kHWW4p_ME = 1.0; kHWW5p_ME = 1.0; kHWW6p_ME = 1.0;
        
		if (doDebug) std::cout << "define branches for the new tree " << std::endl;
		newtree->Branch("EventWeight", &eventweight, "EventWeight/D"); 
		newtree->Branch("Mll", &nMll, "Mll/F");
		newtree->Branch("MT", &nMT, "MT/F");
		newtree->Branch("MT_TrackHWW_Clj", &nMT_TrackHWW_Clj, "MT_TrackHWW_Clj/F");
		newtree->Branch("Mtt", &nMtt, "Mtt/F");
		newtree->Branch("Mtt_TrackHWW_Clj", &nMtt_TrackHWW_Clj, "Mtt_TrackHWW_Clj/F");
		newtree->Branch("Ptll", &nPtll, "Ptll/F");
		newtree->Branch("DPhill", &nDPhill, "DPhill/F");
		newtree->Branch("METRel", &nMETRel, "METRel/F");
		newtree->Branch("MET", &nMET, "MET/F");
		newtree->Branch("MET_phi", &nMET_phi, "MET_phi/F");
		newtree->Branch("METRel_TrackHWW_Clj", &nMETRel_TrackHWW_Clj, "METRel_TrackHWW_Clj/F");
		newtree->Branch("MET_TrackHWW_Clj", &nMET_TrackHWW_Clj, "MET_TrackHWW_Clj/F");
		newtree->Branch("MET_phi_TrackHWW_Clj", &nMET_phi_TrackHWW_Clj, "MET_phi_TrackHWW_Clj/F");
		newtree->Branch("MW0_TrackHWW_Clj", &nMW0_TrackHWW_Clj, "MW0_TrackHWW_Clj/F");
		newtree->Branch("MW1_TrackHWW_Clj", &nMW1_TrackHWW_Clj, "MW1_TrackHWW_Clj/F");
		newtree->Branch("lepPt0", &nlepPt0, "lepPt0/F");
		newtree->Branch("lepPt1", &nlepPt1, "lepPt1/F");
		newtree->Branch("lepEta0", &nlepEta0, "lepEta0/F");
		newtree->Branch("lepEta1", &nlepEta1, "lepEta1/F");
		newtree->Branch("lepPhi0", &nlepPhi0, "lepPhi0/F");
		newtree->Branch("lepPhi1", &nlepPhi1, "lepPhi1/F");
		newtree->Branch("jetPt0", &njetPt0, "jetPt0/F");
		newtree->Branch("jetEta0", &njetEta0, "jetEta0/F");
		newtree->Branch("jetPhi0", &njetPhi0, "jetPhi0/F");
		newtree->Branch("MinDPhi_TrackHWW_Clj", &nMinDPhi_TrackHWW_Clj, "MinDPhi_TrackHWW_Clj/F");
		newtree->Branch("DEtall", &nDEtall, "DEtall/F");

		newtree->Branch("m_el_n", &nel_n, "m_el_n/I");
		newtree->Branch("m_mu_n", &nmu_n, "m_mu_n/I");
		newtree->Branch("m_jet_n", &njet_n,"m_jet_n/I");
		newtree->Branch("nbt", &nBJets, "nBJets/I");
		newtree->Branch("RunNumber", &runnumber, "run/I");
		newtree->Branch("EventNumber", &nEventNumber, "EventNumber/I");
		newtree->Branch("mc_channel_number", &nmc_channel_number, "mc_channel_number/I");
		newtree->Branch("mva_weight", &mva_weight, "mva_weight/D");
		if(runSpin2training) newtree->Branch("spin2_mva_weight", &spin2_mva_weight, "spin2_mva_weight/D");
		if(useXjetsSubtraining) {
			newtree->Branch("xjets_mva_weight", &xjets_mva_weight, "xjets_mva_weight/D");
			if(runSpin2training) newtree->Branch("spin2_xjets_mva_weight", &spin2_xjets_mva_weight, "spin2_xjets_mva_weight/D");
		}

		if( !(TString(fnames[samp_i]).Contains("data") ) ) {
			newtree->Branch("mcevt_pdf_id1", &mcevt_pdf_id1, "mcevt_pdf_id1/I");
			newtree->Branch("mcevt_pdf_id2", &mcevt_pdf_id2, "mcevt_pdf_id2/I");
			newtree->Branch("mcevt_pdf_x1", &mcevt_pdf_x1, "mcevt_pdf_x1/D");
			newtree->Branch("mcevt_pdf_x2", &mcevt_pdf_x2, "mcevt_pdf_x2/D");
			newtree->Branch("mcevt_pdf_scale", &mcevt_pdf_scale, "mcevt_pdf_scale/D");
			newtree->Branch("MTtruth", &nMTtruth, "MTtruth/F");
		}
		newtree->Branch("lepID0", &nlepID0, "lepID0/F");
		newtree->Branch("lepID1", &nlepID1, "lepID1/F");
		if(sys=="Nominal") {
			newtree->Branch("MVAEventWeight", &nMVAEventWeight, "MVAEventWeight/D");
			newtree->Branch("pileupEventWeight", &npileupEventWeight, "pileupEventWeight/D");
			newtree->Branch("pileupEventWeight_090", &npileupEventWeight_090, "pileupEventWeight_090/D");
			newtree->Branch("pileupEventWeight_080", &npileupEventWeight_080, "pileupEventWeight_080/D");
			newtree->Branch("lepTrigSFEventWeight", &nlepTrigSFEventWeight, "lepTrigSFEventWeight/D");
			newtree->Branch("lepTrigSFEventWeightUp", &nlepTrigSFEventWeightUp, "lepTrigSFEventWeightUp/D");
			newtree->Branch("lepTrigSFEventWeightDown", &nlepTrigSFEventWeightDown, "lepTrigSFEventWeightDown/D");
			newtree->Branch("lepSF0Error", &nlepSF0Error, "lepSF0Error/D");
			newtree->Branch("lepSF1Error", &nlepSF1Error, "lepSF1Error/D");
			newtree->Branch("lepSF0ErrorIso", &nlepSF0ErrorIso, "lepSF0ErrorIso/D");
			newtree->Branch("lepSF1ErrorIso", &nlepSF1ErrorIso, "lepSF1ErrorIso/D");
			newtree->Branch("lepSF0EventWeight", &nlepSF0EventWeight, "lepSF0EventWeight/D");
			newtree->Branch("lepSF1EventWeight", &nlepSF1EventWeight, "lepSF1EventWeight/D");
			newtree->Branch("MV120_85_EventWeight", &nMV120_85_EventWeight, "MV120_85_EventWeight/D");
			newtree->Branch("MV120_85_CTJetWeight", &nMV120_85_CTJetWeight, "MV120_85_CTJetWeight/D");
			newtree->Branch("MV120_85_BJetWeight", &nMV120_85_BJetWeight, "MV120_85_BJetWeight/D");
			newtree->Branch("MV120_85_MisTagWeight", &nMV120_85_MisTagWeight, "MV120_85_MisTagWeight/D");
			newtree->Branch("MV120_85_CTJetWeightUp", &nMV120_85_CTJetWeightUp, "MV120_85_CTJetWeightUp/D");
			newtree->Branch("MV120_85_BJetWeightUp", &nMV120_85_BJetWeightUp, "MV120_85_BJetWeightUp/D");
			newtree->Branch("MV120_85_MisTagWeightUp", &nMV120_85_MisTagWeightUp, "MV120_85_MisTagWeightUp/D");
			newtree->Branch("MV120_85_CTJetWeightDown", &nMV120_85_CTJetWeightDown, "MV120_85_CTJetWeightDown/D");
			newtree->Branch("MV120_85_BJetWeightDown", &nMV120_85_BJetWeightDown, "MV120_85_BJetWeightDown/D");
			newtree->Branch("MV120_85_MisTagWeightDown", &nMV120_85_MisTagWeightDown, "MV120_85_MisTagWeightDown/D");
			if(TString(fnames[samp_i]).Contains("Wjets_datadriven") || TString(fnames[samp_i]).Contains("QCD_datadriven") )
				{
					newtree->Branch("SysFakeWeightCorr", &nSysFakeWeight,"SysFakeWeightCorr/D");
					newtree->Branch("SysFakeStatCorr", &nSysFakeStat,"SysFakeStatCorr/D");
					newtree->Branch("SysFakeFlavCorr", &nSysFakeFlav,"SysFakeFlavCorr/D");
					newtree->Branch("SysFakeOtherCorr", &nSysFakeOther,"SysFakeOtherCorr/D");
					newtree->Branch("isFake0", &nisFake0,"isFake0/I");
					newtree->Branch("isFake1", &nisFake1,"isFake1/I");
				}

			if(TString(fnames[samp_i]).Contains("CPMixing_mg5")){
				//-------- weights:
				if(dokHWWScan) {
					newtree->Branch("kHWW0dot1p_ME",&kHWW1p_ME,"kHWW0dot1p_ME/D");
					newtree->Branch("kHWW0dot1m_ME",&kHWW2p_ME,"kHWW0dot1m_ME/D");
					newtree->Branch("kHWW0dot15p_ME",&kHWW3p_ME,"kHWW0dot15p_ME/D");
					newtree->Branch("kHWW0dot15m_ME",&kHWW4p_ME,"kHWW0dot15m_ME/D");
					newtree->Branch("kHWW0dot2p_ME",&kHWW5p_ME,"kHWW0dot2p_ME/D");
					newtree->Branch("kHWW0dot2m_ME",&kHWW6p_ME,"kHWW0dot2m_ME/D");
				} else {
					newtree->Branch("ca08m_ME",&ca00_ME,"ca08m_ME/D");
					newtree->Branch("ca085m_ME",&ca01p_ME,"ca085m_ME/D");
					newtree->Branch("ca09m_ME",&ca02p_ME,"ca09m_ME/D");
					newtree->Branch("ca095_ME",&ca03p_ME,"ca095_ME/D");
					newtree->Branch("ca098_ME",&ca04p_ME,"ca098_ME/D");
					newtree->Branch("ca099_ME",&ca05p_ME,"ca099_ME/D");
				}
				newtree->Branch("source_ME",&source_ME,"source_ME/D");
				newtree->Branch("HPttruth", &HPttruth, "HPttruth/F");
			}

		}
		newtree->Branch("isBlinded", &nisBlinded, "isBlinded/I");
		newtree->Branch("BDPhill",&nBDPhill,"BDPhill/F");
		newtree->Branch("BDPsill",&nBDPsill,"BDPsill/F");
		newtree->Branch("BPLeadLep",&nBPLeadLep,"BPLeadLep/F");
		newtree->Branch("BPSubLeadLep",&nBPSubLeadLep,"BPSubLeadLep/F");
		newtree->Branch("BELeadNeu",&nBELeadNeu,"BELeadNeu/F");
		newtree->Branch("BESubLeadNeu",&nBESubLeadNeu,"BESubLeadNeu/F");
		newtree->Branch("Esum",&nEsum,"Esum/F");
		newtree->Branch("DPt",&nDPt,"DPt/F");
		newtree->Branch("HPt",&nHPt,"HPt/F");
		newtree->Branch("Efun",&Efun,"Efun/F");
    if( sys=="Nominal" && TString(fnames[samp_i]).Contains("ggH125_EF") ) newtree->Branch("higgsPt",&nhiggsPt,"higgsPt/F");  


		if(doDebug) std:: cout << "loop over the entries, " <<t_input->GetEntries()<< std::endl;

		//settings for source
		afloat_t source_ca = 0.3;
		afloat_t source_kSM = 1.0;
		afloat_t source_kHWW = 2.0;
		afloat_t source_kAWW = 2.0;
		afloat_t source_kHdw = 1.0;
		afloat_t source_Lambda = 125.5;
		//settings for target
		afloat_t target_ca = 1.0;
		afloat_t target_kSM = 1.0;
		afloat_t target_kHWW = 0.0;
		afloat_t target_kAWW = 0.0;
		afloat_t target_kHdw = 0.0;
		afloat_t target_Lambda = 125.5;

		if(dokHWWScan){
			if(samp_i==21) {target_kHWW   = 0.1; target_kSM = 1;} //change back kSM to 0 for inf
			if(samp_i==22) {target_kHWW   = -0.1; target_kSM = 1;}
			if(samp_i==23) {target_kHWW   = 0.15;  target_kSM = 1;}
			if(samp_i==24) {target_kHWW   = -0.15;  target_kSM = 1;}
			if(samp_i==25) {target_kHWW   = 0.2;  target_kSM = 1;}
			if(samp_i==26) {target_kHWW   = -0.2;  target_kSM = 1;}

		} else {
			if(samp_i==21) {target_ca     = -0.8; target_kAWW   = 1.0;}
			if(samp_i==22) {target_ca     = -0.85; target_kAWW   = 1.0;}
			if(samp_i==23) {target_ca     = -0.9; target_kAWW   = 1.0;}
			if(samp_i==24) {target_ca     = 0.95; target_kAWW   = 1.0;}
			if(samp_i==25) {target_ca     = 0.98; target_kAWW   = 1.0;}
			if(samp_i==26) {target_ca     = 0.99; target_kAWW   = 1.0;}
		}


		for (Long64_t ev_i = 0; ev_i<t_input->GetEntries(); ev_i++)
			{
				t_input->GetEntry(ev_i);
				if(doDebug) std:: cout << "entry = " << ev_i << std::endl;
        //if(ev_i % 1000 == 0) std::cout <<"Number of events processed " << ev_i << std::endl;

				nMVAEventWeight = MVAEventWeight;
				nMT = MT;
				nMT_TrackHWW_Clj = MT_TrackHWW_Clj;
				nMTtruth=MTtruth;
				nMll = Mll;
				nMtt = Mtt;
				float negMtt=1.;
				Mtt_TrackHWW_Clj=max(negMtt,Mtt_TrackHWW_Clj);
				nMtt_TrackHWW_Clj = Mtt_TrackHWW_Clj;
				nPtll = Ptll;
				nDPhill = DPhill;
				nMETRel = METRel;
				nMET = MET;
				nMET_phi = MET_phi;
				nMETRel_TrackHWW_Clj = METRel_TrackHWW_Clj;
				nMET_TrackHWW_Clj = MET_TrackHWW_Clj;
				nMET_phi_TrackHWW_Clj = MET_phi_TrackHWW_Clj;
				nlepID0 = lepID0;
				nlepID1 = lepID1;
				nlepPt0 = lepPt0;
				nlepPt1 = lepPt1;
				nlepEta0 = lepEta0;
				nlepEta1 = lepEta1;
				nlepPhi0 = lepPhi0;
				nlepPhi1 = lepPhi1;
				njetPt0 = jetPt0;
				njetEta0 = jetEta0;
				njetPhi0 = jetPhi0;
				nBDPhill = BDPhill;
				nBDPsill = BDPsill;
				nBPLeadLep = BPLeadLep;
				nBPSubLeadLep = BPSubLeadLep;
				nBELeadNeu = BELeadNeu;
				nBESubLeadNeu = BESubLeadNeu;
				DPt=fabs(DPt);
				nDPt = DPt;
				nEsum = Esum;
				nHPt=HPt;
				nMinDPhi_TrackHWW_Clj = MinDPhi_TrackHWW_Clj;
				nDEtall=DEtall;
				if(lepPt0>lepPt1) {
					leadLepPt = lepPt0;
					subleadLepPt = lepPt1;
				} else {
					leadLepPt = lepPt1;
					subleadLepPt = lepPt0;
				}

				if(lepPt0>lepPt1) Efun = lepPt0 - 0.5*lepPt1 + 0.5*MET_TrackHWW_Clj;
				else Efun = lepPt1 - 0.5*lepPt0 + 0.5*MET_TrackHWW_Clj;

				nel_n = m_el_n;
				nmu_n = m_mu_n;
				njet_n = m_jet_n;

				if(!(TString(fnames[samp_i]).Contains("Wjets_datadriven") || TString(fnames[samp_i]).Contains("QCD_datadriven"))) {
					runnumber = RunNumber;
					nEventNumber = EventNumber;
				} else {
					runnumber = wjets_RunNumber;
					nEventNumber = wjets_EventNumber;
					nSysFakeWeight = SysFakeWeight;
					nSysFakeStat = SysFakeStat;
					nSysFakeFlav = SysFakeFlav;
					nSysFakeOther = SysFakeOther;
					nisFake0 = isFake0;
					nisFake1 = isFake1;

				}
				nmc_channel_number = mc_channel_number;
				nBJets = nJets_Pt20_MV1_85;
				npileupEventWeight = pileupEventWeight;
				npileupEventWeight_090 = pileupEventWeight_090;
				npileupEventWeight_080 = pileupEventWeight_080;
				nlepTrigSFEventWeight = lepTrigSFEventWeight;
				nlepTrigSFEventWeightUp = lepTrigSFEventWeightUp;
				nlepTrigSFEventWeightDown = lepTrigSFEventWeightDown;
				nlepSF0Error = lepSF0Error;
				nlepSF1Error = lepSF1Error;
				nlepSF0ErrorIso = lepSF0ErrorIso;
				nlepSF1ErrorIso = lepSF1ErrorIso;
				nlepSF0EventWeight = lepSF0EventWeight;
				nlepSF1EventWeight = lepSF1EventWeight;
				nMV120_85_EventWeight = MV120_85_EventWeight;
				nMV120_85_CTJetWeight = MV120_85_CTJetWeight; 
				nMV120_85_BJetWeight = MV120_85_BJetWeight; 
				nMV120_85_MisTagWeight = MV120_85_MisTagWeight;
				nMV120_85_CTJetWeightUp = MV120_85_CTJetWeightUp; 
				nMV120_85_CTJetWeightDown = MV120_85_CTJetWeightDown; 
				nMV120_85_MisTagWeightUp = MV120_85_MisTagWeightUp;
				nMV120_85_MisTagWeightDown = MV120_85_MisTagWeightDown;
				nMV120_85_BJetWeightUp = MV120_85_BJetWeightUp; 
				nMV120_85_BJetWeightDown = MV120_85_BJetWeightDown;
				nisBlinded = isBlinded;

				nMW0_TrackHWW_Clj = sqrt(2*lepPt0*MET_TrackHWW_Clj*(1-cos(TMath::ATan2(MET_y_TrackHWW_Clj,MET_x_TrackHWW_Clj)-lepPhi0)));
				nMW1_TrackHWW_Clj = sqrt(2*lepPt1*MET_TrackHWW_Clj*(1-cos(TMath::ATan2(MET_y_TrackHWW_Clj,MET_x_TrackHWW_Clj)-lepPhi1)));

				if(lepPt0>lepPt1)	leadLepMT = nMW0_TrackHWW_Clj;
				else leadLepMT = nMW1_TrackHWW_Clj;

				if( !(TString(fnames[samp_i]).Contains("data")) ) {
					mcevt_pdf_id1 = m_mcevt_pdf_id1->at(0);
					mcevt_pdf_id2 = m_mcevt_pdf_id1->at(0);
					mcevt_pdf_x1 = m_mcevt_pdf_x1->at(0);
					mcevt_pdf_x2 = m_mcevt_pdf_x2->at(0);
					mcevt_pdf_scale = m_mcevt_pdf_scale->at(0);
				} else {
					mcevt_pdf_id1 = -1;
					mcevt_pdf_id2 = -1;
					mcevt_pdf_x1 = -1;
					mcevt_pdf_x2 = -1;
					mcevt_pdf_scale = -1;
				}

				if(overlapWZ) intoverlapWZ = 1;
				else if(!overlapWZ) intoverlapWZ = 0;

				nhiggsPt=higgsPt;
        bool good_event = true;
				TString channel = "";
				TString jetbin = "";

				//mva ntuples are after the metrel cut
				if(isLowPtCand==1) good_event=false;
				if(!(TString(fnames[samp_i]).Contains("Wjets_datadriven") || TString(fnames[samp_i]).Contains("QCD_datadriven"))) {
					if(!passMueORsimple) good_event=false;
				} else {
					if(!intpassMueORsimple) good_event=false;
				}

				//select lepton channel
				if(m_el_n==2) channel = "ee";
				else if(m_mu_n==2) channel = "mumu";
				else if(m_mu_n==1 && m_el_n==1) channel = "emu";
				if( channel=="" ) good_event = false;
				if( channel=="ee" ) good_event = false;
				if( channel=="mumu" ) good_event = false;

				//metrel cut
				if(doMVAPropertiesCuts) {
					if( MET_TrackHWW_Clj<= 20000. ) good_event = false;
				} else {
					if( m_jet_n < 2 && MET_TrackHWW_Clj<= 25000 ) good_event = false;
					if( (m_el_n==2 || m_mu_n==2)  && MET_TrackHWW_Clj<= 40000) good_event = false;
				}
		
				//additional cuts for Wjets data-driven -- no mva ntuple cuts
				if (TString(fnames[samp_i]).Contains("Wjets_datadriven") || TString(fnames[samp_i]).Contains("QCD_datadriven")) {
					if(lepID0*lepID1 > 0)  good_event = false; // default: OS
					if(m_el_n==1 && m_mu_n==1 && Mll < 10000)  good_event = false;
					if( (m_el_n==2 || m_mu_n==2) && Mll < 12000)  good_event = false;
					if( (m_el_n==2 || m_mu_n==2)  && fabs(Mll - 91187.6) < 15000.)  good_event = false;
				}

				//apply jet cut
				if(m_jet_n==0) jetbin = "0jet";
				else if(m_jet_n==1) jetbin = "1jet";
				else if(m_jet_n > 1) jetbin = "2jet";
				if(jetbin == "") good_event = false;
				if(jetbin == "2jet") good_event = false;

				if(doTopDataDrivenNtuples) {
					if(nBJets==0) good_event = false;
				}

				if(good_event){
					if(doDebug) std:: cout << "good event " << std::endl;
					Double_t reweightFactor = 1.;
					// this is mix ca00
					if( (TString(fnames[samp_i]).Contains("CPMixing_mg5"))) {
						afloat_t event_data[4*6*1];
						TLorentzVector lm, lp, vbar, v;
						TLorentzVector ga, gb;
						ga.SetPtEtaPhiM(mc_pt->at(0)/1000., mc_eta->at(0), mc_phi->at(0), mc_m->at(0)/1000.);
						gb.SetPtEtaPhiM(mc_pt->at(1)/1000., mc_eta->at(1), mc_phi->at(1), mc_m->at(1)/1000.);
            int all_part = 0;
						for(int i=0; i<mc_n; i++){                    
							if((mc_pdgId->at(i) == 11 || mc_pdgId->at(i)== 13 || mc_pdgId->at(i) == 15) && mc_status->at(i)==3){
								lm.SetPtEtaPhiM(mc_pt->at(i)/1000. ,mc_eta->at(i), mc_phi->at(i), mc_m->at(i)/1000.);
								all_part++;
							}
							else if((mc_pdgId->at(i) == -11 || mc_pdgId->at(i)== -13 || mc_pdgId->at(i) == -15) && mc_status->at(i)==3){
								lp.SetPtEtaPhiM(mc_pt->at(i)/1000. ,mc_eta->at(i), mc_phi->at(i), mc_m->at(i)/1000.);
								all_part++;
							}
							else if((mc_pdgId->at(i) == 12 || mc_pdgId->at(i)== 14 || mc_pdgId->at(i)== 16) && mc_status->at(i)==3){
								v.SetPtEtaPhiM(mc_pt->at(i)/1000. ,mc_eta->at(i), mc_phi->at(i), mc_m->at(i)/1000.);
								all_part++;
							}
							else if((mc_pdgId->at(i) == -12 || mc_pdgId->at(i)== -14 || mc_pdgId->at(i)== -16) && mc_status->at(i)==3){
								vbar.SetPtEtaPhiM(mc_pt->at(i)/1000. ,mc_eta->at(i), mc_phi->at(i), mc_m->at(i)/1000.);
								all_part++;
							}
							if(mc_pdgId->at(i) == 5000000 && mc_status->at(i)==3) {
								HPttruth = mc_pt->at(i);
								//std::cout << "status = "<< mc_status->at(i)<< " HPt = " << mc_pt->at(i) << std::endl;
							}
						}// loop over mc_n
             
						if(all_part==4) {
							TLorentzVector particles[6]= {ga, gb, lp, v, lm, vbar};               
							for(int ip=0; ip<6; ++ip){
								event_data[/*jentry*4*6 +*/ ip*4+0] = particles[ip].E();
								event_data[/*jentry*4*6 +*/ ip*4+1] = particles[ip].Px();
								event_data[/*jentry*4*6 +*/ ip*4+2] = particles[ip].Py();
								event_data[/*jentry*4*6 +*/ ip*4+3] = particles[ip].Pz();
							} 

							afloat_t source[1];
							afloat_t result[1];
							
							//source ME
							quadme(source, event_data, 1,0,0,0,0,0,0,0,0,0,0, source_ca, source_kSM, source_kHWW, source_kAWW, source_kHdw, source_Lambda, -1);
							source_ME=source[0];                
							//target ME
							quadme(result, event_data, 1,0,0,0,0,0,0,0,0,0,0, target_ca, target_kSM, target_kHWW, target_kAWW, target_kHdw, target_Lambda, -1);
							if(source[0] != 0.) reweightFactor = result[0]/source[0];
							else reweightFactor = 0.;
							if(dokHWWScan) {
								if(samp_i==21) kHWW1p_ME=result[0];
								else if(samp_i==22) kHWW2p_ME = result[0];
								else if(samp_i==23) kHWW3p_ME = result[0];
								else if(samp_i==24)	kHWW4p_ME = result[0];
								else if(samp_i==25) kHWW5p_ME = result[0];
								else if(samp_i==26) kHWW6p_ME = result[0];
							} else {
								if(samp_i==21) ca00_ME=result[0];
								else if(samp_i==22) ca01p_ME = result[0];
								else if(samp_i==23) ca02p_ME = result[0];
								else if(samp_i==24)	ca03p_ME = result[0];
								else if(samp_i==25) ca04p_ME = result[0];
								else if(samp_i==26) ca05p_ME = result[0];
							}
						} else reweightFactor = 0;
					}//if samp = CPmixing

					if(TString(fnames[samp_i]).Contains("data") && !(TString(fnames[samp_i]).Contains("Wjets_datadriven") ||TString(fnames[samp_i]).Contains("QCD_datadriven")) ) eventweight = 1;
					else if (TString(fnames[samp_i]).Contains("Wjets_datadriven") || TString(fnames[samp_i]).Contains("QCD_datadriven")) {
						eventweight = MVAEventWeight*(1-intoverlapWZ)*WgStarEventWeight;
						if(doTopDataDrivenNtuples) eventweight = -1*eventweight;
					}

					else if(TString(fnames[samp_i]).Contains("CPMixing_mg5")){
						eventweight = invpb*MVAEventWeight*pileupEventWeight_090*(1-intoverlapWZ)*WgStarEventWeight*reweightFactor; //
						//std::cout << "reweighting = " << reweightFactor << std::endl;
						double hpt_reweight=1.;
						if(HPttruth/1000.<5.) hpt_reweight = 0.773397;
						if(HPttruth/1000.>5. && HPttruth/1000.<10.) hpt_reweight = 0.979936;
						if(HPttruth/1000.>10. && HPttruth/1000.<15.) hpt_reweight = 1.04423;
						if(HPttruth/1000.>15. && HPttruth/1000.<20.) hpt_reweight = 1.06319;
						if(HPttruth/1000.>20. && HPttruth/1000.<25.) hpt_reweight = 1.04401;
						if(HPttruth/1000.>25. && HPttruth/1000.<30.) hpt_reweight = 1.03798;
						if(HPttruth/1000.>30. && HPttruth/1000.<35.) hpt_reweight = 1.01462;
						if(HPttruth/1000.>35. && HPttruth/1000.<40.) hpt_reweight = 0.968747;
						if(HPttruth/1000.>40. && HPttruth/1000.<45.) hpt_reweight = 0.959023;
						if(HPttruth/1000.>45. && HPttruth/1000.<50.) hpt_reweight = 0.918;
						if(HPttruth/1000.>50. && HPttruth/1000.<55.) hpt_reweight = 0.99991;
						if(HPttruth/1000.>55. && HPttruth/1000.<60.) hpt_reweight = 0.924358;
						if(HPttruth/1000.>60. && HPttruth/1000.<65.) hpt_reweight = 0.973597;
						if(HPttruth/1000.>65. && HPttruth/1000.<70.) hpt_reweight = 0.851553;
						if(HPttruth/1000.>70. && HPttruth/1000.<75.) hpt_reweight = 1.00951;
						if(HPttruth/1000.>75. && HPttruth/1000.<80.) hpt_reweight = 1.17162;
						if(HPttruth/1000.>80. && HPttruth/1000.<85.) hpt_reweight = 1.0012;
						if(HPttruth/1000.>85. && HPttruth/1000.<90.) hpt_reweight = 1.08455;
						if(HPttruth/1000.>90. && HPttruth/1000.<95.) hpt_reweight = 0.887916;
						eventweight *=  hpt_reweight;
					}
					else{
						eventweight = invpb*MVAEventWeight*pileupEventWeight_090*(1-intoverlapWZ)*WgStarEventWeight;
						if( (TString(fnames[samp_i]).Contains("DrellYan") || TString(fnames[samp_i]).Contains("Zjets")) && !(TString(fnames[samp_i]).Contains("ZjetsEW")) ) eventweight=eventweight*PDFWeight[1];
            
						if(doTopDataDrivenNtuples) eventweight = -1*eventweight;

						if(m_jet_n == 1) eventweight = eventweight*MV120_85_EventWeight;
 
						if(doPtllReweighting && m_jet_n==0 && (TString(fnames[samp_i]).Contains("Zjets") || TString(fnames[samp_i]).Contains("DrellYan") || TString(fnames[samp_i]).Contains("Zgamma"))) {
							double ptllweight = 1.;

							ptllweight = (Ptll_truth/1000. < 2.)*1.14303 
								+ (Ptll_truth/1000. >2. && Ptll_truth/1000. <  4.)*1.08611
								+ (Ptll_truth/1000. >4. && Ptll_truth/1000. <  6.)*1.04975
								+ (Ptll_truth/1000. >6. && Ptll_truth/1000. <  8.)*1.01799
								+ (Ptll_truth/1000. >8. && Ptll_truth/1000. <  10.)*1.00356
								+ (Ptll_truth/1000. >10. && Ptll_truth/1000. <  12.)*0.9887
								+ (Ptll_truth/1000. >12. && Ptll_truth/1000. <  14.)*0.980315
								+ (Ptll_truth/1000. >14. && Ptll_truth/1000. <  16.)*0.977064
								+ (Ptll_truth/1000. >16. && Ptll_truth/1000. <  18.)*0.978649
								+ (Ptll_truth/1000. >18. && Ptll_truth/1000. <  20.)*0.971497
								+ (Ptll_truth/1000. >20. && Ptll_truth/1000. <  25.)*0.992051
								+ (Ptll_truth/1000. >25. && Ptll_truth/1000. <  30.)*1.01695
								+ (Ptll_truth/1000. >30. && Ptll_truth/1000. <  35.)*1.08943
								+ (Ptll_truth/1000. >35. && Ptll_truth/1000. <  40.)*1.13008
								+ (Ptll_truth/1000. >40. && Ptll_truth/1000. <  50.)*1.2204
								+ (Ptll_truth/1000. >50. && Ptll_truth/1000. <  70.)*1.31241
								+ (Ptll_truth/1000. >70. )*1.55979;

							//cout << "eventweight = " << eventweight << ", ptllweight= " << ptllweight << ", prod = " << ptllweight*eventweight<< endl;

							eventweight = ptllweight*eventweight;
						}
            eventweight = isZll*notZgamma*eventweight; 
						// if(  mc_channel_number==167113 || mc_channel_number==167118 || mc_channel_number==167303 || mc_channel_number==167308
						// 		 || mc_channel_number==167114 || mc_channel_number==167119 || mc_channel_number==167304 || mc_channel_number==167309 )
						// 	{
						// 		if( m_mcevt_pdf_id1->at(0) == 21 && m_mcevt_pdf_id1->at(0) == 21 ) eventweight=0*eventweight;
						// 		else eventweight=4*eventweight*pt_weight;
                    
						// 		//std::cout << "higgs pt = " << pt_weight << std::endl;
						// 	}
						if(TString(fnames[samp_i]).Contains("ggH125_EF")) eventweight = eventweight*higgsPtEventWeight; 
					} // calculate eventweight
					if(doDebug) std::cout << "event weight = " << eventweight << std::endl;
					if(doDebug) std::cout << "analysis = " << analysis << std::endl;

					if(jetbin=="0jet") {
						//standard training
						if(doDebug) std::cout << "0 jet, standard training " << analysis << std::endl;
						if(runEvenOdd){
							if(nEventNumber%2) mva_weight = reader_emu_0jet->EvaluateMVA(analysis);
							else mva_weight = reader_evenodd_emu_0jet->EvaluateMVA(analysis);
						} else  mva_weight = reader_emu_0jet->EvaluateMVA(analysis);
						//Xjets sub-training
						if(useXjetsSubtraining) {
							if(doDebug) std::cout << "0 jet, wjets training " << "WjetsBDT_8var_cfg_242" << std::endl;
							if(runEvenOdd){
								if(nEventNumber%2) xjets_mva_weight = reader_subtrain_emu_0jet->EvaluateMVA("WjetsBDT_8var_cfg_242_trainOnEven");
								else xjets_mva_weight = reader_evenodd_subtrain_emu_0jet->EvaluateMVA("WjetsBDT_8var_cfg_242_trainOnOdd");
							} else  xjets_mva_weight = reader_subtrain_emu_0jet->EvaluateMVA("WjetsBDT_8var_cfg_242");
						}
						//Spin2 training
						if(runSpin2training) {
							if(doDebug) std::cout << "0 jet, spin-2 training " << analysis << std::endl;
							if(runEvenOdd){
								if(nEventNumber%2) spin2_mva_weight = reader_train2_emu_0jet->EvaluateMVA(analysis);
								else spin2_mva_weight = reader_evenodd_train2_emu_0jet->EvaluateMVA(analysis);
							} else spin2_mva_weight = reader_train2_emu_0jet->EvaluateMVA(analysis);
						}
					}
					if(jetbin=="1jet") {
						//standard training
						if(runEvenOdd){
							if(doDebug) std::cout << "1 jet, standard training " << analysis << std::endl;
							if(nEventNumber%2) mva_weight = reader_emu_1jet->EvaluateMVA(analysis);
							else mva_weight = reader_evenodd_emu_1jet->EvaluateMVA(analysis);
						} else  mva_weight = reader_emu_1jet->EvaluateMVA(analysis);
						//Xjets sub-training
						if(useXjetsSubtraining) {
							if(doDebug) std::cout << "1 jet, wjets training " << "WjetsBDT_8var_cfg_242" << std::endl;
							if(runEvenOdd){
								if(nEventNumber%2) xjets_mva_weight = reader_subtrain_emu_1jet->EvaluateMVA("WjetsBDT_8var_cfg_242_trainOnEven");
								else xjets_mva_weight = reader_evenodd_subtrain_emu_1jet->EvaluateMVA("WjetsBDT_8var_cfg_242_trainOnOdd");
							} else  xjets_mva_weight = reader_subtrain_emu_1jet->EvaluateMVA("WjetsBDT_8var_cfg_242");
						}
						//Spin2 training
						if(runSpin2training) {
							if(doDebug) std::cout << "1 jet, spin-2 training " << analysis << std::endl;
							if(runEvenOdd){
								if(nEventNumber%2) spin2_mva_weight = reader_train2_emu_1jet->EvaluateMVA(analysis);
								else spin2_mva_weight = reader_evenodd_train2_emu_1jet->EvaluateMVA(analysis);
							} else spin2_mva_weight = reader_train2_emu_1jet->EvaluateMVA(analysis);
						}
					}

					newtree->Fill();

				}//if good event
			}//end loop over events
		target->Write();
		target->Close();
		delete target; 
		
	}//end loop over samples
if(doNtuplesForLimit && runoverallsamples){
		gSystem->Exec( "hadd -f "+output_dir+"/wzzz.root "+output_dir+"/zz.root "+output_dir+"/zzgg.root "+output_dir+"/zz2j.root "+output_dir+"/wz.root "+output_dir+"/wz2j.root " );
		gSystem->Exec( "hadd -f "+output_dir+"/ww.root "+output_dir+"/wwgg.root "+output_dir+"/wwqq.root "+output_dir+"/ww2j.root" );
		gSystem->Exec( "hadd -f "+output_dir+"/zleplep.root "+output_dir+"/zhighjets.root "+output_dir+"/dy.root "+output_dir+"/zg.root "+output_dir+"/zew.root" );
		gSystem->Exec( "hadd -f "+output_dir+"/ztautau.root "+output_dir+"/zhighjetstautau.root "+output_dir+"/dytautau.root "+output_dir+"/zgtautau.root "+output_dir+"/zewtautau.root" );
		gSystem->Exec( "rm "+output_dir+"/zz.root "+output_dir+"/zzgg.root "+output_dir+"/zz2j.root "+output_dir+"/wz.root "+output_dir+"/wz2j.root" );
		gSystem->Exec( "rm "+output_dir+"/wwgg.root "+output_dir+"/wwqq.root "+output_dir+"/ww2j.root" );
		gSystem->Exec( "rm "+output_dir+"/zhighjets.root "+output_dir+"/zhighjetstautau.root "+output_dir+"/dy.root "+output_dir+"/dytautau.root "+output_dir+"/zg.root "+output_dir+"/zgtautau.root "+output_dir+"/zew.root "+output_dir+"/zewtautau.root" );
		//if(sys.Contains("Nominal")) gSystem->Exec( "hadd -f "+output_dir+"/data.root "+output_dir+"/dataUnblind.root "+output_dir+"/dataBlind.root");
	}

	delete reader_emu_0jet;
	delete reader_emu_1jet;
	if(runEvenOdd) {
		delete reader_evenodd_emu_0jet;
		delete reader_evenodd_emu_1jet;
	}
	if(runSpin2training) {
		delete reader_train2_emu_0jet;
		delete reader_train2_emu_1jet;
		if(runEvenOdd) {
			delete reader_evenodd_train2_emu_0jet;
			delete reader_evenodd_train2_emu_1jet;
		}
	}
	if(useXjetsSubtraining) {
		delete reader_subtrain_emu_0jet;
		delete reader_subtrain_emu_1jet;
		if(runEvenOdd) {
			delete reader_evenodd_subtrain_emu_0jet;
			delete reader_evenodd_subtrain_emu_1jet;
		}
	}

	std::cout <<  "done!"<<std::endl;

}
