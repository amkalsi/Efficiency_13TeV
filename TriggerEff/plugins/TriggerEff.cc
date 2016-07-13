// -*- C++ -*-
//
// Package:    NtupleMaker/TriggerEff
// Class:      TriggerEff
// 
/**\class TriggerEff TriggerEff.cc NtupleMaker/TriggerEff/plugins/TriggerEff.cc

Description: [one line class summary]

Implementation:
[Notes on implementation]
*/
//
// Original Author:  kaur amandeepkalsi
//         Created:  Wed, 19 Aug 2015 12:59:54 GMT
//
//


// system include files
#include <memory>
#include "TLorentzVector.h"
#include "THn.h"
#include <THnSparse.h>
#include "TH2.h"
// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/PatCandidates/interface/TriggerObjectStandAlone.h"
#include "DataFormats/PatCandidates/interface/PackedTriggerPrescales.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include <DataFormats/PatCandidates/interface/GenericParticle.h>
#include <DataFormats/PatCandidates/interface/PackedGenParticle.h>
#include <DataFormats/HepMCCandidate/interface/GenStatusFlags.h>
#include <DataFormats/Candidate/interface/Candidate.h>
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Tau.h"
#include "DataFormats/PatCandidates/interface/Photon.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/Math/interface/normalizedPhi.h" 
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "DataFormats/PatCandidates/interface/TriggerPath.h"
#include "DataFormats/PatCandidates/interface/TriggerEvent.h"
#include "DataFormats/PatCandidates/interface/TriggerAlgorithm.h"
#include "DataFormats/HLTReco/interface/TriggerEvent.h"
#include "DataFormats/HLTReco/interface/TriggerFilterObjectWithRefs.h"

#include "FWCore/Common/interface/TriggerResultsByName.h"
#include "DataFormats/HLTReco/interface/TriggerRefsCollections.h"
#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"
#include "TTree.h"
#include "TBranch.h"

using namespace std;
using namespace pat;
using namespace reco;
//
// class declaration
//

class TriggerEff : public edm::EDAnalyzer {
	public:
		explicit TriggerEff(const edm::ParameterSet&);
		~TriggerEff();

		static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


	private:
		virtual void beginJob() override;
		virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
		virtual void endJob() override;

		virtual void beginRun(edm::Run const&, edm::EventSetup const&) override;
		virtual void endRun(edm::Run const&, edm::EventSetup const&) override;
		//virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
		//virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;

		// ----------member data ---------------------------
		edm::EDGetTokenT<edm::ValueMap<bool> > electronVetoIdMapToken_;
		edm::EDGetTokenT<edm::ValueMap<bool> > electronLooseIdMapToken_;
		edm::EDGetTokenT<edm::ValueMap<bool> > electronMediumIdMapToken_;
		edm::EDGetTokenT<edm::ValueMap<bool> > electronTightIdMapToken_;
		edm::EDGetTokenT<edm::ValueMap<bool> > eleHEEPIdMapToken_;

		string HLTPath1, HLTFilter1a, HLTFilter1b;
		string HLTPath2, HLTFilter2a,HLTFilter2b,HLTFilter2c,HLTFilter2d;                                                        
		string HLTPath3 , HLTFilter3a, HLTFilter3b;
		string HLTPath4 , HLTFilter4a, HLTFilter4b;  	
		edm::InputTag triggerResultsTag_;
		TLorentzVector genmuinfo, gentauinfo;
		const reco::Candidate* GetLastCopy (const reco::Candidate* part);
		int GetTauDecay (const reco::Candidate* part);
		int GetTauDecay (const reco::GenParticle& part);
		unsigned int indexmuTau , indexTauh;
		const reco::Candidate* GetFirstCopy (const reco::Candidate* part);
		bool IsFirstCopy (const reco::Candidate* part, const bool checkAbsPdg); 
		reco::GenParticle GetTauHad (const reco::Candidate* part);
		bool isGoodVertex(const reco::Vertex& vtx);	
		float mTCalculation(const pat::Muon& muobject, const pat::MET& metobject);
		float PZetaVis(const pat::Muon& muobject, const pat::Tau& tauobject);
		float PZeta(const pat::Muon& muobject, const pat::Tau& tauobject, const pat::MET& metobject);
		bool TauSelection( const pat::Tau &myTau, const reco::Vertex &vtx );
		bool MuSelection( const pat::Muon &myMuon, const reco::Vertex &vtx );
		bool ExtraMuon(const pat::Muon& muobject, edm::Handle<pat::MuonCollection> muons, const reco::Vertex &vtx);
		bool BJetsinEvent( edm::Handle<pat::JetCollection> myjets,  const pat::Muon &myMuon, const pat::Tau &myTau);
		bool matchedTrigger1, matchedTrigger2 ;
		// matched trigger objects
		pair<bool, TLorentzVector> MuonMatchingToTrigObjects(const pat::Muon& myMuon, edm::Handle<pat::TriggerObjectStandAloneCollection> triggerObjects, const edm::TriggerNames &names,string filterName, string pathname);                 
		pair<vector<bool>,vector<pat::TriggerObjectStandAlone>> preTauMatchingToTrigObjects(const pat::Tau& myTau, edm::Handle<pat::TriggerObjectStandAloneCollection> triggerObjects, const edm::TriggerNames &names,string filterName, string pathname);      
		pair<bool, TLorentzVector> TauMatchingToTrigObjects(const pat::Tau& myTau, edm::Handle<pat::TriggerObjectStandAloneCollection> triggerObjects, const edm::TriggerNames &names,string filterName, string pathname, bool numerator); //numerator boolean is added to not check the trigger path  
		float deltaPhi( float a, float b);
		float dR(float l1eta, float l1phi, float l2eta, float l2phi );       
		edm::Service<TFileService> fs;
		TH1F *hFillPtDen;
		TH2F *h2D_den, *h2D_num; 
		THnSparse *hSumWeight;
		TH1F *hFillPtNum;
		TH1F *hFillEtaDen;
		TH1F *hFillEtaNum;
		TH1F *hNumberOfMuons, *hNumberOftaus;
		TH1F *hFillEventTotal, *hFillEventMuon, *hFillEventTau;
		TH1F *hmuonid, *hmuoniso, *htauid, *hchargereq, *hmtcut, *hpzeta, *hmuonmatch, *htaumatch1, *htaumatch2;
		double TauPtCut,TauEtaCut,TauIsoCutMax,TauIsoCutMin, MuonPtCut, MuonEtaCut, IsoMuonMax, MotherpdgID;
		string TauDMF,TauEleVeto,TauMuVeto,TauIsoString;
		bool isOSCharge;
		bool isMC, isZtau, isZprime, GenReq ;double weightevt, sumweight, DYOthersBG;
		TTree *myTree;
		double totalweight;
		bool iselectron;
		HLTConfigProvider hltConfig_;
		edm::InputTag triggerPrescalesL1max_;
		edm::InputTag triggerPrescalesL1min_;
		float myMT;
		TH1F *tauptdis,*tauetadis,*muptdis,*muetadis,*mutaumassdis, *mtdis;
};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
TriggerEff::TriggerEff(const edm::ParameterSet& iConfig):
	electronVetoIdMapToken_(consumes<edm::ValueMap<bool> >(iConfig.getParameter<edm::InputTag>("electronVetoIdMap"))),
	electronLooseIdMapToken_(consumes<edm::ValueMap<bool> >(iConfig.getParameter<edm::InputTag>("electronLooseIdMap"))),
	electronMediumIdMapToken_(consumes<edm::ValueMap<bool> >(iConfig.getParameter<edm::InputTag>("electronMediumIdMap"))),
	electronTightIdMapToken_(consumes<edm::ValueMap<bool> >(iConfig.getParameter<edm::InputTag>("electronTightIdMap"))),
	eleHEEPIdMapToken_(consumes<edm::ValueMap<bool> >(iConfig.getParameter<edm::InputTag>("eleHEEPIdMap")))
{
	//now do what ever initialization is needed
	triggerResultsTag_  = iConfig.getParameter<edm::InputTag>("triggerResults");
	HLTPath1 = iConfig.getParameter<string>("HLTPath1");
	HLTFilter1a = iConfig.getParameter<string>("HLTFilter1a");
	HLTFilter1b = iConfig.getParameter<string>("HLTFilter1b");
	HLTPath2 = iConfig.getParameter<string>("HLTPath2");
	HLTFilter2a = iConfig.getParameter<string>("HLTFilter2a");
	HLTFilter2b = iConfig.getParameter<string>("HLTFilter2b");
	/////////////////////////
	TauPtCut = iConfig.getParameter<double>("TauPtCut");
	TauEtaCut  = iConfig.getParameter<double>("TauEtaCut");
	TauDMF  = iConfig.getParameter<string>("TauDMF");
	TauEleVeto = iConfig.getParameter<string>("TauEleVeto");
	TauMuVeto = iConfig.getParameter<string>("TauMuVeto");
	TauIsoString = iConfig.getParameter<string>("TauIsoString");
	TauIsoCutMax = iConfig.getParameter<double>("TauIsoCutMax");
	TauIsoCutMin= iConfig.getParameter<double>("TauIsoCutMin");
	isOSCharge = iConfig.getParameter<bool>("isOSCharge");
	MuonPtCut = iConfig.getParameter<double>("MuonPtCut");
	MuonEtaCut = iConfig.getParameter<double>("MuonEtaCut");
	IsoMuonMax  = iConfig.getParameter<double>("IsoMuonMax");

	MotherpdgID = iConfig.getParameter<double>("MotherpdgID");
	isMC = iConfig.getParameter<bool>("isMC");
	isZtau = iConfig.getParameter<bool>("isZtau");
	isZprime = iConfig.getParameter<bool>("isZprime");
	GenReq = iConfig.getParameter<bool>("GenReq");
	DYOthersBG = iConfig.getParameter<bool>("DYOthersBG");
	triggerPrescalesL1max_ = iConfig.getParameter<edm::InputTag>("l1max");
	triggerPrescalesL1min_ = iConfig.getParameter<edm::InputTag>("l1min");
}


TriggerEff::~TriggerEff()
{

	// do anything here that needs to be done at desctruction time
	// (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
	void
TriggerEff::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
	using namespace edm;

	edm::Handle<edm::TriggerResults> triggerBits;
	edm::Handle<pat::TriggerObjectStandAloneCollection> triggerObjects;

	edm::Handle<pat::PackedTriggerPrescales> triggerPrescales;
	edm::Handle<pat::MuonCollection> muons;
	edm::Handle<pat::TauCollection> taus;
	edm::Handle<pat::METCollection> met;
	edm::Handle<pat::JetCollection> jets;
	edm::Handle<edm::View<pat::Electron> > electron_pat;
	iEvent.getByLabel("slimmedElectrons",electron_pat);

	edm::Handle<pat::PackedTriggerPrescales> triggerPrescalesL1max;
	edm::Handle<pat::PackedTriggerPrescales> triggerPrescalesL1min;
	iEvent.getByLabel(triggerPrescalesL1max_, triggerPrescalesL1max);
	iEvent.getByLabel(triggerPrescalesL1min_, triggerPrescalesL1min);

	iEvent.getByLabel("slimmedMuons",muons);
	iEvent.getByLabel("slimmedTaus",taus);
	iEvent.getByLabel("slimmedMETs",met);
	iEvent.getByLabel("slimmedJets", jets);

	iEvent.getByLabel(triggerResultsTag_, triggerBits);
	iEvent.getByLabel("selectedPatTrigger", triggerObjects);
	iEvent.getByLabel("patTrigger", triggerPrescales);
	edm::Handle<reco::VertexCollection> vtx_h;
	iEvent.getByLabel("offlineSlimmedPrimaryVertices", vtx_h);

	edm::Handle<edm::ValueMap<bool> > veto_id_decisions;
	edm::Handle<edm::ValueMap<bool> > loose_id_decisions;
	edm::Handle<edm::ValueMap<bool> > medium_id_decisions;
	edm::Handle<edm::ValueMap<bool> > tight_id_decisions;
	edm::Handle<edm::ValueMap<bool> > heep_id_decisions;
	iEvent.getByToken(electronVetoIdMapToken_,veto_id_decisions);
	iEvent.getByToken(electronLooseIdMapToken_,loose_id_decisions);
	iEvent.getByToken(electronMediumIdMapToken_,medium_id_decisions);
	iEvent.getByToken(electronTightIdMapToken_,tight_id_decisions);  
	iEvent.getByToken(eleHEEPIdMapToken_, heep_id_decisions);
	const pat::MET &MEt = met->front();

	reco::VertexCollection::const_iterator firstGoodVertex = vtx_h->end();
	for (reco::VertexCollection::const_iterator it = vtx_h->begin(); it != firstGoodVertex; it++)
	{
		isGoodVertex(*it);
		firstGoodVertex = it;
		break;
	}
	// require a good vertex 
	if (firstGoodVertex == vtx_h->end()) return;
	if(isMC) {

		hFillEventTotal->Fill(1);
		//////////////////////////////////////// weights ///////////////////////////////
		edm::Handle<GenEventInfoProduct> genEvt;
		iEvent.getByLabel("generator",genEvt);
		//
		//     // event weight
		weightevt=genEvt->weight();
		sumweight = sumweight+weightevt;
		//		cout<<"weightevt"<<weightevt<<endl;
		edm::Handle<edm::View<reco::GenParticle> > pruned;  
		iEvent.getByLabel("prunedGenParticles",pruned); 
		unsigned int Ngen = pruned->size();
		indexmuTau = indexTauh = 0;
		int muhad;
		muhad=0;
		for (unsigned int iGen = 0; iGen < Ngen; iGen++)
		{
			const GenParticle& genP = (*pruned)[iGen];
			//			std::cout<<"genP" << genP.pt()<<std::endl;
			if(fabs(genP.pdgId()) == 13) {
				if(genP.isDirectPromptTauDecayProductFinalState()) {
					const reco::Candidate *Moth = genP.mother(0);

					if(fabs(Moth->mother(0)->pdgId()) == MotherpdgID){
						//const reco::GenParticle &mom = genP.mother(0);				 
						const reco::Candidate *Moth = genP.mother(0);
						if(IsFirstCopy (Moth,true)) { 
							//							std::cout<<"IsFirstCopy:"<<std::endl;
							muhad++;                
							genmuinfo.SetPtEtaPhiM(genP.pt(),genP.eta(),genP.phi(),genP.mass());
						}
					} else {
						const reco::Candidate *Moth1 = GetFirstCopy(Moth);
						if(fabs(Moth1->pdgId()) == 15) {
							if(fabs(Moth1->mother(0)->pdgId()) == MotherpdgID){
								//								std::cout<<"Getting first copy:"<<std::endl; 
								muhad++;		
								genmuinfo.SetPtEtaPhiM(genP.pt(),genP.eta(),genP.phi(),genP.mass());
							} 
						}
					}
				}
			}
		}
		//	std::cout<<"muhad:"<<muhad<<std::endl;
		hNumberOfMuons->Fill(muhad);
//		hFillEventTotal->Fill(1);

		//		if(!(muhad == 1)) return;
		int itauh=0;
		for (unsigned int iGen = 0; iGen < Ngen; iGen++)
		{
			const GenParticle& genP = (*pruned)[iGen];
			if(fabs(genP.pdgId()) == 15) {
				const reco::Candidate* genp = &genP;
				const reco::Candidate *Moth1 = GetFirstCopy(genp);
				if(fabs(Moth1->mother(0)->pdgId()) == MotherpdgID) {
					int decay = GetTauDecay (GetLastCopy(genp));
					if(decay==2) {
						//cout<<"Tau Pt"<<genp->pt()<<endl;
						const  GenParticle& Tauh = GetTauHad (genp);
						//cout<<"Had TauPt"<<Tauh.pt()<<endl;
						itauh++;
						gentauinfo.SetPtEtaPhiM(Tauh.pt(),Tauh.eta(),Tauh.phi(),Tauh.mass());

					}
				}
			}
		}
		hNumberOftaus->Fill(itauh);
		hFillEventMuon->Fill(1);
		if(GenReq) {if(!((muhad == 1) && (itauh == 1))) return;}
		if(DYOthersBG) {if((muhad == 1) && (itauh == 1)) return;} 
	} // MC boolean
	bool trigfired = false;
	const edm::TriggerNames &names = iEvent.triggerNames(*triggerBits);
	//      std::cout << "\n === TRIGGER PATHS === " << std::endl;
	for (unsigned int i = 0, n = triggerBits->size(); i < n; ++i) {
		if((names.triggerName(i) == HLTPath1) && (triggerBits->accept(i) == 1)) {trigfired = true;} 
	}
                                                                                             
	if(!(trigfired)) return;
	hFillEventTau->Fill(1); 
	// NOW TRIGGER LEVEL INFORMATION
	//  bool isData = false;
	matchedTrigger1 = false;
	matchedTrigger2 = false;

	if(!(isMC)) weightevt = 1;
	int muonid = 0;
	int muoniso = 0;
	int tauid = 0;
	int chargereq = 0;
	int mtcut=0;
	int pzeta=0;
	int muonmatch = 0;
	int taumatch1 = 0;
	int taumatch2 = 0;
	for (pat::MuonCollection::const_iterator myMuon = muons->begin(); myMuon != muons->end(); ++myMuon) {

		if((MuSelection( *myMuon, *firstGoodVertex ))){
			muonid++;


			muoniso++;
			// tau loop
			for (pat::TauCollection::const_iterator myTau = taus->begin(); myTau != taus->end(); ++myTau) {                                                                                                   
				if((TauSelection(*myTau, *firstGoodVertex))) {
					tauid++;
					//pat::PackedCandidate const* packedLeadTauCand = dynamic_cast<pat::PackedCandidate const*>(src.leadChargedHadrCand().get());
					if(deltaR(myMuon->p4(), myTau->p4()) <= 0.5) continue; 
					if(isOSCharge) {if(!(myMuon->charge() * myTau->charge() == -1)) continue;}
					if(!isOSCharge) {if(!(myMuon->charge() * myTau->charge() == 1)) continue;}
					chargereq++;
					myMT =  mTCalculation(*myMuon,MEt);
					if(isZtau) {if(myMT >= 40. ) continue; } 
					if(isZprime) {if(!(cos(TMath::Abs(normalizedPhi(myMuon->phi() - myTau->phi())))< -0.95)) continue; }
					mtcut++;
					float PZETA = PZeta(*myMuon, *myTau, MEt);
					float PZETAvis = PZetaVis(*myMuon, *myTau);
					if(isZprime) {if(PZETA - 3.1*PZETAvis < -50) continue; }
					pzeta++;
                                         
					if(BJetsinEvent(jets, *myMuon,*myTau)) continue;
					if(GenReq ) if(!(dR(genmuinfo.Eta(), genmuinfo.Phi(),myMuon->eta(),myMuon->phi()) < 0.5)) continue;
					if(GenReq ) if(!(dR(gentauinfo.Eta(), gentauinfo.Phi(),myTau->eta(),myTau->phi()) < 0.5)) continue;  
					if(ExtraMuon(*myMuon, muons, *firstGoodVertex)) continue; 
	
					if(!(( myTau->p4()+myMuon->p4()).M() > 40.)) continue;
					if(!(( myTau->p4()+myMuon->p4()).M() < 80.)) continue;
					/////////////////////////////// extra electron lepton /////////////////
					iselectron=false;
					int nElec=0;
					for(edm::View<pat::Electron>::const_iterator el = electron_pat->begin(); el != electron_pat->end(); el++) {
						nElec++;

						const Ptr<pat::Electron> elPtr(electron_pat, el - electron_pat->begin() );
						double mhits = el->gsfTrack()->hitPattern().numberOfHits(reco::HitPattern::MISSING_INNER_HITS);

						if( (deltaR(myMuon->p4(), el->p4()) > 0.5) && (deltaR(myTau->p4(), el->p4()) > 0.5) &&  (el->pt() > 10) && (fabs(el->eta()) < 2.5) && (fabs(el->gsfTrack()->dz(firstGoodVertex->position())) < 0.2) && (fabs(el->gsfTrack()->dxy(firstGoodVertex->position())) < 0.045)) { 
							double combRelIsoPF = ((el->pfIsolationVariables().sumChargedHadronPt + max(el->pfIsolationVariables().sumNeutralHadronEt +el->pfIsolationVariables().sumPhotonEt - 0.5 * el->pfIsolationVariables().sumPUPt, 0.0))/el->pt());
							if((*medium_id_decisions)[ elPtr ] == 1 && mhits <= 1) {
								if(combRelIsoPF < 0.3) {
									iselectron = true;
									break;

								}			
							}
						}
					}
					if(iselectron) continue;	

					pair<bool, TLorentzVector> muonmatched;
					muonmatched = MuonMatchingToTrigObjects(*myMuon,triggerObjects, names, HLTFilter1a,HLTPath1); 
					if(!(muonmatched.first)) continue;


					muonmatch++;
					hFillPtDen->Fill(myTau->pt(),weightevt);
					hFillEtaDen->Fill(myTau->eta(),weightevt);
					h2D_den->Fill(myTau->eta(),myTau->pt(),weightevt); 

					tauptdis->Fill(myTau->pt(),weightevt);
					tauetadis->Fill(myTau->eta(),weightevt);
					muptdis->Fill(myMuon->pt(),weightevt); 
					muetadis->Fill(myMuon->eta(),weightevt); 
					mutaumassdis->Fill((myTau->p4()+myMuon->p4()).M(),weightevt); 
					mtdis->Fill(myMT,weightevt); 
					//
					pair<vector<bool>,vector<pat::TriggerObjectStandAlone> > taumatchedNum;
					taumatchedNum  = preTauMatchingToTrigObjects(*myTau,triggerObjects, names, HLTFilter2b, HLTPath2);


					vector<bool> firstPart; vector<pat::TriggerObjectStandAlone> secondPart;
					firstPart.clear(); secondPart.clear();
					firstPart = taumatchedNum.first;
					secondPart = taumatchedNum.second;    
					if(firstPart.size() > 0) {
						//						hFillPtDen->Fill(myTau->pt(),weightevt);
						//						hFillEtaDen->Fill(myTau->eta(),weightevt);
						//						h2D_den->Fill(myTau->eta(),myTau->pt(),weightevt); 
						bool matchedit = false;
						//////////////////// fill histograms//////////////

						//						tauptdis->Fill(myTau->pt(),weightevt);
						//						tauetadis->Fill(myTau->eta(),weightevt);
						//						muptdis->Fill(myMuon->pt(),weightevt); 
						//						muetadis->Fill(myMuon->eta(),weightevt); 
						//						mutaumassdis->Fill((myTau->p4()+myMuon->p4()).M(),weightevt); 
						//						mtdis->Fill(myMT,weightevt); 
						for(unsigned int k =0; k < secondPart.size() ; k++) {

							secondPart.at(k).unpackPathNames(names);         
							for(unsigned int kl =0; kl < (secondPart.at(k).filterLabels()).size(); kl++) {                         
								if((secondPart.at(k).filterLabels())[kl].find(HLTFilter2b) != std::string::npos ) {
									matchedit=true;

								}
							}
						}
						if((matchedit)) {
							taumatch2++; 

							hFillPtNum->Fill(myTau->pt(),weightevt);
							hFillEtaNum->Fill(myTau->eta(),weightevt);                                                                     
							h2D_num->Fill(myTau->eta(),myTau->pt(),weightevt);  
						}
					}
				}
			} // if condition of tau

		}// tau loop objects


		// tau loop starts
		//

	}
	if(muonid >= 1) hmuonid->Fill(1);
	if(muoniso >= 1) hmuoniso->Fill(1);
	if(tauid >= 1) htauid->Fill(1);
	if(chargereq >= 1) hchargereq->Fill(1);
	if(mtcut>=1) hmtcut->Fill(1);
	if(pzeta>=1) hpzeta->Fill(1);
	if(muonmatch >= 1) hmuonmatch->Fill(1);
	if(taumatch1 >= 1) htaumatch1->Fill(1);
	if(taumatch2 >= 1) htaumatch2->Fill(1);
#ifdef THIS_IS_AN_EVENT_EXAMPLE
	Handle<ExampleData> pIn;
	iEvent.getByLabel("example",pIn);
#endif

#ifdef THIS_IS_AN_EVENTSETUP_EXAMPLE
	ESHandle<SetupData> pSetup;
	iSetup.get<SetupRecord>().get(pSetup);
#endif
}


// ------------ method called once each job just before starting event loop  ------------
	void 
TriggerEff::beginJob()
{

	sumweight=0.;
	myTree = fs->make<TTree>("TauTree","TauTree");
	myTree->Branch("totalweight",&totalweight);
	//hSumWeight = fs->make<THnSparseD>("hSumWeight","hSumWeight",1,100000000000000,0,100000000000000);
	double a[]= {0,20,40,60,80,100,120,140,160,180,200,240,300};
	double b[] = {-2.1,-1.566,-1.442,-1.0,-0.5,0,0.5,1.0,1.442,1.566,2.1}; 

	tauptdis = fs->make<TH1F>("tauptdis","tauptdis",500,0,500);
	tauetadis = fs->make<TH1F>("tauetadis","tauetadis",100,-5,5);
	muptdis  = fs->make<TH1F>("muptdis","muptdis",500,0,500);
	muetadis = fs->make<TH1F>("muetadis","muetadis",100,-5,5);

	mutaumassdis = fs->make<TH1F>("mutaumassdis","mutaumassdis",500,0,500);
	mtdis = fs->make<TH1F>("mtdis","mtdis",200,0,200);

	hFillPtDen = fs->make<TH1F>("hFillPtDen","hFillPtDen",500,0,500);
	hFillPtNum = fs->make<TH1F>("hFillPtNum","hFillPtNum",500,0,500);
	hFillEtaDen = fs->make<TH1F>("hFillEtaDen","hFillEtaDen",100,-5,5);
	hFillEtaNum = fs->make<TH1F>("hFillEtaNum","hFillEtaNum",100,-5,5);
	h2D_den = fs->make<TH2F>("h2D_den","h2D_den",10,b,12,a);
	h2D_num = fs->make<TH2F>("h2D_num","h2D_num",10,b,12,a);

	hNumberOfMuons = fs->make<TH1F>("hNumberOfMuons","hNumberOfMuons",10,0,10);
	hNumberOftaus = fs->make<TH1F>("hNumberOftaus","hNumberOftaus",10,0,10);
	hFillEventTotal = fs->make<TH1F>("hFillEventTotal","hFillEventTotal",2,0,2);
	hFillEventMuon = fs->make<TH1F>("hFillEventMuon","hFillEventMuon",2,0,2);
	hFillEventTau = fs->make<TH1F>("hFillEventTau","hFillEventTau",2,0,2);
	hmuonid = fs->make<TH1F>("hmuonid","hmuonid",2,0,2);
	hmuoniso = fs->make<TH1F>("hmuoniso","hmuoniso",2,0,2);
	htauid = fs->make<TH1F>("htauid","htauid",2,0,2);
	hchargereq = fs->make<TH1F>("hchargereq","hchargereq",2,0,2);
	hmtcut = fs->make<TH1F>("hmtcut","hmtcut",2,0,2);
	hpzeta = fs->make<TH1F>("hpzeta","hpzeta",2,0,2);
	hmuonmatch = fs->make<TH1F>("hmuonmatch","hmuonmatch",2,0,2);
	htaumatch1 = fs->make<TH1F>("htaumatch1","htaumatch1",2,0,2);
	htaumatch2 = fs->make<TH1F>("htaumatch2","htaumatch2",2,0,2);
	hFillPtDen->Sumw2();                 
	h2D_den->Sumw2();
	h2D_num->Sumw2();                                                                  
	hFillPtNum->Sumw2();
	hFillEtaDen->Sumw2();
	hFillEtaNum->Sumw2();
}

// ------------ method called once each job just after ending the event loop  ------------
	void 
TriggerEff::endJob() 
{
	//	hSumWeight->Fill(sumweight[0]);
	totalweight = sumweight;
	cout<<"sumweight"<<sumweight<<"totalweight"<<totalweight<<endl;
	myTree->Fill();
}

// ------------ method called when starting to processes a run  ------------

	void 
TriggerEff::beginRun(edm::Run const& iRun, edm::EventSetup const& iSetup)
{
	bool changed(true);
	hltConfig_.init(iRun,iSetup,"HLT",changed);
}


// ------------ method called when ending the processing of a run  ------------

	void 
TriggerEff::endRun(edm::Run const&, edm::EventSetup const&)
{
}


// ------------ method called when starting to processes a luminosity block  ------------
/*
   void 
   TriggerEff::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
   {
   }
   */

// ------------ method called when ending the processing of a luminosity block  ------------
/*
   void 
   TriggerEff::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
   {
   }
   */

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
TriggerEff::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
	//The following says we do not know what parameters are allowed so do no validation
	// Please change this to state exactly what you do use, even if it is no parameters
	edm::ParameterSetDescription desc;
	desc.setUnknown();
	descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(TriggerEff);

bool TriggerEff::ExtraMuon(const pat::Muon& muobject, edm::Handle<pat::MuonCollection> muons, const reco::Vertex &vtx){
	bool ismuon;   
	ismuon = false;
	for (pat::MuonCollection::const_iterator myMuon = muons->begin(); myMuon != muons->end(); ++myMuon) { 
		if(!(myMuon->pt() > 15.)) continue;
		if(fabs(myMuon->eta()) > 2.4) continue;
		if((deltaR(myMuon->p4(), muobject.p4()) < 0.5)) continue;
		if(!(muobject.charge() * myMuon->charge() < 0)) continue;
		if(!(fabs(myMuon->muonBestTrack()->dz(vtx.position())) < 0.2)) continue;
		if(!(fabs(myMuon->muonBestTrack()->dxy(vtx.position())) < 0.045)) continue;
		double iso = ((myMuon->pfIsolationR03().sumChargedHadronPt + max(myMuon->pfIsolationR03().sumNeutralHadronEt  + myMuon->pfIsolationR03().sumPhotonEt - 0.5 * myMuon->pfIsolationR03().sumPUPt, 0.0))/myMuon->pt());
		if(iso < 0.3 && myMuon->isTrackerMuon() && myMuon->isPFMuon() && myMuon->isGlobalMuon()){

			ismuon =  true;

		}  
	}
	          
	return ismuon;
}    


float TriggerEff::mTCalculation(const pat::Muon& muobject, const pat::MET& metobject){
	float mt = -1;
	float pX = muobject.px()+metobject.px();
	float pY = muobject.py()+metobject.py();
	float et = muobject.et() + TMath::Sqrt(metobject.px()*metobject.px() + metobject.py()*metobject.py());
	mt = TMath::Sqrt(et*et-(pX*pX + pY*pY));
	return mt;

}


float TriggerEff::PZetaVis(const pat::Muon& muobject, const pat::Tau& tauobject){
	float pzetavis;
	pzetavis = 999;
	float zetax = TMath::Cos(muobject.phi()) + TMath::Cos(tauobject.phi()) ;
	float zetay = TMath::Sin(muobject.phi()) + TMath::Sin(tauobject.phi()) ;
	float zetaR = TMath::Sqrt(pow(zetax,2) + pow(zetay,2));
	zetax = zetax/zetaR;
	zetay = zetay/zetaR;

	float visPx = muobject.px() + tauobject.px() ;
	float visPy = muobject.py() + tauobject.py() ;

	pzetavis = visPx*zetax + visPy*zetay;
	return pzetavis;

}
float TriggerEff::PZeta(const pat::Muon& muobject, const pat::Tau& tauobject, const pat::MET& metobject){
	float pzeta;
	pzeta = 999;
	float zetax = TMath::Cos(muobject.phi()) + TMath::Cos(tauobject.phi()) ;
	float zetay = TMath::Sin(muobject.phi()) + TMath::Sin(tauobject.phi()) ;
	float zetaR = TMath::Sqrt(pow(zetax,2) + pow(zetay,2));
	zetax = zetax/zetaR;
	zetay = zetay/zetaR;

	float vPx = muobject.px() + tauobject.px()+metobject.px() ;
	float vPy = muobject.py() + tauobject.py()+metobject.py() ;

	pzeta = vPx*zetax + vPy*zetay;
	return pzeta;


}

const reco::Candidate* TriggerEff::GetLastCopy (const reco::Candidate* part)
{

	int cloneInd = -1;
	int id = part->pdgId();
	for (unsigned int iDau = 0; iDau < part->numberOfDaughters(); iDau++)
	{
		const reco::Candidate * Dau = part->daughter( iDau );
		if (id == Dau->pdgId())
		{
			cloneInd = iDau;
			break;
		}
	}

	if (cloneInd == -1) return part;
	else return (GetLastCopy (part->daughter(cloneInd)));

}



int TriggerEff::GetTauDecay (const reco::Candidate* part)
{
	if (abs(part->pdgId()) != 15) return -1; // only on taus
	int decay = -1;
	int nele = 0;
	int nmu = 0;
	for (unsigned int iDau = 0; iDau < part->numberOfDaughters(); iDau++)
	{
		const reco::Candidate * Dau = part->daughter(iDau);
		int dauId = abs(Dau->pdgId());
		if (dauId == 11) nele++;
		if (dauId == 13) nmu++;
	}

	if (nmu == 1 && nele == 0) decay = 0;
	if (nmu == 0 && nele == 1) decay = 1;
	if (nmu == 0 && nele == 0) decay = 2;

	return decay; // -1 if strange things happen
}


int TriggerEff::GetTauDecay (const reco::GenParticle& part)
{
	const reco::Candidate* p = &part;
	return TriggerEff::GetTauDecay(p);
}

bool TriggerEff::IsFirstCopy (const reco::Candidate* part, const bool checkAbsPdg)
{
	bool isFirst = true;
	int thisPdgId = part->pdgId();
	for (unsigned int iMo = 0; iMo < part->numberOfMothers(); iMo++)
	{
		const reco::Candidate * Mo = part->mother(iMo);
		bool pdgMatch = (checkAbsPdg ? (abs(thisPdgId) == abs(Mo->pdgId())) : (thisPdgId == Mo->pdgId()) );
		if (pdgMatch)
		{
			isFirst = false;
			break;
		}
	}
	return isFirst;
}

const reco::Candidate* TriggerEff::GetFirstCopy (const reco::Candidate* part)
{
	int cloneInd = -1;
	int id = part->pdgId();
	for (unsigned int iMot = 0; iMot < part->numberOfMothers(); iMot++)
	{
		const reco::Candidate * Mot = part->mother( iMot );
		if (id == Mot->pdgId())
		{
			cloneInd = iMot;
			break;
		}
	}

	if (cloneInd == -1) return part;
	else return (GetFirstCopy (part->mother(cloneInd)));

}


reco::GenParticle TriggerEff::GetTauHad (const reco::Candidate* part)
{

	reco::Candidate::LorentzVector p4Had (0,0,0,0);
	for (unsigned int iDau = 0; iDau < part->numberOfDaughters(); iDau++)
	{
		const reco::Candidate * Dau = part->daughter( iDau );
		int dauId = abs(Dau->pdgId());
		if (dauId != 12 && dauId != 14 && dauId != 16 && dauId != 11 && dauId != 13) // no neutrinos
			p4Had += Dau->p4();
	}

	int sign = part->pdgId() / abs(part->pdgId());
	reco::GenParticle TauH = reco::GenParticle (part->charge(), p4Had, part->vertex(), sign*66615, part->status(), true);
	return TauH;
}


bool TriggerEff::isGoodVertex(const reco::Vertex& vtx)
{      
	if (vtx.isFake()) return false;
	if (vtx.ndof() < 4.) return false;
	if (vtx.position().Rho() > 2.) return false;
	if (fabs(vtx.position().Z()) > 24) return false;
	return true;
}   

pair<bool, TLorentzVector> TriggerEff::MuonMatchingToTrigObjects(const pat::Muon& myMuon, edm::Handle<pat::TriggerObjectStandAloneCollection> triggerObjects,const edm::TriggerNames &names, string filterName, string pathname) { 
	//	pair<bool, TLorentzVector> mymuonobject;
	bool ismatched = false;
	TLorentzVector matchedtriggerObject(0,0,0,0);
	vector<pat::TriggerObjectStandAlone> ObjMatchedMuon;
	ObjMatchedMuon.clear();
	for(pat::TriggerObjectStandAloneCollection::const_iterator it = triggerObjects->begin() ; it !=triggerObjects->end() ; it++){
		pat::TriggerObjectStandAlone *aObj = const_cast<pat::TriggerObjectStandAlone*>(&(*it));
		aObj->unpackPathNames(names);
		if(aObj->pt() < 18) continue;
		if(!(aObj->hasTriggerObjectType(trigger::TriggerMuon))) continue;
		for(unsigned int k =0; k < (aObj->filterLabels()).size() ; k++){
			if(deltaR( aObj->triggerObject().p4(), myMuon.p4() ) < 0.3 &&((aObj->filterLabels())[k].find(filterName) != std::string::npos )) {
				ObjMatchedMuon.push_back(*aObj);
			}	
		}
	}


	for(unsigned int k =0; k < ObjMatchedMuon.size() ; k++) {

		ObjMatchedMuon.at(k).unpackPathNames(names);	     
		for(unsigned int kl =0; kl < (ObjMatchedMuon.at(k).pathNames()).size(); kl++) {				
			if(deltaR( ObjMatchedMuon.at(k).triggerObject().p4(), myMuon.p4() ) < 0.3 && ((ObjMatchedMuon.at(k).pathNames())[kl].find(pathname) != std::string::npos )) {
				ismatched = true;
				matchedtriggerObject.SetPxPyPzE(ObjMatchedMuon.at(k).triggerObject().px(),ObjMatchedMuon.at(k).triggerObject().py(),ObjMatchedMuon.at(k).triggerObject().pz(),ObjMatchedMuon.at(k).triggerObject().energy());
			}
		}
	}
	pair<bool, TLorentzVector> mymuonobject(ismatched,matchedtriggerObject);
	return mymuonobject;
}


pair<bool, TLorentzVector> TriggerEff::TauMatchingToTrigObjects(const pat::Tau& myTau, edm::Handle<pat::TriggerObjectStandAloneCollection> triggerObjects,const edm::TriggerNames &names, string filterName, string pathname, bool numerator) { 
	//      pair<bool, TLorentzVector> mymuonobject;
	bool ismatched = false;
	TLorentzVector matchedtriggerObject(0,0,0,0);
	vector<pat::TriggerObjectStandAlone> ObjMatchedTau;
	ObjMatchedTau.clear();
	for(pat::TriggerObjectStandAloneCollection::const_iterator it = triggerObjects->begin() ; it !=triggerObjects->end() ; it++){
		pat::TriggerObjectStandAlone *aObj = const_cast<pat::TriggerObjectStandAlone*>(&(*it));
		aObj->unpackPathNames(names);
		for(unsigned int k =0; k < (aObj->filterLabels()).size() ; k++){
			if(deltaR( aObj->triggerObject().p4(), myTau.p4() ) < 0.3 &&((aObj->filterLabels())[k].find(filterName) != std::string::npos )) {
				ObjMatchedTau.push_back(*aObj);
			}       
		}
	}

	if(numerator) {
		for(unsigned int k =0; k < ObjMatchedTau.size() ; k++) {

			ObjMatchedTau.at(k).unpackPathNames(names);
			for(unsigned int kl =0; kl < (ObjMatchedTau.at(k).pathNames()).size(); kl++) {
				if(deltaR( ObjMatchedTau.at(k).triggerObject().p4(), myTau.p4() ) < 0.3 ) {
					ismatched = true;
					matchedtriggerObject.SetPxPyPzE(ObjMatchedTau.at(k).triggerObject().px(),ObjMatchedTau.at(k).triggerObject().py(),ObjMatchedTau.at(k).triggerObject().pz(),ObjMatchedTau.at(k).triggerObject().energy());
				}
			}
		}      




	} else {
		for(unsigned int k =0; k < ObjMatchedTau.size() ; k++) {

			ObjMatchedTau.at(k).unpackPathNames(names);         
			for(unsigned int kl =0; kl < (ObjMatchedTau.at(k).pathNames()).size(); kl++) {                    
				if(deltaR( ObjMatchedTau.at(k).triggerObject().p4(), myTau.p4() ) < 0.3 && ((ObjMatchedTau.at(k).pathNames())[kl].find(pathname) != std::string::npos )) {
					ismatched = true;
					matchedtriggerObject.SetPxPyPzE(ObjMatchedTau.at(k).triggerObject().px(),ObjMatchedTau.at(k).triggerObject().py(),ObjMatchedTau.at(k).triggerObject().pz(),ObjMatchedTau.at(k).triggerObject().energy());
				}
			}
		}
	}
	pair<bool, TLorentzVector> mytauobject(ismatched,matchedtriggerObject);        

	return mytauobject;
}


float TriggerEff::deltaPhi( float a, float b) {
	float result = a-b;
	while (result > M_PI) result -= 2* M_PI;
	while (result <= -M_PI) result += 2* M_PI;
	return (fabs(result));

}        

float TriggerEff::dR(float l1eta, float l1phi, float l2eta, float l2phi ) {
	float deta = l1eta - l2eta;
	float dphi = deltaPhi(l1phi,l2phi);
	return sqrt(deta*deta + dphi*dphi);
}    

bool TriggerEff::TauSelection( const pat::Tau &myTau, const reco::Vertex &vtx){

	if(myTau.pt() <= TauPtCut ) return false;
	if(fabs(myTau.eta()) >= TauEtaCut ) return false;
	if(myTau.tauID(TauDMF) < 0.5) return false;
	if(myTau.tauID(TauMuVeto) < 0.5) return false;
	if(myTau.tauID(TauEleVeto) < 0.5 ) return false;
	if(myTau.tauID(TauIsoString) < 0.5 ) return false;
//	if(myTau.tauID(TauIsoString) >= TauIsoCutMax) return false;
//	if(myTau.tauID(TauIsoString) < TauIsoCutMin) return false; 
	if(!(myTau.leadChargedHadrCand().isNonnull() && myTau.leadChargedHadrCand()->pt() > 5.)) return false;
	pat::PackedCandidate const* packedLeadTauCand = dynamic_cast<pat::PackedCandidate const*>(myTau.leadChargedHadrCand().get());
	if(!(fabs(packedLeadTauCand->dz(vtx.position())) < 0.2)) return false;
	if(!(fabs(packedLeadTauCand->dxy(vtx.position())) < 0.045)) return false;
	return true;
} 


bool TriggerEff::MuSelection( const pat::Muon &myMuon, const reco::Vertex &vtx ) {

	if(myMuon.pt() <= MuonPtCut) return false;
	if(fabs(myMuon.eta()) >=MuonEtaCut) return false;
	if(!(myMuon.isMediumMuon())) return false;
	if(!(myMuon.isPFMuon())) return false;
	if(!(fabs(myMuon.muonBestTrack()->dz(vtx.position())) < 0.2)) return false;
	if(!(fabs(myMuon.muonBestTrack()->dxy(vtx.position())) < 0.045)) return false;
	double iso = ((myMuon.pfIsolationR03().sumChargedHadronPt + max(myMuon.pfIsolationR03().sumNeutralHadronEt  + myMuon.pfIsolationR03().sumPhotonEt - 0.5 * myMuon.pfIsolationR03().sumPUPt, 0.0))/myMuon.pt());

	if(iso >= IsoMuonMax)  return false;
	return true; 

}

pair<vector<bool>,vector<pat::TriggerObjectStandAlone>> TriggerEff::preTauMatchingToTrigObjects(const pat::Tau& myTau, edm::Handle<pat::TriggerObjectStandAloneCollection> triggerObjects, const edm::TriggerNames &names,string filterName, string pathname) {
	vector<bool> ismatchedvector;
	ismatchedvector.clear();
	TLorentzVector matchedtriggerObject(0,0,0,0);
	vector<pat::TriggerObjectStandAlone> ObjMatchedTau;
	ObjMatchedTau.clear(); 
	for(pat::TriggerObjectStandAloneCollection::const_iterator it = triggerObjects->begin() ; it !=triggerObjects->end() ; it++){
		pat::TriggerObjectStandAlone *aObj = const_cast<pat::TriggerObjectStandAlone*>(&(*it));
		aObj->unpackPathNames(names);
		if(!(aObj->hasTriggerObjectType(trigger::TriggerTau))) continue;

		for(unsigned int k =0; k < (aObj->filterLabels()).size() ; k++){
			if(deltaR( aObj->triggerObject().p4(), myTau.p4() ) < 0.3 ) {
				ObjMatchedTau.push_back(*aObj);
				ismatchedvector.push_back(true);
			}       
		}        
	}                
	pair<vector<bool>, vector<pat::TriggerObjectStandAlone>> mytauobject(ismatchedvector,ObjMatchedTau);
	return mytauobject;    



}

bool TriggerEff::BJetsinEvent( edm::Handle<pat::JetCollection> myjets, const pat::Muon &myMuon, const pat::Tau &myTau){
	bool isbjet = false;
	for(pat::JetCollection::const_iterator ijet = myjets->begin() ; ijet != myjets->end() ; ijet++){

		double _bdiscr2 = ijet->bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags");
		//PF jet ID
		float NHF = ijet->neutralHadronEnergyFraction();
		float NEMF = ijet->neutralEmEnergyFraction();
		float CHF = ijet->chargedHadronEnergyFraction();
		float MUF = ijet->muonEnergyFraction();
		float CEMF = ijet->chargedEmEnergyFraction();
		int NumNeutralParticles =ijet->neutralMultiplicity();
		int NumConst = ijet->chargedMultiplicity()+NumNeutralParticles;
		float CHM = ijet->chargedMultiplicity();
		float absjeta = fabs(ijet->eta());

		bool jetid1;

		jetid1 = false;

		if(absjeta<=3.0){

			if( (NHF<0.90 && NEMF<0.90 && NumConst>1) && ((absjeta<=2.4 && CHF>0 && CHM>0 && CEMF<0.99) || absjeta>2.4) ) {

				if( (NHF<0.90 && NEMF<0.90 && NumConst>1 && MUF<0.8) && ((absjeta<=2.4 && CHF>0 && CHM>0 && CEMF<0.90) || absjeta>2.4)){

					jetid1 = true;

				}

			}

		}else{

			if(NEMF<0.90 && NumNeutralParticles>10 ){

				jetid1 = true;      

			}

		}   



		if(ijet->pt() > 30. && fabs(ijet->eta()) < 2.4 &&( jetid1 ) && _bdiscr2 > 0.890 && (deltaR(ijet->p4(),myTau.p4()) > 0.5) && (deltaR(ijet->p4(),myMuon.p4()) > 0.5)) {

			isbjet = true;
		}
	}
	return isbjet;
}
