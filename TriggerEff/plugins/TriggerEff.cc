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

		//virtual void beginRun(edm::Run const&, edm::EventSetup const&) override;
		//virtual void endRun(edm::Run const&, edm::EventSetup const&) override;
		//virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
		//virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;

		// ----------member data ---------------------------

		string HLTPath1, HLTFilter1;
		string HLTPath2, HLTFilter2a,HLTFilter2b,HLTFilter2c,HLTFilter2d;                                                        
		string HLTPath3 , HLTFilter3a, HLTFilter3b;	
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

		bool matchedTrigger1, matchedTrigger2 ;
		// matched trigger objects
		pair<bool, TLorentzVector> MuonMatchingToTrigObjects(const pat::Muon& myMuon, edm::Handle<pat::TriggerObjectStandAloneCollection> triggerObjects, const edm::TriggerNames &names,string filterName, string pathname);                 
		pair<bool, TLorentzVector> TauMatchingToTrigObjects(const pat::Tau& myTau, edm::Handle<pat::TriggerObjectStandAloneCollection> triggerObjects, const edm::TriggerNames &names,string filterName, string pathname, bool numerator); //numerator boolean is added to not check the trigger path  
		float deltaPhi( float a, float b);
		float dR(float l1eta, float l1phi, float l2eta, float l2phi );       
		edm::Service<TFileService> fs;
		TH1F *hFillPtDen;
		TH1F *hFillPtNum;
		TH1F *hFillEtaDen;
		TH1F *hFillEtaNum; 
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
TriggerEff::TriggerEff(const edm::ParameterSet& iConfig)

{
	//now do what ever initialization is needed
	triggerResultsTag_  = iConfig.getParameter<edm::InputTag>("triggerResults");
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
	edm::Handle<edm::View<reco::GenParticle> > pruned;
	edm::Handle<pat::MuonCollection> muons;
	edm::Handle<pat::TauCollection> taus;
	edm::Handle<pat::METCollection> met;

	iEvent.getByLabel("slimmedMuons",muons);
	iEvent.getByLabel("slimmedTaus",taus);
	iEvent.getByLabel("slimmedMETs",met);

	iEvent.getByLabel(triggerResultsTag_, triggerBits);
	iEvent.getByLabel("selectedPatTrigger", triggerObjects);
	iEvent.getByLabel("patTrigger", triggerPrescales);
	iEvent.getByLabel("prunedGenParticles",pruned);
	edm::Handle<reco::VertexCollection> vtx_h;
	iEvent.getByLabel("offlineSlimmedPrimaryVertices", vtx_h);
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

	const edm::TriggerNames &names = iEvent.triggerNames(*triggerBits);
//	std::cout << "\n === TRIGGER PATHS === " << std::endl;
//	for (unsigned int i = 0, n = triggerBits->size(); i < n; ++i) {
//	}


	bool isMC = true;
	if(isMC) {

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
					//const reco::GenParticle &mom = genP.mother(0);				 
					const reco::Candidate *Moth = genP.mother(0);
					if(IsFirstCopy (Moth,15)) { 
						muhad++;                
						genmuinfo.SetPtEtaPhiM(genP.pt(),genP.eta(),genP.phi(),genP.mass());

					} else {
						const reco::Candidate *Moth1 = GetFirstCopy(Moth);
						if(fabs(Moth1->pdgId()) == 15) {
							muhad++;		
							genmuinfo.SetPtEtaPhiM(genP.pt(),genP.eta(),genP.phi(),genP.mass());
						} 
					}
				}
			}
		}
		if(!(muhad == 1)) return;
		int itauh=0;
		for (unsigned int iGen = 0; iGen < Ngen; iGen++)
		{
			const GenParticle& genP = (*pruned)[iGen];
			if(fabs(genP.pdgId()) == 15) {
				const reco::Candidate* genp = &genP;
				const reco::Candidate *Moth1 = GetFirstCopy(genp);
				if(fabs(Moth1->mother(0)->pdgId()) == 23) {
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
//		if(!(itauh == 1)) return;
	}



	// NOW TRIGGER LEVEL INFORMATION

	bool isData = false;
	matchedTrigger1 = false;
	matchedTrigger2 = false;

	if(isMC){

		HLTPath1 = "HLT_IsoMu24_eta2p1_v1" ;
		HLTFilter1= "hltL3crIsoL1sMu20Eta2p1L1f0L2f10QL3f24QL3trkIsoFiltered0p09" ;

		HLTPath2 = "HLT_IsoMu17_eta2p1_LooseIsoPFTau20_v1";
		HLTFilter2a= "hltL3crIsoL1sMu16erTauJet20erL1f0L2f10QL3f17QL3trkIsoFiltered0p09";
		HLTFilter2b="hltOverlapFilterIsoMu17LooseIsoPFTau20";
		HLTFilter2c= "hltPFTau20TrackLooseIsoAgainstMuon";
		HLTFilter2d= "hltOverlapFilterIsoMu17LooseIsoPFTau20";

		HLTPath3="HLT_DoubleMediumIsoPFTau40_Trk1_eta2p1_Reg_v1";
		HLTFilter3a="hltDoublePFTau40TrackPt1MediumIsolationDz02Reg";
		HLTFilter3b="hltDoublePFTau40TrackPt1MediumIsolationReg";

	}
	else if(isData) {


		HLTPath1 = "HLT_IsoMu24_eta2p1_v2" ;                         
		HLTFilter1= "hltL3crIsoL1sMu20Eta2p1L1f0L2f10QL3f24QL3trkIsoFiltered0p09" ;      

		HLTPath2 = "HLT_IsoMu17_eta2p1_LooseIsoPFTau20_v2" ;          
		HLTFilter2a= "hltL3crIsoL1sMu16erTauJet20erL1f0L2f10QL3f17QL3trkIsoFiltered0p09";
		HLTFilter2b="hltOverlapFilterIsoMu17LooseIsoPFTau20";        
		HLTFilter2c= "hltPFTau20TrackLooseIsoAgainstMuon";           
		HLTFilter2d= "hltOverlapFilterIsoMu17LooseIsoPFTau20";       

		HLTPath3="HLT_DoubleMediumIsoPFTau40_Trk1_eta2p1_Reg_v2";   
		HLTFilter3a="hltDoublePFTau40TrackPt1MediumIsolationDz02Reg";
		HLTFilter3b="hltDoublePFTau40TrackPt1MediumIsolationReg";


	}

	for (pat::MuonCollection::const_iterator myMuon = muons->begin(); myMuon != muons->end(); ++myMuon) {

		if(myMuon->pt() > 17. && fabs(myMuon->eta()) < 2.1 && myMuon->isTightMuon(*firstGoodVertex)) { 
			float iso = ((myMuon->pfIsolationR03().sumChargedHadronPt + max(myMuon->pfIsolationR03().sumNeutralHadronEt  + myMuon->pfIsolationR03().sumPhotonEt - 0.5 * myMuon->pfIsolationR03().sumPUPt, 0.0))/myMuon->pt());


			if(iso >= 0.1)  continue;
			// tau loop
			for (pat::TauCollection::const_iterator myTau = taus->begin(); myTau != taus->end(); ++myTau) {                                                                                                   

				if(myTau->pt() > 20. && fabs(myTau->eta()) < 2.3 && (myTau->tauID("decayModeFinding") > 0.5) && (myTau->tauID("againstMuonTight3") > 0.5) && (myTau->tauID("againstElectronVLooseMVA5") > 0.5 ) && (myTau->tauID("byCombinedIsolationDeltaBetaCorrRaw3Hits") < 1.5))  { 

					//pat::PackedCandidate const* packedLeadTauCand = dynamic_cast<pat::PackedCandidate const*>(src.leadChargedHadrCand().get());
					if(deltaR(myMuon->p4(), myTau->p4()) <= 0.5) continue; 
					if(!(myMuon->charge() * myTau->charge() == -1)) continue;
					float myMT =  mTCalculation(*myMuon,MEt);
					if(myMT >= 40. ) continue; 
					float PZETA = PZeta(*myMuon, *myTau, MEt);
					float PZETAvis = PZetaVis(*myMuon, *myTau);
					if(PZETA - 3.1*PZETAvis < -50) continue; 
					if(!(dR(genmuinfo.Eta(), genmuinfo.Phi(),myMuon->eta(),myMuon->phi()) < 0.5)) continue;
					if(!(dR(gentauinfo.Eta(), gentauinfo.Phi(),myTau->eta(),myTau->phi()) < 0.5)) continue;   
					pair<bool, TLorentzVector> muonmatched;
					muonmatched = MuonMatchingToTrigObjects(*myMuon,triggerObjects, names, HLTFilter2a, HLTPath2); 
					if(!(muonmatched.first)) continue;
					//cout<<"got a muon"<<endl; 
					pair<bool, TLorentzVector> taumatched; 
					taumatched  = TauMatchingToTrigObjects(*myTau,triggerObjects, names, HLTFilter2c, HLTPath2,false); 
					if(!(taumatched.first)) continue;
					//cout<<"got a tau"<<endl;
					hFillPtDen->Fill(myTau->pt());
					hFillEtaDen->Fill(myTau->eta());
					pair<bool, TLorentzVector> taumatchedNum;
					taumatchedNum  = TauMatchingToTrigObjects(*myTau,triggerObjects, names, HLTFilter3b, HLTPath2,true);
//					//cout<<"result of numertaor"<<taumatchedNum.first<<endl;  
                                        if(!(taumatchedNum.first)) continue;
					hFillPtNum->Fill(myTau->pt());
					hFillEtaNum->Fill(myTau->eta());                                                                      

				} // if condition of tau

			}// tau loop objects


			// tau loop starts
			//

		}
	}
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
	hFillPtDen = fs->make<TH1F>("hFillPtDen","hFillPtDen",500,0,500);
	hFillPtNum = fs->make<TH1F>("hFillPtNum","hFillPtNum",500,0,500);
	hFillEtaDen = fs->make<TH1F>("hFillEtaDen","hFillEtaDen",100,-5,5);
	hFillEtaNum = fs->make<TH1F>("hFillEtaNum","hFillEtaNum",100,-5,5);

}

// ------------ method called once each job just after ending the event loop  ------------
	void 
TriggerEff::endJob() 
{
}

// ------------ method called when starting to processes a run  ------------
/*
   void 
   TriggerEff::beginRun(edm::Run const&, edm::EventSetup const&)
   {
   }
   */

// ------------ method called when ending the processing of a run  ------------
/*
   void 
   TriggerEff::endRun(edm::Run const&, edm::EventSetup const&)
   {
   }
   */

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
