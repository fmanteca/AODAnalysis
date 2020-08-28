//  ______                          _____          _                  _             
//  |  _  \                        |_   _|        | |                | |            
//  | | | | ___   ___  _ __    _ __  | |    _ __  | |_  _   _  _ __  | |  ___  _ __ 
//  | | | |/ _ \ / _ \| '_ \  | '_ \ | |   | '_ \ | __|| | | || '_ \ | | / _ \| '__|
//  | |/ /|  __/|  __/| |_) | | |_) || |   | | | || |_ | |_| || |_) || ||  __/| |   
//  |___/  \___| \___|| .__/  | .__/ \_/   |_| |_| \__| \__,_|| .__/ |_| \___||_|   
//                    | |     | |                             | |                   
//                    |_|     |_|                             |_|                   



#include <memory>
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "RecoMuon/TrackingTools/interface/MuonPatternRecoDumper.h"
#include "RecoMuon/TrackingTools/interface/MuonServiceProxy.h"

#include "Geometry/Records/interface/MuonGeometryRecord.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/MuonDetId/interface/MuonSubdetId.h"
#include "DataFormats/MuonDetId/interface/DTWireId.h"
#include "DataFormats/MuonDetId/interface/CSCDetId.h"
#include "DataFormats/MuonDetId/interface/RPCDetId.h"
#include "DataFormats/SiPixelDetId/interface/PXBDetId.h"
#include "DataFormats/SiPixelDetId/interface/PXFDetId.h"
#include "DataFormats/SiStripDetId/interface/SiStripDetId.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/GenRunInfoProduct.h"
#include "DataFormats/MuonReco/interface/MuonShower.h"
#include "RecoMuon/MuonIdentification/interface/MuonShowerInformationFiller.h"
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Photon.h"
#include "DataFormats/PatCandidates/interface/IsolatedTrack.h"
#include "DataFormats/TrackReco/interface/HitPattern.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/Common/interface/RefToBase.h"
#include "DataFormats/TrackReco/interface/HitPattern.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/PatCandidates/interface/PFIsolation.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "Geometry/CommonDetUnit/interface/GeomDet.h"
#include "Geometry/CommonDetUnit/interface/GeomDetEnumerators.h"

#include "Geometry/CommonDetUnit/interface/GlobalTrackingGeometry.h"
#include "RecoMuon/DetLayers/interface/MuonDetLayerGeometry.h"
#include "RecoMuon/TransientTrackingRecHit/interface/MuonTransientTrackingRecHitBuilder.h"
#include "RecoMuon/TransientTrackingRecHit/interface/MuonTransientTrackingRecHit.h"
#include "RecoMuon/Records/interface/MuonRecoGeometryRecord.h"
#include "RecoMuon/GlobalTrackingTools/interface/StateSegmentMatcher.h"

#include "Geometry/CommonDetUnit/interface/GeomDetType.h"
#include "Geometry/CommonDetUnit/interface/GeomDet.h"
#include "Geometry/CSCGeometry/interface/CSCGeometry.h"
#include "Geometry/DTGeometry/interface/DTGeometry.h"
#include "Geometry/RPCGeometry/interface/RPCGeometry.h"

#include "MagneticField/Engine/interface/MagneticField.h"
#include "RecoMuon/TrackingTools/interface/MuonServiceProxy.h"
#include "TrackingTools/GeomPropagators/interface/Propagator.h"
#include "TrackingTools/Records/interface/TrackingComponentsRecord.h"

#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "TrackingTools/PatternTools/interface/Trajectory.h"
#include "TrackingTools/TrackFitters/interface/TrajectoryFitter.h"
#include "FWCore/Framework/interface/ConsumesCollector.h"

#include <string>
#include <iostream>
#include <vector>
#include <algorithm>

#include "TLorentzVector.h"
#include "TTree.h"
#include "TFile.h"
#include "TMath.h"

//=======================================================================================================================================================================================================================//

///////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////// FUNCTIONS ///////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////

std::string getdetid(std::string subdet, int id1, int id2, int id3){

  // Returns the detector element in the standard CMS nomenclature
  // if DT: id1 = wheel, id2 = station, id3 = sector
  // if CSC: id1 = endcap, id2 = disk, id3 = ring number
  
  std::string id;
  if(subdet == "DT"){
    id = "MB" + std::to_string(id1) + "/" +  std::to_string(id2) + "/" + std::to_string(id3);
  }else if(subdet == "CSC"){
    if(id1 == 1){
      id = "ME+" +  std::to_string(id2) + "/" + std::to_string(id3);
    }else if(id1==2){
      id = "ME-" +  std::to_string(id2) + "/" + std::to_string(id3);
    }
  }

  return id;
}

float dist3d(GlobalPoint gp1, GlobalPoint gp2){

  // Returns the distance between two points (x1,y1,z1) and (x2,y2,z2) in the 3d space

  return std::sqrt(std::pow((gp1.x()-gp2.x()),2) + std::pow((gp1.y()-gp2.y()),2) + std::pow((gp1.z()-gp2.z()),2));

}

float dist2d_xz(GlobalPoint gp1, GlobalPoint gp2){

  return std::sqrt(std::pow((gp1.x()-gp2.x()),2) + std::pow((gp1.z()-gp2.z()),2));

}

float dist2d_yz(GlobalPoint gp1, GlobalPoint gp2){

  return std::sqrt(std::pow((gp1.y()-gp2.y()),2) + std::pow((gp1.z()-gp2.z()),2));

}

float dist2d_xy(GlobalPoint gp1, GlobalPoint gp2){

  return std::sqrt(std::pow((gp1.x()-gp2.x()),2) + std::pow((gp1.y()-gp2.y()),2));

}

/////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////// DATA DEFINITION //////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////// BRANCHES /////////////////////////////////////

//-> EVENT INFO

Int_t Event_id;
Int_t Event_luminosityBlock;
Int_t Event_run;
Int_t Event_nMuons;

//-> MUON INFO

std::vector<float> Muon_Genpt;
std::vector<float> Muon_GlbTrack_pt;
std::vector<float> Muon_GlbTrack_eta;
std::vector<float> Muon_GlbTrack_phi;
std::vector<float> Muon_GlbTrack_charge;
std::vector<float> Muon_GlbTrack_ptErr;
std::vector<float> Muon_GlbTrack_Chindf;
std::vector<float> Muon_InnerTrack_pt;
std::vector<float> Muon_InnerTrack_eta;
std::vector<float> Muon_InnerTrack_phi;
std::vector<float> Muon_InnerTrack_charge;
std::vector<float> Muon_InnerTrack_ptErr;
std::vector<float> Muon_InnerTrack_Chindf;
std::vector<float> Muon_TunePTrack_pt;
std::vector<float> Muon_TunePTrack_eta;
std::vector<float> Muon_TunePTrack_phi;
std::vector<float> Muon_TunePTrack_charge;
std::vector<float> Muon_TunePTrack_ptErr;
std::vector<float> Muon_TunePTrack_Chindf;
std::vector<int> Muon_Muonid;
std::vector<int> Muon_Eventid;
std::vector<int> Muon_EventluminosityBlock;
std::vector<int> Muon_nGeomDets;
std::vector<int> Muon_nHits;
std::vector<int> Muon_nShowers;
std::vector<int> Muon_hasShowerInStation_DT_1;
std::vector<int> Muon_hasShowerInStation_DT_2;
std::vector<int> Muon_hasShowerInStation_DT_3;
std::vector<int> Muon_hasShowerInStation_DT_4;
std::vector<int> Muon_hasShowerInStation_CSC_1;
std::vector<int> Muon_hasShowerInStation_CSC_2;
std::vector<int> Muon_hasShowerInStation_CSC_3;
std::vector<int> Muon_hasShowerInStation_CSC_4;
std::vector<int> Muon_nDigisInStation_DT_1;
std::vector<int> Muon_nDigisInStation_DT_2;
std::vector<int> Muon_nDigisInStation_DT_3;
std::vector<int> Muon_nDigisInStation_DT_4;
std::vector<int> Muon_nDigisInStation_CSC_1;
std::vector<int> Muon_nDigisInStation_CSC_2;
std::vector<int> Muon_nDigisInStation_CSC_3;
std::vector<int> Muon_nDigisInStation_CSC_4;

//-> PROPAGATION INFO

std::vector<float> Prop_x;
std::vector<float> Prop_y;
std::vector<float> Prop_z;
std::vector<unsigned int> Prop_Detid;
std::vector<int> Prop_Muonid;
std::vector<int> Prop_Eventid;
std::vector<int> Prop_EventluminosityBlock;
std::vector<std::string> Prop_DetElement;
std::vector<int> Prop_isDT;
std::vector<int> Prop_isCSC;
std::vector<int> Prop_DTstation;
std::vector<int> Prop_CSCstation;

//-> HITS INFO

std::vector<float> Hit_x;
std::vector<float> Hit_y;
std::vector<float> Hit_z;
std::vector<float> Hit_distToProp;
std::vector<std::string> Hit_DetElement;
std::vector<unsigned int> Hit_Detid;
std::vector<int> Hit_Hitid;
std::vector<int> Hit_Muonid;
std::vector<int> Hit_Eventid;
std::vector<int> Hit_EventluminosityBlock;
std::vector<int> Hit_isDT;
std::vector<int> Hit_isCSC;
std::vector<int> Hit_DTstation;
std::vector<int> Hit_CSCstation;
std::vector<float> Hit_Compatibility;
std::vector<float> Hit_dirx;
std::vector<float> Hit_diry;
std::vector<float> Hit_dirz;
std::vector<float> Hit_chi2;
std::vector<int> Hit_ndof;

/////////////////////////////////////// OUTPUT //////////////////////////////////////

TFile *file_out;
TTree *tree_out;

class MuonServiceProxy;

//=======================================================================================================================================================================================================================//
class RECOAnalysis : public edm::one::EDAnalyzer<edm::one::SharedResources>  {

   public:
      explicit RECOAnalysis(const edm::ParameterSet&);
      ~RECOAnalysis();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

   private:
      virtual void beginJob() override;
      virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
      virtual void endJob() override;

      std::string output_filename;

      edm::InputTag theDTSegmentLabel, theCSCSegmentLabel,inputMuonShowerInformationValueMap;
      edm::EDGetTokenT<DTRecSegment4DCollection> dtSegmentsToken;
      edm::EDGetTokenT<CSCSegmentCollection> cscSegmentsToken;
      edm::EDGetTokenT<edm::ValueMap<reco::MuonShower> > inputMuonShowerInformationValueMapToken;
      edm::EDGetTokenT<edm::View<reco::Muon> > theMuonCollection;
      edm::EDGetTokenT<edm::View<reco::GenParticle> > theGenParticleCollection;

      std::string propagator_;
      edm::ParameterSet parameters;
      MuonServiceProxy *theService;
  
};

//=======================================================================================================================================================================================================================//
RECOAnalysis::RECOAnalysis(const edm::ParameterSet& iConfig)
{

   usesResource("TFileService");
   
   parameters = iConfig;

   theMuonCollection = consumes<edm::View<reco::Muon> >  (parameters.getParameter<edm::InputTag>("MuonCollection"));
   theGenParticleCollection = consumes<edm::View<reco::GenParticle> >  (parameters.getParameter<edm::InputTag>("genParticleCollection"));

   propagator_ = iConfig.getParameter<std::string>("Propagator");
   
   edm::ParameterSet serviceParameters = iConfig.getParameter<edm::ParameterSet>("ServiceParameters");
   theService = new MuonServiceProxy(serviceParameters);

   theDTSegmentLabel = iConfig.getParameter<edm::InputTag>("segmentsDt");
   theCSCSegmentLabel = iConfig.getParameter<edm::InputTag>("segmentsCSC");
   dtSegmentsToken = consumes<DTRecSegment4DCollection>(theDTSegmentLabel);
   cscSegmentsToken = consumes<CSCSegmentCollection>(theCSCSegmentLabel);

   inputMuonShowerInformationValueMap = iConfig.getParameter<edm::InputTag>("inputMuonShowerInformationValueMap");
   inputMuonShowerInformationValueMapToken =  consumes<edm::ValueMap<reco::MuonShower>>(inputMuonShowerInformationValueMap);

}

//=======================================================================================================================================================================================================================//
RECOAnalysis::~RECOAnalysis()
{

}

//=======================================================================================================================================================================================================================//
void RECOAnalysis::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  
   /////////////////////////////////////////////////////////////////////////////////////
   ///////////////////////////////////// MAIN CODE /////////////////////////////////////
   /////////////////////////////////////////////////////////////////////////////////////

   //////////////////////////////// GET THE COLLECTIONS ////////////////////////////////

   edm::Handle<edm::View<reco::Muon> > muons;
   iEvent.getByToken(theMuonCollection, muons);

   edm::Handle<edm::View<reco::GenParticle> > genparticles;
   iEvent.getByToken(theGenParticleCollection, genparticles);

   edm::Handle<DTRecSegment4DCollection> dtSegments;
   iEvent.getByToken(dtSegmentsToken, dtSegments);

   edm::Handle<CSCSegmentCollection> cscSegments;
   iEvent.getByToken(cscSegmentsToken, cscSegments);

   edm::Handle<edm::ValueMap<reco::MuonShower>> muonShowerInformation;
   iEvent.getByToken(inputMuonShowerInformationValueMapToken, muonShowerInformation);

   theService->update(iSetup);

   /////////////////////////////////// EVENT INFO //////////////////////////////////////

   Event_id = iEvent.id().event();
   Event_run = iEvent.id().run();
   Event_luminosityBlock = iEvent.id().luminosityBlock();

   /////////////////////////////////////////////////////////////////////////////////////
   ///////////////////////////////// GET MUON VARIABLES ////////////////////////////////
   /////////////////////////////////////////////////////////////////////////////////////
   int iMuon = 0;
   unsigned int idx = 0;

   for (auto itmuon=muons->begin(); itmuon != muons->end(); itmuon++){

     if(!(itmuon->innerTrack().isNonnull())) continue;
     if(itmuon->innerTrack()->pt() < 200.) continue; // High-pT muons

     // MC Truth Matching + gen pT storage
     bool GenMatch = false;
     for (auto itgenparticle=genparticles->begin(); itgenparticle != genparticles->end(); itgenparticle++){

       if(!(std::abs(itgenparticle->pdgId()) == 13 && itgenparticle->status() == 1)) continue;
       
       if(deltaR(*itmuon, *itgenparticle) < 0.3){
     	 GenMatch = true;
     	 Muon_Genpt.push_back(itgenparticle->pt());
     	 break;
       }else{continue;}
     }
     
     if(!GenMatch) continue;

     iMuon++;     
     
     // Store Muon_* variables
     (itmuon->globalTrack().isNonnull())?Muon_GlbTrack_pt.push_back(itmuon->globalTrack()->pt()):Muon_GlbTrack_pt.push_back(-9999.);
     (itmuon->globalTrack().isNonnull())?Muon_GlbTrack_ptErr.push_back(itmuon->globalTrack()->ptError()):Muon_GlbTrack_ptErr.push_back(-9999.);
     (itmuon->globalTrack().isNonnull())?Muon_GlbTrack_eta.push_back(itmuon->globalTrack()->eta()):Muon_GlbTrack_eta.push_back(-9999.);
     (itmuon->globalTrack().isNonnull())?Muon_GlbTrack_phi.push_back(itmuon->globalTrack()->phi()):Muon_GlbTrack_phi.push_back(-9999.);
     (itmuon->globalTrack().isNonnull())?Muon_GlbTrack_charge.push_back(itmuon->globalTrack()->charge()):Muon_GlbTrack_charge.push_back(-9999.);
     (itmuon->globalTrack().isNonnull())?Muon_GlbTrack_Chindf.push_back(itmuon->globalTrack()->chi2()/(float)itmuon->globalTrack()->ndof()):Muon_GlbTrack_Chindf.push_back(-9999.);

     (itmuon->innerTrack().isNonnull())?Muon_InnerTrack_pt.push_back(itmuon->innerTrack()->pt()):Muon_InnerTrack_pt.push_back(-9999.);
     (itmuon->innerTrack().isNonnull())?Muon_InnerTrack_ptErr.push_back(itmuon->innerTrack()->ptError()):Muon_InnerTrack_ptErr.push_back(-9999.);
     (itmuon->innerTrack().isNonnull())?Muon_InnerTrack_eta.push_back(itmuon->innerTrack()->eta()):Muon_InnerTrack_eta.push_back(-9999.);
     (itmuon->innerTrack().isNonnull())?Muon_InnerTrack_phi.push_back(itmuon->innerTrack()->phi()):Muon_InnerTrack_phi.push_back(-9999.);
     (itmuon->innerTrack().isNonnull())?Muon_InnerTrack_charge.push_back(itmuon->innerTrack()->charge()):Muon_InnerTrack_charge.push_back(-9999.);
     (itmuon->innerTrack().isNonnull())?Muon_InnerTrack_Chindf.push_back(itmuon->innerTrack()->chi2()/(float)itmuon->innerTrack()->ndof()):Muon_InnerTrack_Chindf.push_back(-9999.);

     (itmuon->tunePMuonBestTrack().isNonnull())?Muon_TunePTrack_pt.push_back(itmuon->tunePMuonBestTrack()->pt()):Muon_TunePTrack_pt.push_back(-9999.);
     (itmuon->tunePMuonBestTrack().isNonnull())?Muon_TunePTrack_ptErr.push_back(itmuon->tunePMuonBestTrack()->ptError()):Muon_TunePTrack_ptErr.push_back(-9999.);
     (itmuon->tunePMuonBestTrack().isNonnull())?Muon_TunePTrack_eta.push_back(itmuon->tunePMuonBestTrack()->eta()):Muon_TunePTrack_eta.push_back(-9999.);
     (itmuon->tunePMuonBestTrack().isNonnull())?Muon_TunePTrack_phi.push_back(itmuon->tunePMuonBestTrack()->phi()):Muon_TunePTrack_phi.push_back(-9999.);
     (itmuon->tunePMuonBestTrack().isNonnull())?Muon_TunePTrack_charge.push_back(itmuon->tunePMuonBestTrack()->charge()):Muon_TunePTrack_charge.push_back(-9999.);
     (itmuon->tunePMuonBestTrack().isNonnull())?Muon_TunePTrack_Chindf.push_back(itmuon->tunePMuonBestTrack()->chi2()/(float)itmuon->tunePMuonBestTrack()->ndof()):Muon_TunePTrack_Chindf.push_back(-9999.);

     // Build the transient track of the tracker track and get its final state on surface
     reco::TransientTrack trackinner(itmuon->innerTrack(), &*theService->magneticField(), theService->trackingGeometry());     
     TrajectoryStateOnSurface outerTSOS = trackinner.outermostMeasurementState(); 

     std::map<const GeomDet*, std::vector<TrackingRecHit *> > DetAllSegmentsMap; 
     
     // Loop over dtSegments

     for (auto itHit = dtSegments->begin(); itHit != dtSegments->end(); itHit++) {
       
       //Only valid segments & hasZ & hasPhi
       if(!itHit->isValid() || !itHit->hasPhi()) continue;
       DetId myDet = itHit->geographicalId();
       const GeomDet *geomDet = theService->trackingGeometry()->idToDet(myDet);

       if(geomDet->geographicalId().subdetId()  == MuonSubdetId::DT){
	 DTWireId id(geomDet->geographicalId().rawId());
	 if(id.station() != 4 && !itHit->hasZed()){continue;}
	 if(id.station() != 4 && itHit->dimension()!=4){continue;}
       }
       


       //Get the GeomDet associated to this DetId
       std::map<const GeomDet*, std::vector<TrackingRecHit *> >::iterator it = DetAllSegmentsMap.find(geomDet);

       if(it == DetAllSegmentsMap.end()) {
	 //No -> we create a pair of GeomDet and vector of hits, and put the hit in the vector.
	 std::vector<TrackingRecHit *> trhit;
	 TrackingRecHit *rechitref = (TrackingRecHit *)&(*itHit);
	 trhit.push_back(rechitref);
	 DetAllSegmentsMap.insert(std::pair<const GeomDet*, std::vector<TrackingRecHit *> > (geomDet, trhit));
       } else {
	 //Yes -> we just put the hit in the corresponding hit vector.
	 TrackingRecHit *rechitref = (TrackingRecHit *) &(*itHit);
	 it->second.push_back(rechitref);
       }

     }

     // Loop over csc segments

     for (auto itHit = cscSegments->begin(); itHit != cscSegments->end(); itHit++) {

       //Only valid hits     
       if(!itHit->isValid()) continue;
       DetId myDet = itHit->geographicalId();
       const GeomDet *geomDet = theService->trackingGeometry()->idToDet(myDet);

       //Get the GeomDet associated to this DetId
       std::map<const GeomDet*, std::vector<TrackingRecHit *> >::iterator it = DetAllSegmentsMap.find(geomDet);

       if(it == DetAllSegmentsMap.end()) {
	 //No -> we create a pair of GeomDet and vector of hits, and put the hit in the vector.
	 std::vector<TrackingRecHit *> trhit;
	 TrackingRecHit *rechitref = (TrackingRecHit *)&(*itHit);
	 trhit.push_back(rechitref);
	 DetAllSegmentsMap.insert(std::pair<const GeomDet*, std::vector<TrackingRecHit *> > (geomDet, trhit));
       } else {
	 //Yes -> we just put the hit in the corresponding hit vector.
	 TrackingRecHit *rechitref = (TrackingRecHit *) &(*itHit);
	 it->second.push_back(rechitref);
       }

     }
     
     int iHit = 0;
     int iGeomDet = 0;
     //Now we do the extrapolations 

     for(auto it = DetAllSegmentsMap.begin(); it != DetAllSegmentsMap.end(); it++) {

         //Propagate
         std::pair<TrajectoryStateOnSurface, double> muonState = theService->propagator(propagator_)->propagateWithPath(outerTSOS, it->first->surface());
	 
	 if(muonState.second < 0.) continue; //backwards extrapolation

     	 // Store the hit if the extrapolation is valid
         if(muonState.first.isValid()){

	   GlobalPoint prop_gp = muonState.first.globalPosition();

	   //edm::RefToBase<reco::Muon> muRefTmp = muons->refAt(idx);
	   //reco::CandidateBaseRef muonBaseRef(muRefTmp);
	   //reco::MuonShower MuonShowerInfo = (*muonShowerInformation)[muonBaseRef];

	   if(it->first->geographicalId().subdetId()  == MuonSubdetId::DT){
	     DTWireId id(it->first->geographicalId().rawId());
	     if(id.station() == 4){
	       //if(dist2d_xy(prop_gp, it->first->surface().toGlobal(Local3DPoint(0.,0.,0.))) > 350) continue;
	       if(dist2d_xy(prop_gp, it->first->surface().toGlobal(Local3DPoint(0.,0.,0.))) > 260) continue;
	       Prop_isDT.push_back(1);
	       Prop_isCSC.push_back(0);
	       Prop_DTstation.push_back(id.station());
	       Prop_CSCstation.push_back(-9999);
	       Prop_DetElement.push_back(getdetid("DT", id.wheel(), id.station(), id.sector()));
	     }else if(id.station() == 3){ 
	       //if(dist2d_xz(prop_gp, it->first->surface().toGlobal(Local3DPoint(0.,0.,0.))) > 235. && MuonShowerInfo.stationShowerSizeT.at(id.station() - 1) < 235.) continue;
	       if(dist2d_xz(prop_gp, it->first->surface().toGlobal(Local3DPoint(0.,0.,0.))) > 235.) continue;
	       //if(dist2d_xz(prop_gp, it->first->surface().toGlobal(Local3DPoint(0.,0.,0.))) > 300.) continue;
	       Prop_isDT.push_back(1);
	       Prop_isCSC.push_back(0);
	       Prop_DTstation.push_back(id.station());
	       Prop_CSCstation.push_back(-9999);
	       Prop_DetElement.push_back(getdetid("DT", id.wheel(), id.station(), id.sector()));
	       //std::cout << "DetElement: " <<  getdetid("DT", id.wheel(), id.station(), id.sector()) << std::endl;
	       //std::cout << "Coords: (" << prop_gp.x() << "," << prop_gp.y() << "," << prop_gp.z() << ")" << std::endl; 
	     }else{
	       //if(dist2d_xz(prop_gp, it->first->surface().toGlobal(Local3DPoint(0.,0.,0.))) > 160. && MuonShowerInfo.stationShowerSizeT.at(id.station() - 1) < 160.) continue;
	       if(dist2d_xz(prop_gp, it->first->surface().toGlobal(Local3DPoint(0.,0.,0.))) > 160.) continue;
	       //if(dist2d_xz(prop_gp, it->first->surface().toGlobal(Local3DPoint(0.,0.,0.))) > 200.) continue;
	       Prop_isDT.push_back(1);
	       Prop_isCSC.push_back(0);
	       Prop_DTstation.push_back(id.station());
	       Prop_CSCstation.push_back(-9999);
	       Prop_DetElement.push_back(getdetid("DT", id.wheel(), id.station(), id.sector()));
	     }
	   }else if(it->first->geographicalId().subdetId()  == MuonSubdetId::CSC){
	     CSCDetId id(it->first->geographicalId().rawId());
	     if((id.station() == 2 || id.station() == 3 || id.station() ==4) && id.ring() == 2){ 
	       //if(dist2d_xy(prop_gp, it->first->surface().toGlobal(Local3DPoint(0.,0.,0.))) > 185. && MuonShowerInfo.stationShowerSizeT.at(id.station() - 1) < 185.) continue;
	       if(dist2d_xy(prop_gp, it->first->surface().toGlobal(Local3DPoint(0.,0.,0.))) > 185.) continue;
	       //if(dist2d_xy(prop_gp, it->first->surface().toGlobal(Local3DPoint(0.,0.,0.))) > 200.) continue;
	       Prop_isDT.push_back(0);
	       Prop_isCSC.push_back(1);
	       Prop_DTstation.push_back(-9999);
	       Prop_CSCstation.push_back(id.station());
	       Prop_DetElement.push_back(getdetid("CSC", id.endcap(), id.station(), id.ring()));
	     }else{
	       //if(dist2d_xy(prop_gp, it->first->surface().toGlobal(Local3DPoint(0.,0.,0.))) > 120. && MuonShowerInfo.stationShowerSizeT.at(id.station() - 1) < 120.) continue;
	       if(dist2d_xy(prop_gp, it->first->surface().toGlobal(Local3DPoint(0.,0.,0.))) > 120.) continue;
	       //if(dist2d_xy(prop_gp, it->first->surface().toGlobal(Local3DPoint(0.,0.,0.))) > 200.) continue;
	       Prop_isDT.push_back(0);
	       Prop_isCSC.push_back(1);
	       Prop_DTstation.push_back(-9999);
	       Prop_CSCstation.push_back(id.station());
	       Prop_DetElement.push_back(getdetid("CSC", id.endcap(), id.station(), id.ring()));
	     }
	   }

	   iGeomDet++;

     	   Prop_x.push_back(prop_gp.x());
     	   Prop_y.push_back(prop_gp.y());
     	   Prop_z.push_back(prop_gp.z());
	   Prop_Detid.push_back((*it).first->geographicalId().rawId());
	   Prop_Muonid.push_back(iMuon);
	   Prop_Eventid.push_back(iEvent.id().event());
	   Prop_EventluminosityBlock.push_back(iEvent.id().luminosityBlock());

	   for(int i=0; i<(int)(*it).second.size(); i++){
	     
	     LocalPoint hit_lp = (*it).second.at(i)->localPosition();
	     GlobalPoint hit_gp = it->first->surface().toGlobal(hit_lp);

	     float distToProp = dist3d(hit_gp,prop_gp);

	     iHit++;

	     if((*it).second.at(i)->geographicalId().subdetId() == MuonSubdetId::DT){
	       DTWireId id((*it).second.at(i)->geographicalId().rawId());
	       Hit_DetElement.push_back(getdetid("DT", id.wheel(), id.station(), id.sector()));
	       Hit_isDT.push_back(1);
	       Hit_isCSC.push_back(0);
	       Hit_DTstation.push_back(id.station());
	       Hit_CSCstation.push_back(-9999);
	       if(id.station()==4){distToProp = dist2d_xy(hit_gp,prop_gp);} // set xy distance for DT station4 (no Z coord)
	       DTRecSegment4D *mySegment = dynamic_cast<DTRecSegment4D *>((*it).second.at(i));
	       StateSegmentMatcher SegmentComp(outerTSOS, *mySegment, mySegment->localDirectionError());
	       Hit_Compatibility.push_back(SegmentComp.value());
	       GlobalVector gv = it->first->surface().toGlobal(mySegment->localDirection());
	       Hit_dirx.push_back(gv.x());
	       Hit_diry.push_back(gv.y());
	       Hit_dirz.push_back(gv.z());
	       Hit_chi2.push_back(mySegment->chi2());
	       Hit_ndof.push_back(mySegment->degreesOfFreedom());
	       // if(distToProp > 110 && distToProp < 130 && id.station()!=4){
	       // 	 //if(distToProp < 20 && id.station()!=4){
	       // 	 std::cout << "#################################################" << std::endl;
	       // 	 std::cout << "Event: " << iEvent.id().event() << std::endl; 
	       // 	 std::cout << "Muon: " << iMuon << std::endl;
	       // 	 std::cout << "Segment has phi: " << mySegment->hasPhi() << std::endl;
	       // 	 std::cout << "Segment has Z: " << mySegment->hasZed() << std::endl;
	       // 	 //std::cout << "GeomDet position: " << it->first->surface().position() << std::endl;
	       // }

	     }else if((*it).second.at(i)->geographicalId().subdetId() == MuonSubdetId::CSC){
	       CSCDetId id((*it).second.at(i)->geographicalId().rawId());
	       Hit_DetElement.push_back(getdetid("CSC", id.endcap(), id.station(), id.ring()));
	       Hit_isDT.push_back(0);
	       Hit_isCSC.push_back(1);
	       Hit_DTstation.push_back(-9999);
	       Hit_CSCstation.push_back(id.station());
	       CSCSegment *mySegment = dynamic_cast<CSCSegment *>((*it).second.at(i));
	       StateSegmentMatcher SegmentComp(outerTSOS, *mySegment, mySegment->localDirectionError());
	       Hit_Compatibility.push_back(SegmentComp.value());
	       GlobalVector gv = it->first->surface().toGlobal(mySegment->localDirection());
	       Hit_dirx.push_back(gv.x());
	       Hit_diry.push_back(gv.y());
	       Hit_dirz.push_back(gv.z());
	       Hit_chi2.push_back(mySegment->chi2());
	       Hit_ndof.push_back(mySegment->degreesOfFreedom());
	     }

	     Hit_x.push_back(hit_gp.x()); 
	     Hit_y.push_back(hit_gp.y()); 
	     Hit_z.push_back(hit_gp.z()); 
	     Hit_Detid.push_back((*it).first->geographicalId().rawId()); 
	     Hit_Hitid.push_back(iHit);
	     Hit_distToProp.push_back(distToProp);
	     Hit_Muonid.push_back(iMuon);
	     Hit_Eventid.push_back(iEvent.id().event());
	     Hit_EventluminosityBlock.push_back(iEvent.id().luminosityBlock());
	   }
	   
     	 }

     }    
     
     Muon_nGeomDets.push_back(iGeomDet);
     Muon_nHits.push_back(iHit);
     Muon_nShowers.push_back(itmuon->numberOfShowers());
     Muon_hasShowerInStation_DT_1.push_back(itmuon->hasShowerInStation(1,MuonSubdetId::DT));
     Muon_hasShowerInStation_DT_2.push_back(itmuon->hasShowerInStation(2,MuonSubdetId::DT));
     Muon_hasShowerInStation_DT_3.push_back(itmuon->hasShowerInStation(3,MuonSubdetId::DT));
     Muon_hasShowerInStation_DT_4.push_back(itmuon->hasShowerInStation(4,MuonSubdetId::DT));
     Muon_hasShowerInStation_CSC_1.push_back(itmuon->hasShowerInStation(1,MuonSubdetId::CSC));
     Muon_hasShowerInStation_CSC_2.push_back(itmuon->hasShowerInStation(2,MuonSubdetId::CSC));
     Muon_hasShowerInStation_CSC_3.push_back(itmuon->hasShowerInStation(3,MuonSubdetId::CSC));
     Muon_hasShowerInStation_CSC_4.push_back(itmuon->hasShowerInStation(4,MuonSubdetId::CSC));
     Muon_nDigisInStation_DT_1.push_back(itmuon->nDigisInStation(1,MuonSubdetId::DT));
     Muon_nDigisInStation_DT_2.push_back(itmuon->nDigisInStation(2,MuonSubdetId::DT));
     Muon_nDigisInStation_DT_3.push_back(itmuon->nDigisInStation(3,MuonSubdetId::DT));
     Muon_nDigisInStation_DT_4.push_back(itmuon->nDigisInStation(4,MuonSubdetId::DT));
     Muon_nDigisInStation_CSC_1.push_back(itmuon->nDigisInStation(1,MuonSubdetId::CSC));
     Muon_nDigisInStation_CSC_2.push_back(itmuon->nDigisInStation(2,MuonSubdetId::CSC));
     Muon_nDigisInStation_CSC_3.push_back(itmuon->nDigisInStation(3,MuonSubdetId::CSC));
     Muon_nDigisInStation_CSC_4.push_back(itmuon->nDigisInStation(4,MuonSubdetId::CSC));
     Muon_Muonid.push_back(iMuon);
     Muon_Eventid.push_back(iEvent.id().event());
     Muon_EventluminosityBlock.push_back(iEvent.id().luminosityBlock());

     idx++;
   }

   Event_nMuons = iMuon;
     
   /////////////////////////////////////////////////////////////////////////////////////
   /////////////////////////////////// FILL THE TREE ///////////////////////////////////
   /////////////////////////////////////////////////////////////////////////////////////

   tree_out->Fill();

   Muon_Muonid.clear();
   Muon_Eventid.clear();
   Muon_EventluminosityBlock.clear();
   Muon_nGeomDets.clear();
   Muon_nHits.clear();
   Muon_hasShowerInStation_DT_1.clear();
   Muon_hasShowerInStation_DT_2.clear();
   Muon_hasShowerInStation_DT_3.clear();
   Muon_hasShowerInStation_DT_4.clear();
   Muon_hasShowerInStation_CSC_1.clear();
   Muon_hasShowerInStation_CSC_2.clear();
   Muon_hasShowerInStation_CSC_3.clear();
   Muon_hasShowerInStation_CSC_4.clear();
   Muon_nDigisInStation_DT_1.clear();
   Muon_nDigisInStation_DT_2.clear();
   Muon_nDigisInStation_DT_3.clear();
   Muon_nDigisInStation_DT_4.clear();
   Muon_nDigisInStation_CSC_1.clear();
   Muon_nDigisInStation_CSC_2.clear();
   Muon_nDigisInStation_CSC_3.clear();
   Muon_nDigisInStation_CSC_4.clear();
   Muon_nShowers.clear();
   Muon_Genpt.clear();
   Muon_GlbTrack_pt.clear();
   Muon_GlbTrack_eta.clear();
   Muon_GlbTrack_phi.clear();
   Muon_GlbTrack_charge.clear();
   Muon_GlbTrack_ptErr.clear();
   Muon_GlbTrack_Chindf.clear();
   Muon_InnerTrack_pt.clear();
   Muon_InnerTrack_eta.clear();
   Muon_InnerTrack_phi.clear();
   Muon_InnerTrack_charge.clear();
   Muon_InnerTrack_ptErr.clear();
   Muon_InnerTrack_Chindf.clear();
   Muon_TunePTrack_pt.clear();
   Muon_TunePTrack_eta.clear();
   Muon_TunePTrack_phi.clear();
   Muon_TunePTrack_charge.clear();
   Muon_TunePTrack_ptErr.clear();
   Muon_TunePTrack_Chindf.clear();
   Prop_x.clear();
   Prop_y.clear();
   Prop_z.clear();
   Prop_Detid.clear();
   Prop_Muonid.clear();
   Prop_Eventid.clear();
   Prop_EventluminosityBlock.clear();
   Prop_isDT.clear();
   Prop_isCSC.clear();
   Prop_DTstation.clear();
   Prop_CSCstation.clear();
   Prop_DetElement.clear();

   Hit_x.clear();
   Hit_y.clear();
   Hit_z.clear();
   Hit_distToProp.clear();
   Hit_Detid.clear();
   Hit_Hitid.clear();
   Hit_Muonid.clear();
   Hit_Eventid.clear();
   Hit_EventluminosityBlock.clear();
   Hit_DetElement.clear();
   Hit_isDT.clear();
   Hit_isCSC.clear();
   Hit_DTstation.clear();
   Hit_CSCstation.clear();
   Hit_Compatibility.clear();
   Hit_dirx.clear();
   Hit_diry.clear();
   Hit_dirz.clear();
   Hit_chi2.clear();
   Hit_ndof.clear();

}

//=======================================================================================================================================================================================================================//
void RECOAnalysis::beginJob()
{

    // Output file definition
    output_filename = parameters.getParameter<std::string>("nameOfOutput");
    file_out = new TFile(output_filename.c_str(), "RECREATE");
    file_out->cd();

    std::cout << "the file is created" << std::endl;
    
    // Output Tree definition
    tree_out = new TTree("Events", "Events");

    // ///////////////////////////////// EVENT INFO BRANCHES ///////////////////////////////

    tree_out->Branch("Event_id", &Event_id, "Event_id/I");
    tree_out->Branch("Event_run", &Event_run, "Event_run/I");
    tree_out->Branch("Event_luminosityBlock", &Event_luminosityBlock, "Event_luminosityBlock/I");
    tree_out->Branch("Event_nMuons", &Event_nMuons, "Event_nMuons/I");

    // ////////////////////////////// MUON BRANCHES //////////////////////////////

    tree_out->Branch("Muon_Eventid", "vector<int>", &Muon_Eventid);
    tree_out->Branch("Muon_EventluminosityBlock", "vector<int>", &Muon_EventluminosityBlock);
    tree_out->Branch("Muon_Muonid",  "vector<int>", &Muon_Muonid);
    tree_out->Branch("Muon_Genpt", "vector<float>", &Muon_Genpt);
    tree_out->Branch("Muon_GlbTrack_pt", "vector<float>", &Muon_GlbTrack_pt);
    tree_out->Branch("Muon_GlbTrack_eta", "vector<float>", &Muon_GlbTrack_eta);
    tree_out->Branch("Muon_GlbTrack_phi", "vector<float>", &Muon_GlbTrack_phi);
    tree_out->Branch("Muon_GlbTrack_charge", "vector<float>", &Muon_GlbTrack_charge);
    tree_out->Branch("Muon_GlbTrack_ptErr", "vector<float>", &Muon_GlbTrack_ptErr);
    tree_out->Branch("Muon_GlbTrack_Chindf", "vector<float>", &Muon_GlbTrack_Chindf);
    tree_out->Branch("Muon_InnerTrack_pt", "vector<float>", &Muon_InnerTrack_pt);
    tree_out->Branch("Muon_InnerTrack_eta", "vector<float>", &Muon_InnerTrack_eta);
    tree_out->Branch("Muon_InnerTrack_phi", "vector<float>", &Muon_InnerTrack_phi);
    tree_out->Branch("Muon_InnerTrack_charge", "vector<float>", &Muon_InnerTrack_charge);
    tree_out->Branch("Muon_InnerTrack_ptErr", "vector<float>", &Muon_InnerTrack_ptErr);
    tree_out->Branch("Muon_InnerTrack_Chindf", "vector<float>", &Muon_InnerTrack_Chindf);
    tree_out->Branch("Muon_TunePTrack_pt", "vector<float>", &Muon_TunePTrack_pt);
    tree_out->Branch("Muon_TunePTrack_eta", "vector<float>", &Muon_TunePTrack_eta);
    tree_out->Branch("Muon_TunePTrack_phi", "vector<float>", &Muon_TunePTrack_phi);
    tree_out->Branch("Muon_TunePTrack_charge", "vector<float>", &Muon_TunePTrack_charge);
    tree_out->Branch("Muon_TunePTrack_ptErr", "vector<float>", &Muon_TunePTrack_ptErr);
    tree_out->Branch("Muon_TunePTrack_Chindf", "vector<float>", &Muon_TunePTrack_Chindf);
    tree_out->Branch("Muon_nGeomDets", "vector<int>", &Muon_nGeomDets);
    tree_out->Branch("Muon_nHits", "vector<int>", &Muon_nHits);
    tree_out->Branch("Muon_nShowers", "vector<int>", &Muon_nShowers);
    tree_out->Branch("Muon_hasShowerInStation_DT_1", "vector<int>", &Muon_hasShowerInStation_DT_1);
    tree_out->Branch("Muon_hasShowerInStation_DT_2", "vector<int>", &Muon_hasShowerInStation_DT_2);
    tree_out->Branch("Muon_hasShowerInStation_DT_3", "vector<int>", &Muon_hasShowerInStation_DT_3);
    tree_out->Branch("Muon_hasShowerInStation_DT_4", "vector<int>", &Muon_hasShowerInStation_DT_4);
    tree_out->Branch("Muon_hasShowerInStation_CSC_1", "vector<int>", &Muon_hasShowerInStation_CSC_1);
    tree_out->Branch("Muon_hasShowerInStation_CSC_2", "vector<int>", &Muon_hasShowerInStation_CSC_2);
    tree_out->Branch("Muon_hasShowerInStation_CSC_3", "vector<int>", &Muon_hasShowerInStation_CSC_3);
    tree_out->Branch("Muon_hasShowerInStation_CSC_4", "vector<int>", &Muon_hasShowerInStation_CSC_4);
    tree_out->Branch("Muon_nDigisInStation_DT_1", "vector<int>", &Muon_nDigisInStation_DT_1);
    tree_out->Branch("Muon_nDigisInStation_DT_2", "vector<int>", &Muon_nDigisInStation_DT_2);
    tree_out->Branch("Muon_nDigisInStation_DT_3", "vector<int>", &Muon_nDigisInStation_DT_3);
    tree_out->Branch("Muon_nDigisInStation_DT_4", "vector<int>", &Muon_nDigisInStation_DT_4);
    tree_out->Branch("Muon_nDigisInStation_CSC_1", "vector<int>", &Muon_nDigisInStation_CSC_1);
    tree_out->Branch("Muon_nDigisInStation_CSC_2", "vector<int>", &Muon_nDigisInStation_CSC_2);
    tree_out->Branch("Muon_nDigisInStation_CSC_3", "vector<int>", &Muon_nDigisInStation_CSC_3);
    tree_out->Branch("Muon_nDigisInStation_CSC_4", "vector<int>", &Muon_nDigisInStation_CSC_4);

    // ////////////////////////////// HIT & EXTRAPOLATION BRANCHES //////////////////////////////

    tree_out->Branch("Prop_x", "vector<float>", &Prop_x);
    tree_out->Branch("Prop_y", "vector<float>", &Prop_y);
    tree_out->Branch("Prop_z", "vector<float>", &Prop_z);
    tree_out->Branch("Prop_Detid", "vector<unsigned int>", &Prop_Detid);
    tree_out->Branch("Prop_Muonid", "vector<int>", &Prop_Muonid);
    tree_out->Branch("Prop_Eventid", "vector<int>", &Prop_Eventid);
    tree_out->Branch("Prop_EventluminosityBlock", "vector<int>", &Prop_EventluminosityBlock);
    tree_out->Branch("Prop_DetElement", "vector<string>", &Prop_DetElement);
    tree_out->Branch("Prop_isDT", "vector<int>", &Prop_isDT);
    tree_out->Branch("Prop_isCSC", "vector<int>", &Prop_isCSC);
    tree_out->Branch("Prop_DTstation", "vector<int>", &Prop_DTstation);
    tree_out->Branch("Prop_CSCstation", "vector<int>", &Prop_CSCstation);

    tree_out->Branch("Hit_x", "vector<float>", &Hit_x);
    tree_out->Branch("Hit_y", "vector<float>", &Hit_y);
    tree_out->Branch("Hit_z", "vector<float>", &Hit_z);
    tree_out->Branch("Hit_distToProp", "vector<float>", &Hit_distToProp);
    tree_out->Branch("Hit_Detid", "vector<unsigned int>", &Hit_Detid);
    tree_out->Branch("Hit_Hitid", "vector<int>", &Hit_Hitid);
    tree_out->Branch("Hit_Muonid", "vector<int>", &Hit_Muonid);
    tree_out->Branch("Hit_Eventid", "vector<int>", &Hit_Eventid);
    tree_out->Branch("Hit_EventluminosityBlock", "vector<int>", &Hit_EventluminosityBlock);
    tree_out->Branch("Hit_DetElement", "vector<string>", &Hit_DetElement);
    tree_out->Branch("Hit_isDT", "vector<int>", &Hit_isDT);
    tree_out->Branch("Hit_isCSC", "vector<int>", &Hit_isCSC);
    tree_out->Branch("Hit_DTstation", "vector<int>", &Hit_DTstation);
    tree_out->Branch("Hit_CSCstation", "vector<int>", &Hit_CSCstation);
    tree_out->Branch("Hit_Compatibility", "vector<float>", &Hit_Compatibility);
    tree_out->Branch("Hit_dirx", "vector<float>", &Hit_dirx);
    tree_out->Branch("Hit_diry", "vector<float>", &Hit_diry);
    tree_out->Branch("Hit_dirz", "vector<float>", &Hit_dirz);
    tree_out->Branch("Hit_chi2", "vector<float>", &Hit_chi2);
    tree_out->Branch("Hit_ndof", "vector<int>", &Hit_ndof);


}

//=======================================================================================================================================================================================================================//
void RECOAnalysis::endJob() 
{

    std::cout << "The event is writen" << std::endl;
    file_out->cd();
    tree_out->Write();
    file_out->Close();

}

//=======================================================================================================================================================================================================================//
void RECOAnalysis::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

DEFINE_FWK_MODULE(RECOAnalysis);
