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
#include "DataFormats/MuonDetId/interface/MuonSubdetId.h"
#include "DataFormats/MuonDetId/interface/DTWireId.h"
#include "DataFormats/MuonDetId/interface/CSCDetId.h"
#include "DataFormats/MuonDetId/interface/RPCDetId.h"
#include "DataFormats/SiPixelDetId/interface/PXBDetId.h"
#include "DataFormats/SiPixelDetId/interface/PXFDetId.h"
#include "DataFormats/SiStripDetId/interface/SiStripDetId.h"

#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/GenRunInfoProduct.h"

#include "DataFormats/Common/interface/Handle.h"
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

// muon info
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


/////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////// DATA DEFINITION //////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////


////////////////////////////////////// BRANCHES /////////////////////////////////////

//-> EVENT INFO
Int_t Event_event;
Int_t Event_luminosityBlock;
Int_t Event_run;


//-> MUON INFO
Int_t   nMuons;

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
// std::vector<float> Muon_Pickypt;
// std::vector<float> Muon_PickyptErr;
// std::vector<float> Muon_PickyChindf;
// std::vector<float> Muon_Dytpt;
// std::vector<float> Muon_DytptErr;
// std::vector<float> Muon_DytChindf;
std::vector<float> Muon_TunePTrack_pt;
std::vector<float> Muon_TunePTrack_eta;
std::vector<float> Muon_TunePTrack_phi;
std::vector<float> Muon_TunePTrack_charge;
std::vector<float> Muon_TunePTrack_ptErr;
std::vector<float> Muon_TunePTrack_Chindf;

std::vector<std::vector<float> > Hit_Track_x;
std::vector<std::vector<float> > Hit_Track_y;
std::vector<std::vector<float> > Hit_Track_z;
std::vector<std::vector<int> > Hit_Track_subdetid;
std::vector<std::vector<int> > Hit_Track_DT_station;
std::vector<std::vector<int> > Hit_Track_DT_layer;
std::vector<std::vector<int> > Hit_Track_DT_superlayer;
std::vector<std::vector<int> > Hit_Track_DT_wheel;
std::vector<std::vector<int> > Hit_Track_DT_sector;
std::vector<std::vector<int> > Hit_Track_CSC_endcap;
std::vector<std::vector<int> > Hit_Track_CSC_station;
std::vector<std::vector<int> > Hit_Track_CSC_ringN;
std::vector<std::vector<int> > Hit_Track_CSC_chamber;
std::vector<std::vector<int> > Hit_Track_CSC_layer;
std::vector<std::vector<float> > Hit_prop_x;
std::vector<std::vector<float> > Hit_prop_y;
std::vector<std::vector<float> > Hit_prop_z;

std::vector<std::vector<float> > Hit_DetAll_x;
std::vector<std::vector<float> > Hit_DetAll_y;
std::vector<std::vector<float> > Hit_DetAll_z;
std::vector<std::vector<int> > Hit_DetAll_subdetid;
std::vector<std::vector<int> > Hit_DetAll_DT_station;
std::vector<std::vector<int> > Hit_DetAll_DT_layer;
std::vector<std::vector<int> > Hit_DetAll_DT_superlayer;
std::vector<std::vector<int> > Hit_DetAll_DT_wheel;
std::vector<std::vector<int> > Hit_DetAll_DT_sector;
std::vector<std::vector<int> > Hit_DetAll_CSC_endcap;
std::vector<std::vector<int> > Hit_DetAll_CSC_station;
std::vector<std::vector<int> > Hit_DetAll_CSC_ringN;
std::vector<std::vector<int> > Hit_DetAll_CSC_chamber;
std::vector<std::vector<int> > Hit_DetAll_CSC_layer;

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

      edm::InputTag theDTSegmentLabel, theCSCSegmentLabel;
      edm::EDGetTokenT<DTRecSegment4DCollection> dtSegmentsToken;
      edm::EDGetTokenT<CSCSegmentCollection> cscSegmentsToken;

      edm::EDGetTokenT<edm::View<reco::Muon> > theMuonCollection;
      
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
   propagator_ = iConfig.getParameter<std::string>("Propagator");
   edm::ParameterSet serviceParameters = iConfig.getParameter<edm::ParameterSet>("ServiceParameters");
   theService = new MuonServiceProxy(serviceParameters);

   theDTSegmentLabel = iConfig.getParameter<edm::InputTag>("segmentsDt");
   theCSCSegmentLabel = iConfig.getParameter<edm::InputTag>("segmentsCSC");
   dtSegmentsToken = consumes<DTRecSegment4DCollection>(theDTSegmentLabel);
   cscSegmentsToken = consumes<CSCSegmentCollection>(theCSCSegmentLabel);

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

   edm::Handle<DTRecSegment4DCollection> dtSegments;
   iEvent.getByToken(dtSegmentsToken, dtSegments);

   edm::Handle<CSCSegmentCollection> cscSegments;
   iEvent.getByToken(cscSegmentsToken, cscSegments);

   //-------------------------------------------------------------------------//
   //------------------------Added by Pablo-----------------------------------//
   //-------------------------------------------------------------------------//
   theService->update(iSetup);
   //-------------------------------------------------------------------------//
   //-------------------------------------------------------------------------//

   /////////////////////////////////// EVENT INFO //////////////////////////////////////

   Event_event = iEvent.id().event();
   Event_run = iEvent.id().run();
   Event_luminosityBlock = iEvent.id().luminosityBlock();

 
   /////////////////////////////////////////////////////////////////////////////////////
   ///////////////////////////////// GET MUON VARIABLES ////////////////////////////////
   /////////////////////////////////////////////////////////////////////////////////////
   int iMuon = 0;

   for (auto itmuon=muons->begin(); itmuon != muons->end(); itmuon++){

     if(!(itmuon->innerTrack().isNonnull() && itmuon->globalTrack().isNonnull())){continue;}
     if(itmuon->globalTrack()->pt() < 200.){continue;} // High-pT muons

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

     std::vector<float> temp_x;
     std::vector<float> temp_y;
     std::vector<float> temp_z;
     std::vector<float> temp_prop_x;
     std::vector<float> temp_prop_y;
     std::vector<float> temp_prop_z;
     std::vector<int> temp_subdetid;
     std::vector<int> temp_DT_station;
     std::vector<int> temp_DT_layer;
     std::vector<int> temp_DT_superlayer;
     std::vector<int> temp_DT_wheel;
     std::vector<int> temp_DT_sector;
     std::vector<int> temp_CSC_endcap;
     std::vector<int> temp_CSC_station;
     std::vector<int> temp_CSC_ringN;
     std::vector<int> temp_CSC_chamber;
     std::vector<int> temp_CSC_layer;

     std::vector<float> temp_all_x;
     std::vector<float> temp_all_y;
     std::vector<float> temp_all_z;
     std::vector<float> temp_all_prop_x;
     std::vector<float> temp_all_prop_y;
     std::vector<float> temp_all_prop_z;
     std::vector<int> temp_all_subdetid;
     std::vector<int> temp_all_DT_station;
     std::vector<int> temp_all_DT_layer;
     std::vector<int> temp_all_DT_superlayer;
     std::vector<int> temp_all_DT_wheel;
     std::vector<int> temp_all_DT_sector;
     std::vector<int> temp_all_CSC_endcap;
     std::vector<int> temp_all_CSC_station;
     std::vector<int> temp_all_CSC_ringN;
     std::vector<int> temp_all_CSC_chamber;
     std::vector<int> temp_all_CSC_layer;



     //-------------------------------------------------------------------------//
     //------------------------Added by Pablo-----------------------------------//
     //-------------------------------------------------------------------------//

     //We build the transient track of the tracker track and get its final state on surface
     reco::TransientTrack trackinner(itmuon->innerTrack(), &*theService->magneticField(), theService->trackingGeometry());     
     TrajectoryStateOnSurface outerTSOS = trackinner.outermostMeasurementState(); 

     //We also build the transient track of the global track
     reco::TransientTrack track(itmuon->globalTrack(), &*theService->magneticField(), theService->trackingGeometry());     

     //We make a map in which to store a vector of trackingrechits per detector in the muon system
     std::map<const GeomDet*, std::vector<TrackingRecHit *> > DetRecHitMap; 

     for (auto itHit = track.recHitsBegin(); itHit != track.recHitsEnd(); itHit++) {
       //Hit_muonLink.push_back(iMuon); //meter en el std::map el muonLink como third element
        //Only valid hits     
        if(!(*itHit)->isValid()) continue;
        DetId myDet = (*itHit)->geographicalId();
        //Only if the hit is in the muon system
        if(!(myDet.det() == DetId::Muon)) continue;
        //Only if it's a DT or CSC
        if(!(myDet.subdetId() == MuonSubdetId::DT || myDet.subdetId() == MuonSubdetId::CSC)) continue;
        //Get the GeomDet associated to this DetIt 
        const GeomDet *geomDet = theService->trackingGeometry()->idToDet(myDet);
        //Is it already in the map?
        std::map<const GeomDet*, std::vector<TrackingRecHit *> >::iterator it = DetRecHitMap.find(geomDet);
        if(it == DetRecHitMap.end()) {
            //No -> we create a pair of GeomDet and vector of hits, and put the hit in the vector.
            std::vector<TrackingRecHit *> trhit;
            trhit.push_back(*itHit);
            DetRecHitMap.insert(std::pair<const GeomDet*, std::vector<TrackingRecHit *> > (geomDet, trhit));

	    temp_x.push_back((*itHit)->localPosition().x()); 
	    temp_y.push_back((*itHit)->localPosition().y()); 
            temp_z.push_back((*itHit)->localPosition().z()); 
            temp_subdetid.push_back(geomDet->geographicalId().subdetId());
	    if((*itHit)->geographicalId().subdetId() == MuonSubdetId::CSC){
	      CSCDetId id((*itHit)->geographicalId().rawId());
	      temp_CSC_endcap.push_back(id.endcap());
	      temp_CSC_station.push_back(id.station());
	      temp_CSC_ringN.push_back(id.ring());
	      temp_CSC_chamber.push_back(id.chamber());
	      temp_CSC_layer.push_back(id.layer());
	      temp_DT_station.push_back(-9999);
	      temp_DT_layer.push_back(-9999);
	      temp_DT_superlayer.push_back(-9999);
	      temp_DT_wheel.push_back(-9999);
	      temp_DT_sector.push_back(-9999);
	      
	    }else if((*itHit)->geographicalId().subdetId() == MuonSubdetId::DT){
	      DTWireId id((*itHit)->geographicalId().rawId());
	      temp_CSC_endcap.push_back(-9999);
	      temp_CSC_station.push_back(-9999);
	      temp_CSC_ringN.push_back(-9999);
	      temp_CSC_chamber.push_back(-9999);
	      temp_CSC_layer.push_back(-9999);
	      temp_DT_station.push_back(id.station());
	      temp_DT_layer.push_back(id.layer());
	      temp_DT_superlayer.push_back(id.superLayer());
	      temp_DT_wheel.push_back(id.wheel());
	      temp_DT_sector.push_back(id.sector());
	    }

        } else { 
            //Yes -> we just put the hit in the corresponding hit vector.
	  it->second.push_back(*itHit);

	  temp_x.push_back((*itHit)->localPosition().x()); 
	  temp_y.push_back((*itHit)->localPosition().y()); 
	  temp_z.push_back((*itHit)->localPosition().z()); 
	  temp_subdetid.push_back(geomDet->geographicalId().subdetId());
	  if((*itHit)->geographicalId().subdetId() == MuonSubdetId::CSC){
	    CSCDetId id((*itHit)->geographicalId().rawId());
	    temp_CSC_endcap.push_back(id.endcap());
	    temp_CSC_station.push_back(id.station());
	    temp_CSC_ringN.push_back(id.ring());
	    temp_CSC_chamber.push_back(id.chamber());
	    temp_CSC_layer.push_back(id.layer());
	    temp_DT_station.push_back(-9999);
	    temp_DT_layer.push_back(-9999);
	    temp_DT_superlayer.push_back(-9999);
	    temp_DT_wheel.push_back(-9999);
	    temp_DT_sector.push_back(-9999);
	    
	  }else if((*itHit)->geographicalId().subdetId() == MuonSubdetId::DT){
	    DTWireId id((*itHit)->geographicalId().rawId());
	    temp_CSC_endcap.push_back(-9999);
	    temp_CSC_station.push_back(-9999);
	    temp_CSC_ringN.push_back(-9999);
	    temp_CSC_chamber.push_back(-9999);
	    temp_CSC_layer.push_back(-9999);
	    temp_DT_station.push_back(id.station());
	    temp_DT_layer.push_back(id.layer());
	    temp_DT_superlayer.push_back(id.superLayer());
	    temp_DT_wheel.push_back(id.wheel());
	    temp_DT_sector.push_back(id.sector());
	  }

	}
     }


     Hit_Track_x.push_back(temp_x);
     Hit_Track_y.push_back(temp_y);
     Hit_Track_z.push_back(temp_z);
     Hit_Track_subdetid.push_back(temp_subdetid);
     Hit_Track_DT_station.push_back(temp_DT_station);
     Hit_Track_DT_layer.push_back(temp_DT_layer);
     Hit_Track_DT_superlayer.push_back(temp_DT_superlayer);
     Hit_Track_DT_wheel.push_back(temp_DT_wheel);
     Hit_Track_DT_sector.push_back(temp_DT_sector);
     Hit_Track_CSC_endcap.push_back(temp_CSC_endcap);
     Hit_Track_CSC_station.push_back(temp_CSC_station);
     Hit_Track_CSC_ringN.push_back(temp_CSC_ringN);
     Hit_Track_CSC_chamber.push_back(temp_CSC_chamber);
     Hit_Track_CSC_layer.push_back(temp_CSC_layer);


     std::map<const GeomDet*, std::vector<const TrackingRecHit *> > DetAllSegmentsMap; 

     for(auto ittrack = DetRecHitMap.begin(); ittrack != DetRecHitMap.end(); ittrack++) {

        DetId originalDet = ittrack->first->geographicalId();
        for (auto itHit = dtSegments->begin(); itHit != dtSegments->end(); itHit++) {
           //Only valid hits     
           if(!itHit->isValid()) continue;
           DetId myDet = itHit->geographicalId();
           if(myDet != originalDet) continue;
           //Get the GeomDet associated to this DetIt 
           std::map<const GeomDet*, std::vector<const TrackingRecHit *> >::iterator it = DetAllSegmentsMap.find(ittrack->first);
           if(it == DetAllSegmentsMap.end()) {
               //No -> we create a pair of GeomDet and vector of hits, and put the hit in the vector.
               std::vector<const TrackingRecHit *> trhit;
               const TrackingRecHit *rechitref = (const TrackingRecHit *)&(*itHit);
               trhit.push_back(rechitref);
               DetAllSegmentsMap.insert(std::pair<const GeomDet*, std::vector<const TrackingRecHit *> > (ittrack->first, trhit));
           } else {
               //Yes -> we just put the hit in the corresponding hit vector.
               const TrackingRecHit *rechitref = (const TrackingRecHit *) &(*itHit);
               it->second.push_back(rechitref);
           }
        }
        for (auto itHit = cscSegments->begin(); itHit != cscSegments->end(); itHit++) {
           //Only valid hits     
           if(!itHit->isValid()) continue;
           DetId myDet = itHit->geographicalId();
           if(myDet != originalDet) continue;
           //Get the GeomDet associated to this DetIt 
           std::map<const GeomDet*, std::vector<const TrackingRecHit *> >::iterator it = DetAllSegmentsMap.find(ittrack->first);
           if(it == DetAllSegmentsMap.end()) {
               //No -> we create a pair of GeomDet and vector of hits, and put the hit in the vector.
               std::vector<const TrackingRecHit *> trhit;
               const TrackingRecHit *rechitref = (const TrackingRecHit *)&(*itHit);
               trhit.push_back(rechitref);
               DetAllSegmentsMap.insert(std::pair <const GeomDet*, std::vector<const TrackingRecHit *> > (ittrack->first, trhit));

	       temp_all_x.push_back(rechitref->localPosition().x()); 
	       temp_all_y.push_back(rechitref->localPosition().y()); 
	       temp_all_z.push_back(rechitref->localPosition().z()); 
	       temp_all_subdetid.push_back(ittrack->first->geographicalId().subdetId());
	       if(rechitref->geographicalId().subdetId() == MuonSubdetId::CSC){
		 CSCDetId id(rechitref->geographicalId().rawId());
		 temp_all_CSC_endcap.push_back(id.endcap());
		 temp_all_CSC_station.push_back(id.station());
		 temp_all_CSC_ringN.push_back(id.ring());
		 temp_all_CSC_chamber.push_back(id.chamber());
		 temp_all_CSC_layer.push_back(id.layer());
		 temp_all_DT_station.push_back(-9999);
		 temp_all_DT_layer.push_back(-9999);
		 temp_all_DT_superlayer.push_back(-9999);
		 temp_all_DT_wheel.push_back(-9999);
		 temp_all_DT_sector.push_back(-9999);
		 
	       }else if(rechitref->geographicalId().subdetId() == MuonSubdetId::DT){
		 DTWireId id(rechitref->geographicalId().rawId());
		 temp_all_CSC_endcap.push_back(-9999);
		 temp_all_CSC_station.push_back(-9999);
		 temp_all_CSC_ringN.push_back(-9999);
		 temp_all_CSC_chamber.push_back(-9999);
		 temp_all_CSC_layer.push_back(-9999);
		 temp_all_DT_station.push_back(id.station());
		 temp_all_DT_layer.push_back(id.layer());
		 temp_all_DT_superlayer.push_back(id.superLayer());
		 temp_all_DT_wheel.push_back(id.wheel());
		 temp_all_DT_sector.push_back(id.sector());
	       }
	       
           } else {
               //Yes -> we just put the hit in the corresponding hit vector.
               const TrackingRecHit *rechitref = (const TrackingRecHit *) &(*itHit);
               it->second.push_back(rechitref);

	       temp_all_x.push_back(rechitref->localPosition().x()); 
	       temp_all_y.push_back(rechitref->localPosition().y()); 
	       temp_all_z.push_back(rechitref->localPosition().z()); 
	       temp_all_subdetid.push_back(ittrack->first->geographicalId().subdetId());
	       if(rechitref->geographicalId().subdetId() == MuonSubdetId::CSC){
		 CSCDetId id(rechitref->geographicalId().rawId());
		 temp_all_CSC_endcap.push_back(id.endcap());
		 temp_all_CSC_station.push_back(id.station());
		 temp_all_CSC_ringN.push_back(id.ring());
		 temp_all_CSC_chamber.push_back(id.chamber());
		 temp_all_CSC_layer.push_back(id.layer());
		 temp_all_DT_station.push_back(-9999);
		 temp_all_DT_layer.push_back(-9999);
		 temp_all_DT_superlayer.push_back(-9999);
		 temp_all_DT_wheel.push_back(-9999);
		 temp_all_DT_sector.push_back(-9999);
		 
	       }else if(rechitref->geographicalId().subdetId() == MuonSubdetId::DT){
		 DTWireId id(rechitref->geographicalId().rawId());
		 temp_all_CSC_endcap.push_back(-9999);
		 temp_all_CSC_station.push_back(-9999);
		 temp_all_CSC_ringN.push_back(-9999);
		 temp_all_CSC_chamber.push_back(-9999);
		 temp_all_CSC_layer.push_back(-9999);
		 temp_all_DT_station.push_back(id.station());
		 temp_all_DT_layer.push_back(id.layer());
		 temp_all_DT_superlayer.push_back(id.superLayer());
		 temp_all_DT_wheel.push_back(id.wheel());
		 temp_all_DT_sector.push_back(id.sector());
	       }
	       
           }
        }
     }


     Hit_DetAll_x.push_back(temp_all_x);
     Hit_DetAll_y.push_back(temp_all_y);
     Hit_DetAll_z.push_back(temp_all_z);
     Hit_DetAll_subdetid.push_back(temp_all_subdetid);
     Hit_DetAll_DT_station.push_back(temp_all_DT_station);
     Hit_DetAll_DT_layer.push_back(temp_all_DT_layer);
     Hit_DetAll_DT_superlayer.push_back(temp_all_DT_superlayer);
     Hit_DetAll_DT_wheel.push_back(temp_all_DT_wheel);
     Hit_DetAll_DT_sector.push_back(temp_all_DT_sector);
     Hit_DetAll_CSC_endcap.push_back(temp_all_CSC_endcap);
     Hit_DetAll_CSC_station.push_back(temp_all_CSC_station);
     Hit_DetAll_CSC_ringN.push_back(temp_all_CSC_ringN);
     Hit_DetAll_CSC_chamber.push_back(temp_all_CSC_chamber);
     Hit_DetAll_CSC_layer.push_back(temp_all_CSC_layer);

 
     //Now we do the extrapolations 
     for(auto it = DetRecHitMap.begin(); it != DetRecHitMap.end(); it++) {
       //std::cout << "Hit local position" << (*it).second.at(0)->localPosition() << std::endl;
       //std::cout << "Det ID" << (*it).first->geographicalId().subdetId() << std::endl;
         //Propagate
         std::pair<TrajectoryStateOnSurface, double> muonState = theService->propagator(propagator_)->propagateWithPath(outerTSOS, it->first->surface());
         if(muonState.first.isValid()){
	   temp_prop_x.push_back(muonState.first.localPosition().x());
	   temp_prop_y.push_back(muonState.first.localPosition().y());
	   temp_prop_z.push_back(muonState.first.localPosition().z());
	 }else{
	   temp_prop_x.push_back(-9999.);
	   temp_prop_y.push_back(-9999.);
	   temp_prop_z.push_back(-9999.);
	 }

	 //Global coords: muonState.first.globalPosition()

         //Here we have everything that we need:
         //1.- The geomDet in order to the local/global transformations
         //2.- The extrapolated state at the geomdet
         //3.- The vector of hits
     }    

     Hit_prop_x.push_back(temp_prop_x);
     Hit_prop_y.push_back(temp_prop_y);
     Hit_prop_z.push_back(temp_prop_z);

     temp_x.clear();
     temp_y.clear();
     temp_z.clear();
     temp_prop_x.clear();
     temp_prop_y.clear();
     temp_prop_z.clear();
     temp_subdetid.clear();
     temp_DT_station.clear();
     temp_DT_layer.clear();
     temp_DT_superlayer.clear();
     temp_DT_wheel.clear();
     temp_DT_sector.clear();
     temp_CSC_endcap.clear();
     temp_CSC_station.clear();
     temp_CSC_ringN.clear();
     temp_CSC_chamber.clear();
     temp_CSC_layer.clear();
     temp_all_x.clear();
     temp_all_y.clear();
     temp_all_z.clear();
     temp_all_prop_x.clear();
     temp_all_prop_y.clear();
     temp_all_prop_z.clear();
     temp_all_subdetid.clear();
     temp_all_DT_station.clear();
     temp_all_DT_layer.clear();
     temp_all_DT_superlayer.clear();
     temp_all_DT_wheel.clear();
     temp_all_DT_sector.clear();
     temp_all_CSC_endcap.clear();
     temp_all_CSC_station.clear();
     temp_all_CSC_ringN.clear();
     temp_all_CSC_chamber.clear();
     temp_all_CSC_layer.clear();


     // for(auto it = DetAllSegmentsMap.begin(); it != DetAllSegmentsMap.end(); it++) {
     //   for(int i=0; i<(int)(*it).second.size(); i++){
     // 	 temp_x.push_back((*it).second.at(i)->localPosition().x()); 
     // 	 temp_y.push_back((*it).second.at(i)->localPosition().y()); 
     // 	 temp_z.push_back((*it).second.at(i)->localPosition().z()); 
     // 	 temp_subdetid.push_back((*it).first->geographicalId().subdetId());
     // 	 if((*it).first->geographicalId().subdetId() == MuonSubdetId::CSC){
     // 	   CSCDetId id((*it).first->geographicalId().rawId());
     // 	   temp_CSC_endcap.push_back(id.endcap());
     // 	   temp_CSC_station.push_back(id.station());
     // 	   temp_CSC_ringN.push_back(id.ring());
     // 	   temp_CSC_chamber.push_back(id.chamber());
     // 	   temp_CSC_layer.push_back(id.layer());
     // 	   temp_DT_station.push_back(-9999);
     // 	   temp_DT_layer.push_back(-9999);
     // 	   temp_DT_superlayer.push_back(-9999);
     // 	   temp_DT_wheel.push_back(-9999);
     // 	   temp_DT_sector.push_back(-9999);

     // 	 }else if((*it).first->geographicalId().subdetId() == MuonSubdetId::DT){
     // 	   DTWireId id((*it).first->geographicalId().rawId());
     // 	   temp_CSC_endcap.push_back(-9999);
     // 	   temp_CSC_station.push_back(-9999);
     // 	   temp_CSC_ringN.push_back(-9999);
     // 	   temp_CSC_chamber.push_back(-9999);
     // 	   temp_CSC_layer.push_back(-9999);
     // 	   temp_DT_station.push_back(id.station());
     // 	   temp_DT_layer.push_back(id.layer());
     // 	   temp_DT_superlayer.push_back(id.superLayer());
     // 	   temp_DT_wheel.push_back(id.wheel());
     // 	   temp_DT_sector.push_back(id.sector());
     // 	 }
     //   }
     // }     



     /*	 
	 if (detid.det() == DetId::Tracker) {

	   const PixelGeomDetUnit* theGeomDet = dynamic_cast<const PixelGeomDetUnit*> (theTracker.idToDet(detid) );

	    //const StripGeomDetUnit* theGeomDet = dynamic_cast<const StripGeomDetUnit*>( theTracker.idToDet(detid) );
	   //GlobalPoint GP = theGeomDet->surface().toGlobal(Local3DPoint(lp));
	   GlobalPoint gp = theGeomDet->surface().toGlobal(Local3DPoint(lp.x(),lp.y(),lp.z()));

	   Hit_gpX.push_back(gp.x());
	   Hit_gpY.push_back(gp.y());
	   Hit_gpZ.push_back(gp.z());

   */
   }

     
   /////////////////////////////////////////////////////////////////////////////////////
   /////////////////////////////////// FILL THE TREE ///////////////////////////////////
   /////////////////////////////////////////////////////////////////////////////////////
   nMuons = iMuon;
   tree_out->Fill();
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
   Hit_Track_x.clear();
   Hit_Track_y.clear();
   Hit_Track_z.clear();
   Hit_prop_x.clear(); 
   Hit_prop_y.clear();
   Hit_prop_z.clear();
   Hit_Track_subdetid.clear();
   Hit_Track_DT_station.clear();
   Hit_Track_DT_layer.clear();
   Hit_Track_DT_superlayer.clear();
   Hit_Track_DT_wheel.clear();
   Hit_Track_DT_sector.clear();
   Hit_Track_CSC_endcap.clear();
   Hit_Track_CSC_station.clear();
   Hit_Track_CSC_ringN.clear();
   Hit_Track_CSC_chamber.clear();
   Hit_Track_CSC_layer.clear();
   Hit_DetAll_x.clear();
   Hit_DetAll_y.clear();
   Hit_DetAll_z.clear();
   Hit_DetAll_subdetid.clear();
   Hit_DetAll_DT_station.clear();
   Hit_DetAll_DT_layer.clear();
   Hit_DetAll_DT_superlayer.clear();
   Hit_DetAll_DT_wheel.clear();
   Hit_DetAll_DT_sector.clear();
   Hit_DetAll_CSC_endcap.clear();
   Hit_DetAll_CSC_station.clear();
   Hit_DetAll_CSC_ringN.clear();
   Hit_DetAll_CSC_chamber.clear();
   Hit_DetAll_CSC_layer.clear();
    
}

//=======================================================================================================================================================================================================================//
void RECOAnalysis::beginJob()
{

    // Output file definition
    // edm::Service<TFileService> *file_out;
    // output_filename = parameters.getParameter<std::string>("nameOfOutput");
    // tree_out = file_out->make<TTree>("Events","Events");
    file_out = new TFile(output_filename.c_str(), "RECREATE");
    file_out->cd();

    std::cout << "the file is created" << std::endl;
    
    // Output Tree definition
    tree_out = new TTree("Events", "Events");

    // ///////////////////////////////// EVENT INFO BRANCHES ///////////////////////////////

    tree_out->Branch("Event_event", &Event_event, "Event_event/I");
    tree_out->Branch("Event_run", &Event_run, "Event_run/I");
    tree_out->Branch("Event_luminosityBlock", &Event_luminosityBlock, "Event_luminosityBlock/I");

    // ////////////////////////////// MUON BRANCHES //////////////////////////////

    tree_out->Branch("nMuons", &nMuons, "nMuons/I");
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

    // ////////////////////////////// HIT BRANCHES //////////////////////////////

   tree_out->Branch("Hit_Track_x", "vector<vector<float> >", &Hit_Track_x);
   tree_out->Branch("Hit_Track_y", "vector<vector<float> >", &Hit_Track_y);
   tree_out->Branch("Hit_Track_z", "vector<vector<float> >", &Hit_Track_z);
   tree_out->Branch("Hit_Track_subdetid", "vector<vector<int> >", &Hit_Track_subdetid);
   tree_out->Branch("Hit_Track_DT_station", "vector<vector<int> >", &Hit_Track_DT_station);
   tree_out->Branch("Hit_Track_DT_layer", "vector<vector<int> >", &Hit_Track_DT_layer);
   tree_out->Branch("Hit_Track_DT_superlayer", "vector<vector<int> >", &Hit_Track_DT_superlayer);
   tree_out->Branch("Hit_Track_DT_wheel", "vector<vector<int> >", &Hit_Track_DT_wheel);
   tree_out->Branch("Hit_Track_DT_sector", "vector<vector<int> >", &Hit_Track_DT_sector);
   tree_out->Branch("Hit_Track_CSC_endcap", "vector<vector<int> >", &Hit_Track_CSC_endcap);
   tree_out->Branch("Hit_Track_CSC_station", "vector<vector<int> >", &Hit_Track_CSC_station);
   tree_out->Branch("Hit_Track_CSC_ringN", "vector<vector<int> >", &Hit_Track_CSC_ringN);
   tree_out->Branch("Hit_Track_CSC_chamber", "vector<vector<int> >", &Hit_Track_CSC_chamber);
   tree_out->Branch("Hit_Track_CSC_layer", "vector<vector<int> >", &Hit_Track_CSC_layer);

   tree_out->Branch("Hit_prop_x", "vector<vector<float> >", &Hit_prop_x);
   tree_out->Branch("Hit_prop_y", "vector<vector<float> >", &Hit_prop_y);
   tree_out->Branch("Hit_prop_z", "vector<vector<float> >", &Hit_prop_z);

   tree_out->Branch("Hit_DetAll_x", "vector<vector<float> >", &Hit_DetAll_x);
   tree_out->Branch("Hit_DetAll_y", "vector<vector<float> >", &Hit_DetAll_y);
   tree_out->Branch("Hit_DetAll_z", "vector<vector<float> >", &Hit_DetAll_z);
   tree_out->Branch("Hit_DetAll_subdetid", "vector<vector<int> >", &Hit_DetAll_subdetid);
   tree_out->Branch("Hit_DetAll_DT_station", "vector<vector<int> >", &Hit_DetAll_DT_station);
   tree_out->Branch("Hit_DetAll_DT_layer", "vector<vector<int> >", &Hit_DetAll_DT_layer);
   tree_out->Branch("Hit_DetAll_DT_superlayer", "vector<vector<int> >", &Hit_DetAll_DT_superlayer);
   tree_out->Branch("Hit_DetAll_DT_wheel", "vector<vector<int> >", &Hit_DetAll_DT_wheel);
   tree_out->Branch("Hit_DetAll_DT_sector", "vector<vector<int> >", &Hit_DetAll_DT_sector);
   tree_out->Branch("Hit_DetAll_CSC_endcap", "vector<vector<int> >", &Hit_DetAll_CSC_endcap);
   tree_out->Branch("Hit_DetAll_CSC_station", "vector<vector<int> >", &Hit_DetAll_CSC_station);
   tree_out->Branch("Hit_DetAll_CSC_ringN", "vector<vector<int> >", &Hit_DetAll_CSC_ringN);
   tree_out->Branch("Hit_DetAll_CSC_chamber", "vector<vector<int> >", &Hit_DetAll_CSC_chamber);
   tree_out->Branch("Hit_DetAll_CSC_layer", "vector<vector<int> >", &Hit_DetAll_CSC_layer);
 
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
