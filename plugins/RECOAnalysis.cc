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
#include "RecoMuon/Navigation/interface/MuonNavigationSchool.h"
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

std::vector<int> nHits_Track;
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

std::vector<int> nHits_DetAll;
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
   int iHit_Track = 0;
   int iHit_DetAll = 0;

   for (auto itmuon=muons->begin(); itmuon != muons->end(); itmuon++){

     if(!(itmuon->innerTrack().isNonnull())){continue;}
     if(itmuon->innerTrack()->pt() < 200.){continue;} // High-pT muons

     iMuon++;
     if(!(itmuon->globalTrack().isNonnull())){continue;}


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



     //-------------------------------------------------------------------------//
     //------------------------Added by Pablo-----------------------------------//
     //-------------------------------------------------------------------------//

     //We build the transient track of the tracker track and get its final state on surface
     reco::TransientTrack trackinner(itmuon->innerTrack(), &*theService->magneticField(), theService->trackingGeometry());     
     TrajectoryStateOnSurface outerTSOS = trackinner.outermostMeasurementState(); 

     //We also build the transient track of the global track
     //reco::TransientTrack track(itmuon->globalTrack(), &*theService->magneticField(), theService->trackingGeometry());     

     trackingRecHit_iterator lastHit = trackinner.recHitsEnd() - 1;       
     DetId outerDetId((*lastHit)->geographicalId());
       
     // Get the layer on which the seed relies
     const DetLayer *initialLayer = theService->detLayerGeometry()->idToLayer(outerDetId);
       
     PropagationDirection detLayerOrder = oppositeToMomentum;
       
       // ask for compatible layers
     std::vector<const DetLayer *> detLayers;
     detLayers = theService->muonNavigationSchool()->compatibleLayers(*initialLayer, *outerTSOS.freeState(), detLayerOrder);


     std::vector<const GeomDet*> ExtrapolationDets;


     for(auto it = detLayers.begin(); it != detLayers.end(); it++){

       if(!(GEOMDET->geographicalId().det() == DetId::Muon)) continue;
       if(!(GEOMDET->geographicalId().subdetId() == MuonSubdetId::DT || GEOMDET->geographicalId().subdetId() == MuonSubdetId::CSC)) continue;

       //std::pair<TrajectoryStateOnSurface, double> muonState = theService->propagator(propagator_)->propagateWithPath(outerTSOS, (*it)->surface());
       std::pair<TrajectoryStateOnSurface, double> muonState = theService->propagator(propagator_)->propagateWithPath(outerTSOS, (*GEOMDET)->surface());

       if(muonState.first.isValid()){

	 std::cout << muonState.first.globalPosition() << std::endl;
	 ExtrapolationDets.push_back(GEOMDET);
       }
     }



     std::map<const GeomDet*, std::vector<const TrackingRecHit *> > DetAllSegmentsMap; 
     
     for(auto itGeomDet = ExtrapolationDets.begin(); itGeomDet!= ExtrapolationDets.end(); itGeomDet++){

       for (auto itHit = dtSegments->begin(); itHit != dtSegments->end(); itHit++) {
	 //Only valid hits     
	 if(!itHit->isValid()) continue;
	 DetId myDet = itHit->geographicalId();
	 if(myDet != itGeomDet) continue;
           //Get the GeomDet associated to this DetIt 
           std::map<const GeomDet*, std::vector<const TrackingRecHit *> >::iterator it = DetAllSegmentsMap.find(itGeomDet);
           if(it == DetAllSegmentsMap.end()) {
               //No -> we create a pair of GeomDet and vector of hits, and put the hit in the vector.
               std::vector<const TrackingRecHit *> trhit;
               const TrackingRecHit *rechitref = (const TrackingRecHit *)&(*itHit);
               trhit.push_back(rechitref);
               DetAllSegmentsMap.insert(std::pair<const GeomDet*, std::vector<const TrackingRecHit *> > (itGeomDet, trhit));
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
           if(myDet != itGeomDet) continue;

           //Get the GeomDet associated to this DetIt 
           std::map<const GeomDet*, std::vector<const TrackingRecHit *> >::iterator it = DetAllSegmentsMap.find(itGeomDet);
           if(it == DetAllSegmentsMap.end()) {
               //No -> we create a pair of GeomDet and vector of hits, and put the hit in the vector.
               std::vector<const TrackingRecHit *> trhit;
               const TrackingRecHit *rechitref = (const TrackingRecHit *)&(*itHit);
               trhit.push_back(rechitref);
               DetAllSegmentsMap.insert(std::pair <const GeomDet*, std::vector<const TrackingRecHit *> > (itGeomDet, trhit));
	       
           } else {
               //Yes -> we just put the hit in the corresponding hit vector.
               const TrackingRecHit *rechitref = (const TrackingRecHit *) &(*itHit);
               it->second.push_back(rechitref);
	       
           }
        }
       


     }
   }

     
   /////////////////////////////////////////////////////////////////////////////////////
   /////////////////////////////////// FILL THE TREE ///////////////////////////////////
   /////////////////////////////////////////////////////////////////////////////////////
   nMuons = iMuon;
   tree_out->Fill();
   // Muon_GlbTrack_pt.clear();
   // Muon_GlbTrack_eta.clear();
   // Muon_GlbTrack_phi.clear();
   // Muon_GlbTrack_charge.clear();
   // Muon_GlbTrack_ptErr.clear();
   // Muon_GlbTrack_Chindf.clear();
   // Muon_InnerTrack_pt.clear();
   // Muon_InnerTrack_eta.clear();
   // Muon_InnerTrack_phi.clear();
   // Muon_InnerTrack_charge.clear();
   // Muon_InnerTrack_ptErr.clear();
   // Muon_InnerTrack_Chindf.clear();
   // Muon_TunePTrack_pt.clear();
   // Muon_TunePTrack_eta.clear();
   // Muon_TunePTrack_phi.clear();
   // Muon_TunePTrack_charge.clear();
   // Muon_TunePTrack_ptErr.clear();
   // Muon_TunePTrack_Chindf.clear();
   // nHits_Track.clear();
   // Hit_Track_x.clear();
   // Hit_Track_y.clear();
   // Hit_Track_z.clear();
   // Hit_prop_x.clear(); 
   // Hit_prop_y.clear();
   // Hit_prop_z.clear();
   // Hit_Track_subdetid.clear();
   // Hit_Track_DT_station.clear();
   // Hit_Track_DT_layer.clear();
   // Hit_Track_DT_superlayer.clear();
   // Hit_Track_DT_wheel.clear();
   // Hit_Track_DT_sector.clear();
   // Hit_Track_CSC_endcap.clear();
   // Hit_Track_CSC_station.clear();
   // Hit_Track_CSC_ringN.clear();
   // Hit_Track_CSC_chamber.clear();
   // Hit_Track_CSC_layer.clear();
   // nHits_DetAll.clear();
   // Hit_DetAll_x.clear();
   // Hit_DetAll_y.clear();
   // Hit_DetAll_z.clear();
   // Hit_DetAll_subdetid.clear();
   // Hit_DetAll_DT_station.clear();
   // Hit_DetAll_DT_layer.clear();
   // Hit_DetAll_DT_superlayer.clear();
   // Hit_DetAll_DT_wheel.clear();
   // Hit_DetAll_DT_sector.clear();
   // Hit_DetAll_CSC_endcap.clear();
   // Hit_DetAll_CSC_station.clear();
   // Hit_DetAll_CSC_ringN.clear();
   // Hit_DetAll_CSC_chamber.clear();
   // Hit_DetAll_CSC_layer.clear();
    
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

   // tree_out->Branch("nHits_Track", "vector<int>", &nHits_Track);
   // tree_out->Branch("Hit_Track_x", "vector<vector<float> >", &Hit_Track_x);
   // tree_out->Branch("Hit_Track_y", "vector<vector<float> >", &Hit_Track_y);
   // tree_out->Branch("Hit_Track_z", "vector<vector<float> >", &Hit_Track_z);
   // tree_out->Branch("Hit_Track_subdetid", "vector<vector<int> >", &Hit_Track_subdetid);
   // tree_out->Branch("Hit_Track_DT_station", "vector<vector<int> >", &Hit_Track_DT_station);
   // tree_out->Branch("Hit_Track_DT_layer", "vector<vector<int> >", &Hit_Track_DT_layer);
   // tree_out->Branch("Hit_Track_DT_superlayer", "vector<vector<int> >", &Hit_Track_DT_superlayer);
   // tree_out->Branch("Hit_Track_DT_wheel", "vector<vector<int> >", &Hit_Track_DT_wheel);
   // tree_out->Branch("Hit_Track_DT_sector", "vector<vector<int> >", &Hit_Track_DT_sector);
   // tree_out->Branch("Hit_Track_CSC_endcap", "vector<vector<int> >", &Hit_Track_CSC_endcap);
   // tree_out->Branch("Hit_Track_CSC_station", "vector<vector<int> >", &Hit_Track_CSC_station);
   // tree_out->Branch("Hit_Track_CSC_ringN", "vector<vector<int> >", &Hit_Track_CSC_ringN);
   // tree_out->Branch("Hit_Track_CSC_chamber", "vector<vector<int> >", &Hit_Track_CSC_chamber);
   // tree_out->Branch("Hit_Track_CSC_layer", "vector<vector<int> >", &Hit_Track_CSC_layer);

   // tree_out->Branch("Hit_prop_x", "vector<vector<float> >", &Hit_prop_x);
   // tree_out->Branch("Hit_prop_y", "vector<vector<float> >", &Hit_prop_y);
   // tree_out->Branch("Hit_prop_z", "vector<vector<float> >", &Hit_prop_z);

   // tree_out->Branch("nHits_DetAll", "vector<int>", &nHits_DetAll);
   // tree_out->Branch("Hit_DetAll_x", "vector<vector<float> >", &Hit_DetAll_x);
   // tree_out->Branch("Hit_DetAll_y", "vector<vector<float> >", &Hit_DetAll_y);
   // tree_out->Branch("Hit_DetAll_z", "vector<vector<float> >", &Hit_DetAll_z);
   // tree_out->Branch("Hit_DetAll_subdetid", "vector<vector<int> >", &Hit_DetAll_subdetid);
   // tree_out->Branch("Hit_DetAll_DT_station", "vector<vector<int> >", &Hit_DetAll_DT_station);
   // tree_out->Branch("Hit_DetAll_DT_layer", "vector<vector<int> >", &Hit_DetAll_DT_layer);
   // tree_out->Branch("Hit_DetAll_DT_superlayer", "vector<vector<int> >", &Hit_DetAll_DT_superlayer);
   // tree_out->Branch("Hit_DetAll_DT_wheel", "vector<vector<int> >", &Hit_DetAll_DT_wheel);
   // tree_out->Branch("Hit_DetAll_DT_sector", "vector<vector<int> >", &Hit_DetAll_DT_sector);
   // tree_out->Branch("Hit_DetAll_CSC_endcap", "vector<vector<int> >", &Hit_DetAll_CSC_endcap);
   // tree_out->Branch("Hit_DetAll_CSC_station", "vector<vector<int> >", &Hit_DetAll_CSC_station);
   // tree_out->Branch("Hit_DetAll_CSC_ringN", "vector<vector<int> >", &Hit_DetAll_CSC_ringN);
   // tree_out->Branch("Hit_DetAll_CSC_chamber", "vector<vector<int> >", &Hit_DetAll_CSC_chamber);
   // tree_out->Branch("Hit_DetAll_CSC_layer", "vector<vector<int> >", &Hit_DetAll_CSC_layer);
 
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
