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
std::vector<float> Muon_TunePTrack_pt;
std::vector<float> Muon_TunePTrack_eta;
std::vector<float> Muon_TunePTrack_phi;
std::vector<float> Muon_TunePTrack_charge;
std::vector<float> Muon_TunePTrack_ptErr;
std::vector<float> Muon_TunePTrack_Chindf;

std::vector<std::vector<int> > nHits_GeomDet;
std::vector<int> nGeomDets;

std::vector<std::vector<int> > GeomDet_subdetid;
std::vector<std::vector<int> > GeomDet_subdet_place_one;
std::vector<std::vector<int> > GeomDet_subdet_place_two;
std::vector<std::vector<int> > GeomDet_subdet_place_three;
std::vector<std::vector<int> > GeomDet_subdet_place_four;
std::vector<std::vector<int> > GeomDet_subdet_place_five;

std::vector<std::vector<float> > Prop_x;
std::vector<std::vector<float> > Prop_y;
std::vector<std::vector<float> > Prop_z;

std::vector<std::vector<std::vector<float> > > Hit_x;
std::vector<std::vector<std::vector<float> > > Hit_y;
std::vector<std::vector<std::vector<float> > > Hit_z;



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

     if(!(itmuon->innerTrack().isNonnull())){continue;}
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

     std::vector<float> temp_prop_x;
     std::vector<float> temp_prop_y;
     std::vector<float> temp_prop_z;
     std::vector<int> temp_subdetid;
     std::vector<int> temp_subdet_place_one;
     std::vector<int> temp_subdet_place_two;
     std::vector<int> temp_subdet_place_three;
     std::vector<int> temp_subdet_place_four;
     std::vector<int> temp_subdet_place_five;

     std::vector<std::vector<float> > temp_hit_x_in_allGeomDet;
     std::vector<std::vector<float> > temp_hit_y_in_allGeomDet;
     std::vector<std::vector<float> > temp_hit_z_in_allGeomDet;

     std::vector<float> temp_hit_x_in_eachGeomDet;
     std::vector<float> temp_hit_y_in_eachGeomDet;
     std::vector<float> temp_hit_z_in_eachGeomDet;
     std::vector<int> temp_nHits_geomdet;


     // Build the transient track of the tracker track and get its final state on surface
     reco::TransientTrack trackinner(itmuon->innerTrack(), &*theService->magneticField(), theService->trackingGeometry());     
     TrajectoryStateOnSurface outerTSOS = trackinner.outermostMeasurementState(); 

     std::map<const GeomDet*, std::vector<TrackingRecHit *> > DetAllSegmentsMap; 

     
     // Loop over dtSegments

     for (auto itHit = dtSegments->begin(); itHit != dtSegments->end(); itHit++) {

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
     
     int nGeomDet = 0;
 
     //Now we do the extrapolations 

     for(auto it = DetAllSegmentsMap.begin(); it != DetAllSegmentsMap.end(); it++) {

         //Propagate
         std::pair<TrajectoryStateOnSurface, double> muonState = theService->propagator(propagator_)->propagateWithPath(outerTSOS, it->first->surface());

     	 // Store the hit if the extrapolation is valid
         if(muonState.first.isValid()){

	   nGeomDet++;

     	   temp_prop_x.push_back(muonState.first.localPosition().x());
     	   temp_prop_y.push_back(muonState.first.localPosition().y());
     	   temp_prop_z.push_back(muonState.first.localPosition().z());
	   temp_subdetid.push_back((*it).first->geographicalId().subdetId());

	   if((*it).first->geographicalId().subdetId() == MuonSubdetId::CSC){
	     CSCDetId id((*it).first->geographicalId().rawId());
	     temp_subdet_place_one.push_back(id.endcap());
	     temp_subdet_place_two.push_back(id.station());
	     temp_subdet_place_three.push_back(id.ring());
	     temp_subdet_place_four.push_back(id.chamber());
	     temp_subdet_place_five.push_back(id.layer());
	     
	   }else if((*it).first->geographicalId().subdetId() == MuonSubdetId::DT){
	     DTWireId id((*it).first->geographicalId().rawId());
	     temp_subdet_place_one.push_back(id.station());
	     temp_subdet_place_two.push_back(id.layer());
	     temp_subdet_place_three.push_back(id.superLayer());
	     temp_subdet_place_four.push_back(id.wheel());
	     temp_subdet_place_five.push_back(id.sector());
	   }
	   
	   int nHits_geomdet = 0;

	   for(int i=0; i<(int)(*it).second.size(); i++){
	     nHits_geomdet++;
	     std::cout<<(*it).second.at(i)->localPosition().x()<<std::endl;
	     temp_hit_x_in_eachGeomDet.push_back((*it).second.at(i)->localPosition().x()); 
	     temp_hit_y_in_eachGeomDet.push_back((*it).second.at(i)->localPosition().y()); 
	     temp_hit_z_in_eachGeomDet.push_back((*it).second.at(i)->localPosition().z()); 
	   }

	   temp_nHits_geomdet.push_back(nHits_geomdet);
	   temp_hit_x_in_allGeomDet.push_back(temp_hit_x_in_eachGeomDet);
	   temp_hit_y_in_allGeomDet.push_back(temp_hit_y_in_eachGeomDet);
	   temp_hit_z_in_allGeomDet.push_back(temp_hit_z_in_eachGeomDet);
	   temp_hit_x_in_eachGeomDet.clear();
	   temp_hit_y_in_eachGeomDet.clear();
	   temp_hit_z_in_eachGeomDet.clear();
	   
	   
     	 }
	 
     	 //Global coords: muonState.first.globalPosition()

         //Here we have everything that we need:
         //1.- The geomDet in order to the local/global transformations
         //2.- The extrapolated state at the geomdet
         //3.- The vector of hits
     }    
     
     // Counters:
     nGeomDets.push_back(nGeomDet);
     nHits_GeomDet.push_back(temp_nHits_geomdet);
     temp_nHits_geomdet.clear();

     // Hit & propagation info
     Hit_x.push_back(temp_hit_x_in_allGeomDet);
     Hit_y.push_back(temp_hit_y_in_allGeomDet);
     Hit_z.push_back(temp_hit_z_in_allGeomDet);
     temp_hit_x_in_allGeomDet.clear();
     temp_hit_y_in_allGeomDet.clear();
     temp_hit_z_in_allGeomDet.clear();

     Prop_x.push_back(temp_prop_x);
     Prop_y.push_back(temp_prop_y);
     Prop_z.push_back(temp_prop_z);
     temp_prop_x.clear();
     temp_prop_y.clear();
     temp_prop_z.clear();

     GeomDet_subdetid.push_back(temp_subdetid);
     GeomDet_subdet_place_one.push_back(temp_subdet_place_one);
     GeomDet_subdet_place_two.push_back(temp_subdet_place_two);
     GeomDet_subdet_place_three.push_back(temp_subdet_place_three);
     GeomDet_subdet_place_four.push_back(temp_subdet_place_four);
     GeomDet_subdet_place_five.push_back(temp_subdet_place_five);
     temp_subdetid.clear();
     temp_subdet_place_one.clear();
     temp_subdet_place_two.clear();
     temp_subdet_place_three.clear();
     temp_subdet_place_four.clear();
     temp_subdet_place_five.clear();

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

   nHits_GeomDet.clear();
   nGeomDets.clear();
   Hit_x.clear();
   Hit_y.clear();
   Hit_z.clear();
   Prop_x.clear();
   Prop_y.clear();
   Prop_z.clear();
   GeomDet_subdetid.clear();
   GeomDet_subdet_place_one.clear();
   GeomDet_subdet_place_two.clear();
   GeomDet_subdet_place_three.clear();
   GeomDet_subdet_place_four.clear();
   GeomDet_subdet_place_five.clear();
   
}

//=======================================================================================================================================================================================================================//
void RECOAnalysis::beginJob()
{
    gInterpreter->GenerateDictionary("vector<vector<vector<float> > >", "vector");

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

    // ////////////////////////////// HIT & EXTRAPOLATION BRANCHES //////////////////////////////

    tree_out->Branch("nGeomDets", "vector<int>", &nGeomDets);
    tree_out->Branch("nHits_GeomDet", "vector<vector<int> >", &nHits_GeomDet);
    tree_out->Branch("GeomDet_subdetid", "vector<vector<int> >", &GeomDet_subdetid);
    tree_out->Branch("GeomDet_subdet_place_one", "vector<vector<int> >", &GeomDet_subdet_place_one);
    tree_out->Branch("GeomDet_subdet_place_two", "vector<vector<int> >", &GeomDet_subdet_place_two);
    tree_out->Branch("GeomDet_subdet_place_three", "vector<vector<int> >", &GeomDet_subdet_place_three);
    tree_out->Branch("GeomDet_subdet_place_four", "vector<vector<int> >", &GeomDet_subdet_place_four);
    tree_out->Branch("GeomDet_subdet_place_five", "vector<vector<int> >", &GeomDet_subdet_place_five);
    
    tree_out->Branch("Prop_x", "vector<vector<float> >", &Prop_x);
    tree_out->Branch("Prop_y", "vector<vector<float> >", &Prop_y);
    tree_out->Branch("Prop_z", "vector<vector<float> >", &Prop_z);
    
    tree_out->Branch("Hit_x", "vector<vector<vector<float> > >", &Hit_x);
    tree_out->Branch("Hit_y", "vector<vector<vector<float> > >", &Hit_y);
    tree_out->Branch("Hit_z", "vector<vector<vector<float> > >", &Hit_z);
 
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
