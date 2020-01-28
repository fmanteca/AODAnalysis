#include <memory>
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

//framework includes
//#include "DataFormats/Common/interface/EDProduct.h"

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
#include "Geometry/Records/interface/MuonGeometryRecord.h"
#include "Geometry/CSCGeometry/interface/CSCGeometry.h"
#include "Geometry/DTGeometry/interface/DTGeometry.h"
#include "Geometry/RPCGeometry/interface/RPCGeometry.h"
#include "DataFormats/MuonDetId/interface/MuonSubdetId.h"
#include "DataFormats/MuonDetId/interface/RPCDetId.h"
#include "DataFormats/MuonDetId/interface/DTWireId.h"

#include "MagneticField/Engine/interface/MagneticField.h"
#include "RecoMuon/TrackingTools/interface/MuonServiceProxy.h"
#include "TrackingTools/GeomPropagators/interface/Propagator.h"
#include "TrackingTools/Records/interface/TrackingComponentsRecord.h"

#include <string>
#include <iostream>
#include <vector>
#include <algorithm>

#include "TLorentzVector.h"
#include "TTree.h"
#include "TFile.h"

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


//-> MUON SELECTION
Int_t   nMuons;

std::vector<float> Muon_pt;
std::vector<float> Muon_eta;
std::vector<float> Muon_phi;
Int_t   nHits;
Float_t Hit_x;
Float_t Hit_y;
Float_t Hit_z;


/////////////////////////////////////// OUTPUT //////////////////////////////////////

TFile *file_out;
TTree *tree_out;


//=======================================================================================================================================================================================================================//
class AODAnalysis : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
  static const int N_MAX_STORED = 10;
  static const int N_MAX_STORED_HIT = 1000;

   public:
      typedef struct {
    
       int n_;
       int muonLink_[N_MAX_STORED_HIT];
    
       int system_[N_MAX_STORED_HIT];
       int endcap_[N_MAX_STORED_HIT];
       int station_[N_MAX_STORED_HIT];
       int ring_[N_MAX_STORED_HIT];
       int chamber_[N_MAX_STORED_HIT];
       int layer_[N_MAX_STORED_HIT];
       int superLayer_[N_MAX_STORED_HIT];
       int wheel_[N_MAX_STORED_HIT];
       int sector_[N_MAX_STORED_HIT];
    
       float gpX_[N_MAX_STORED_HIT];
       float gpY_[N_MAX_STORED_HIT];
       float gpZ_[N_MAX_STORED_HIT];
       // float gpEta_[N_MAX_STORED_HIT];
       // float gpPhi_[N_MAX_STORED_HIT];
       float lpX_[N_MAX_STORED_HIT];
       float lpY_[N_MAX_STORED_HIT];
       float lpZ_[N_MAX_STORED_HIT];
    

     } storage_hit;

      AODAnalysis::storage_hit storageRecMuon_;
      explicit AODAnalysis(const edm::ParameterSet&);
      ~AODAnalysis();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


   private:
      virtual void beginJob() override;
      virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
      virtual void endJob() override;

      std::string output_filename;
      edm::ParameterSet parameters;

      edm::EDGetTokenT<edm::View<reco::Muon> > theMuonCollection;
  


};
//=======================================================================================================================================================================================================================//




//=======================================================================================================================================================================================================================//
AODAnalysis::AODAnalysis(const edm::ParameterSet& iConfig)
{

   usesResource("TFileService");
   
   parameters = iConfig;

   theMuonCollection = consumes<edm::View<reco::Muon> >  (parameters.getParameter<edm::InputTag>("MuonCollection"));
   //theGeneralTrackCollection = consumes<edm::View<reco::Track> > (parameters.getParameter<edm::InputTag>("GeneralTrackCollection"));
   
}

//=======================================================================================================================================================================================================================//
AODAnalysis::~AODAnalysis()
{

}

//=======================================================================================================================================================================================================================//
void AODAnalysis::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{

  //magnetic field information
  // edm::ESHandle<MagneticField> field;
  // edm::ESHandle<GlobalTrackingGeometry> globalTrackingGeometry;
  // iSetup.get<IdealMagneticFieldRecord>().get(field);
  // iSetup.get<GlobalTrackingGeometryRecord>().get(globalTrackingGeometry);
  // iSetup.get<TrackingComponentsRecord>().get( PropagatorSource_, thePropagator );
  // theField = &*field;


   /////////////////////////////////////////////////////////////////////////////////////
   ///////////////////////////////////// MAIN CODE /////////////////////////////////////
   /////////////////////////////////////////////////////////////////////////////////////



   //////////////////////////////// GET THE COLLECTIONS ////////////////////////////////
   
   edm::Handle<edm::View<reco::Muon> > muons;
   //edm::Handle<edm::View<reco::Track> > generalTracks;



  iEvent.getByToken(theMuonCollection, muons);
   //   iEvent.getByToken(theGeneralTrackCollection, generalTracks);



   /////////////////////////////////// EVENT INFO //////////////////////////////////////


   Event_event = iEvent.id().event();
   Event_run = iEvent.id().run();
   Event_luminosityBlock = iEvent.id().luminosityBlock();



   ////////////////////////////// MUON FEATURES //////////////////////////////


   nMuons = muons->size();
 
   /////////////////////////////////////////////////////////////////////////////////////
   ///////////////////////////////// GET MUON VARIABLES ////////////////////////////////
   /////////////////////////////////////////////////////////////////////////////////////

   int iMuonHit = 0;
   //int iTrackHit= 0;
   //int numTracks= 0;

   for (auto itmuon=muons->begin(); itmuon != muons->end(); itmuon++){

     int iMuon = itmuon - muons->begin();

     std::cout << "HOLA MUON " << iMuon << std::endl;

     Muon_pt.push_back(itmuon->pt());
     Muon_eta.push_back(itmuon->eta());
     Muon_phi.push_back(itmuon->phi());
     
     for (auto itHit = itmuon->outerTrack()->recHitsBegin(); itHit != itmuon->outerTrack()->recHitsEnd(); itHit++) {

       int iHit = itHit - itmuon->outerTrack()->recHitsBegin();

       DetId detid = (*itHit)->geographicalId(); 

       //if (detid.det() != DetId::Muon && detid.det() != DetId::Tracker) {
       //REMOVE MUON ID = TRACKER
       if (detid.det() != DetId::Muon) {
	 continue;
       }

       LocalPoint lp = (*itHit)->localPosition();
       float lpX= lp.x();
       float lpY = lp.y();
       float lpZ = lp.z();
       float gpX;
       float gpY;
       float gpZ;

       if (detid.det() == DetId::Muon) {


       std::cout << " - HOLA HIT" << std::endl;
       std::cout << "   THIS DETECTOR: " << detid.det() << " (MuonDet=2)" << std::endl;
       std::cout << "   THIS SUBDETECTOR: " << detid.subdetId() << " (DT=1, CSC=2, RPC=3)" << std::endl;
       std::cout << "   HITS LOCAL COORDS (x,y,z): " << "(" << lpX << "," << lpY << "," << lpZ << ")" << std::endl; 
	 
     	 int systemMuon  = detid.subdetId(); // 1 DT; 2 CSC; 3 RPC
     	 int endcap= -999;
     	 int station= -999;
     	 int ring= -999;
     	 int chamber= -999;
     	 int layer= -999;
     	 int superLayer  = -999;
     	 int wheel = -999;
     	 int sector = -999;
     	 if ( systemMuon == MuonSubdetId::CSC) {
     	   CSCDetId id(detid.rawId());
     	   endcap= id.endcap();
     	   station= id.station();
     	   ring= id.ring();
     	   chamber= id.chamber();
     	   layer= id.layer();

	   edm::ESHandle<CSCGeometry> theCSCGeometry;
	   iSetup.get<MuonGeometryRecord>().get(theCSCGeometry);

	   const CSCGeometry& theCSCMuon(*theCSCGeometry);
	   //DetId theDetUnitId((*itHit)->detUnitId());
	   //const GeomDetUnit *theDet = theCSCMuon.idToDetUnit(theDetUnitId);
	   const GeomDetUnit *theDet = theCSCMuon.idToDetUnit(detid);
	   const BoundPlane& bSurface = theDet->surface();
	   gpX = bSurface.toGlobal((*itHit)->localPosition()).x();
	   gpY = bSurface.toGlobal((*itHit)->localPosition()).y();
	   gpZ = bSurface.toGlobal((*itHit)->localPosition()).z();

     	   printf("   CSC\t[endcap][station][ringN][chamber][layer]:[%d][%d][%d][%d][%d]\t",
     			     endcap, station, ring, chamber, layer);

     	 }
     	 else if ( systemMuon == MuonSubdetId::DT ) {
     	   DTWireId id(detid.rawId());
     	   station= id.station();
     	   layer= id.layer();
     	   superLayer= id.superLayer();
     	   wheel= id.wheel();
     	   sector= id.sector();
     	   printf("   DT \t[station][layer][superlayer]:[%d][%d][%d]\n", station,layer,superLayer);
	   
     	 }
     	 else if ( systemMuon == MuonSubdetId::RPC) {
     	   RPCDetId id(detid.rawId());
     	   station= id.station();
     	   printf("   RPC\t[station]:[%d]\n", station);
     	 }
       
       
	 
       	 storageRecMuon_.muonLink_[iMuonHit]= iMuon;
       	 storageRecMuon_.system_[iMuonHit]= systemMuon;
       	 storageRecMuon_.endcap_[iMuonHit]= endcap;
       	 storageRecMuon_.station_[iMuonHit]= station;
       	 storageRecMuon_.ring_[iMuonHit]= ring;
       	 storageRecMuon_.chamber_[iMuonHit]= chamber;
       	 storageRecMuon_.layer_[iMuonHit]= layer;
       	 storageRecMuon_.superLayer_[iMuonHit]= superLayer;
       	 storageRecMuon_.wheel_[iMuonHit]= wheel;
       	 storageRecMuon_.sector_[iMuonHit]= sector;
	 
       	 storageRecMuon_.gpX_[iMuonHit]= gpX;
       	 storageRecMuon_.gpY_[iMuonHit]= gpY;
       	 storageRecMuon_.gpZ_[iMuonHit]= gpZ;
       	 // storageRecMuon_.gpEta_[iMuonHit]= gpRecEta;
       	 // storageRecMuon_.gpPhi_[iMuonHit]= gpRecPhi;
       	 storageRecMuon_.lpX_[iMuonHit]= lpX;
       	 storageRecMuon_.lpY_[iMuonHit]= lpY;
       	 storageRecMuon_.lpZ_[iMuonHit]= lpZ;
       	 iMuonHit++;
	 
       }
	 else{continue;}


       // else if (detid.det() == DetId::Tracker) {
	 
       // 	 if (debug_) printf("Tracker\n");

       // 	 StoreTrackerRecHits(detid, iTrack, iTrackHit);

       // 	 storageTrackHit_.gpX_[iTrackHit]= gpRecX;
       // 	 storageTrackHit_.gpY_[iTrackHit]= gpRecY;
       // 	 storageTrackHit_.gpZ_[iTrackHit]= gpRecZ;
       // 	 storageTrackHit_.gpEta_[iTrackHit]= gpRecEta;
       // 	 storageTrackHit_.gpPhi_[iTrackHit]= gpRecPhi;
       // 	 storageTrackHit_.lpX_[iTrackHit]= lpX;
       // 	 storageTrackHit_.lpY_[iTrackHit]= lpY;
       // 	 storageTrackHit_.lpZ_[iTrackHit]= lpZ;
       // 	 iTrackHit++;
       // }
       // else printf("THIS CAN NOT HAPPEN\n");
       
       // trkExtrap(detid, numTracks, iTrack, iRec, recoStart, lp, trackExtrap);
       // numTracks++;

       // if (debug_) printf("\tLocal Positon:  \tx = %2.2f\ty = %2.2f\tz = %2.2f\n",lpX, lpY, lpZ);
       // if (debug_) printf("\tGlobal Position: \tx = %6.2f\ty = %6.2f\tz = %6.2f\teta = %4.2f\tphi = %3.2f\n",
       // 			  gpRecX,gpRecY,gpRecZ,gpRecEta,gpRecPhi);
       

     }
   }
     
   storageRecMuon_.n_ = iMuonHit; 
   //storageTrackHit_.n_= iTrackHit;
   //trackExtrap.n_ = numTracks;


 




   

   
   /////////////////////////////////////////////////////////////////////////////////////
   /////////////////////////////////// FILL THE TREE ///////////////////////////////////
   /////////////////////////////////////////////////////////////////////////////////////
   tree_out->Fill();
   Muon_pt.clear();
   Muon_eta.clear();
   Muon_phi.clear();


}

//=======================================================================================================================================================================================================================//
void AODAnalysis::beginJob()
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



    // ////////////////////////////// PRIMARY VERTEX BRANCHES //////////////////////////////

    tree_out->Branch("nMuons", &nMuons, "nMuons/I");
    tree_out->Branch("Muon_pt", "vector<float>", &Muon_pt);
    tree_out->Branch("Muon_eta", "vector<float>", &Muon_eta);
    tree_out->Branch("Muon_phi", "vector<float>", &Muon_phi);

    tree_out -> Branch("recHits", &storageRecMuon_, 
    		       "n_/I:"
    		       "muonLink_[1000]/I:"
	   
    		       "system_[1000]/I:"
    		       "endcap_[1000]/I:"
    		       "station_[1000]/I:"
    		       "ring_[1000]/I:"
    		       "chamber_[1000]/I:"
    		       "layer_[1000]/I:"
    		       "superLayer_[1000]/I:"
    		       "wheel_[1000]/I:"
    		       "sector_[1000]/I:"
	   
    		       "gpX_[1000]/F:"
    		       "gpY_[1000]/F:"
    		       "gpZ_[1000]/F:"
    		       // "gpEta_[1000]/F:"
    		       // "gpPhi_[1000]/F:"
    		       "lpX_[1000]/F:"
    		       "lpY_[1000]/F:"
    		       "lpZ_[1000]/F"
    					  );


}

//=======================================================================================================================================================================================================================//
void AODAnalysis::endJob() 
{


    std::cout << "The event is writen" << std::endl;
    file_out->cd();
    tree_out->Write();
    file_out->Close();

}

//=======================================================================================================================================================================================================================//
void AODAnalysis::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

DEFINE_FWK_MODULE(AODAnalysis);

