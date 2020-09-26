/*
 * CMSSW plugin that performs a Window-based inference of networks using RecHits and produced PF candidates.
 *
 * Authors: Gerrit Van Onsem <Gerrit.Van.Onsem@cern.ch>
 *          Marcel Rieger <marcel.rieger@cern.ch>
 */

#include <memory>
#include <iostream>
#include <fstream>

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/Exception.h"

#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/HGCRecHit/interface/HGCRecHitCollections.h"
#include "DataFormats/CaloRecHit/interface/CaloClusterFwd.h"
#include "SimDataFormats/CaloAnalysis/interface/SimCluster.h"

#include "RecoLocalCalo/HGCalRecAlgos/interface/RecHitTools.h"

#include "RecoHGCal/GraphReco/interface/InferenceWindow.h"

#include "DataFormats/Common/interface/View.h"

#include "DataFormats/ParticleFlowCandidate/interface/PFCandidateFwd.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h" 
#include "DataFormats/Candidate/interface/LeafCandidate.h"

#include "DataFormats/HGCalReco/interface/Trackster.h"

using namespace ticl;

// macros for simplified logs
// message logger disabled for the moment
// #define INFO edm::LogInfo("peprCandidateFromHitProducer")
// #define WARNING edm::LogWarning("peprCandidateFromHitProducer")
// #define ERROR edm::LogError("peprCandidateFromHitProducer")
#define INFO std::cout << "peprCandidateFromHitProducer INFO   : "
#define WARNING std::cout << "peprCandidateFromHitProducer WARNING: "
#define ERROR std::cout << "peprCandidateFromHitProducer ERROR  : "


class peprCandidateFromHitProducer: public edm::stream::EDProducer<> {
 public:
    explicit peprCandidateFromHitProducer(const edm::ParameterSet&);
    ~peprCandidateFromHitProducer();

 private:
    void beginStream(edm::StreamID);
    void endStream();
    void produce(edm::Event&, const edm::EventSetup&) override;

    void fillWindows(const edm::Event&);
    void writeInputArrays(const std::vector<std::vector<float> >&);
    void readOutputArrays();

    // options
    std::vector<edm::InputTag> recHitCollections_;

    // tokens
    std::vector<edm::EDGetTokenT<HGCRecHitCollection> > recHitTokens_;  
    edm::EDGetTokenT<edm::View<reco::Track>> tracksToken_;
    edm::EDGetTokenT<std::vector<SimCluster>> simClusterToken_;

    //FIXME, to be replaced
    bool batchedModel_;
    std::string tritonPath_;
    std::string inpipeName_;
    std::string outpipeName_;

    // rechit tools
    hgcal::RecHitTools recHitTools_;

    // windows
    std::vector<InferenceWindow> windows_;

    double minEta_;
    double maxEta_;
    double etaFrameWidth_;
    double phiFrameWidth_;
    size_t nEtaSegments_;
    size_t nPhiSegments_;


};


peprCandidateFromHitProducer::peprCandidateFromHitProducer(const edm::ParameterSet& config) :
        recHitCollections_(config.getParameter<std::vector<edm::InputTag> >("recHitCollections")), 
        tracksToken_(consumes<edm::View<reco::Track>>(config.getParameter<edm::InputTag>("tracks"))),
        simClusterToken_(consumes<std::vector<SimCluster>>(config.getParameter<edm::InputTag>("simClusters"))), 
        batchedModel_(config.getParameter<bool>("batchedModel")), 
        tritonPath_(config.getParameter<std::string>("tritonPath")), 
        inpipeName_(config.getParameter<std::string>("inpipeName")),
        outpipeName_(config.getParameter<std::string>("outpipeName")),
        //FIXME: actually these are all not needed if windows are created in the constructor!
        minEta_(config.getParameter<double>("minEta")),
        maxEta_(config.getParameter<double>("maxEta")),
        etaFrameWidth_(config.getParameter<double>("etaFrameWidth")),
        phiFrameWidth_(config.getParameter<double>("phiFrameWidth")),
        nEtaSegments_((size_t)config.getParameter<uint32_t>("nEtaSegments")),
        nPhiSegments_((size_t)config.getParameter<uint32_t>("nPhiSegments"))
        {

    // sanity checks for sliding windows

    // get tokens
    for (edm::InputTag& recHitCollection : recHitCollections_) {
        recHitTokens_.push_back(consumes<HGCRecHitCollection>(recHitCollection));
    }
   

    produces<reco::PFCandidateCollection>();
    //produces<std::vector<Trackster>>();


    //launch scripts to work with Triton clients
    //https://github.com/cms-pepr/HGCalML/tree/master/triton
    
    std::string str0 = "cd " + tritonPath_;
    const char *cdcommand = str0.c_str();
    system(cdcommand);

    std::string str1 = "./cmssw_oc_server.sh > serverlog.txt &";
    const char *clientcommand1 = str1.c_str();
    std::cout << "Executing: " << clientcommand1 << std::endl;
    system(clientcommand1); 

    std::cout << "Sleeping 30s... " << std::endl;
    system("sleep 30");

    std::string str2 = "./cmssw_oc_forward_client.sh " + inpipeName_;
    const char *clientcommand2 = str2.c_str();
    std::cout << "Executing: " << clientcommand2 << std::endl;
    system(clientcommand2);

    system("cd -");

}

peprCandidateFromHitProducer::~peprCandidateFromHitProducer() {

    //FIXME: find process ID of scripts run in the constructor
    //system("killtask <...> 15");
    //system("killtask <...> 15");
}


void peprCandidateFromHitProducer::beginStream(edm::StreamID streamId) {
    windows_ = InferenceWindow::createWindows(nPhiSegments_,nEtaSegments_,minEta_,maxEta_,etaFrameWidth_,phiFrameWidth_);

    // FIXME, make configurable?
    for(auto& w: windows_)
        w.setMode(WindowBase::useRechits);
}

void peprCandidateFromHitProducer::endStream() {

    //windows_.clear();
}

void peprCandidateFromHitProducer::produce(edm::Event& event, const edm::EventSetup& setup) {

    recHitTools_.getEventSetup(setup);

    auto simclusters = event.get(simClusterToken_);
    std::cout << "The size of the simclusters is " << simclusters.size() << std::endl;


    // fill rechits into windows
    fillWindows(event);


    //FIXME: check window treatment
    //std::vector<tensorflow::Tensor> windowoutputs;
    // run the evaluation per window
    std::vector<std::vector<float> >  hitFeatures;
    for (auto & window : windows_) {

        hitFeatures = window.getHitFeatures();
        std::cout << "  hitFeatures size = " << hitFeatures.size() << std::endl;
        writeInputArrays(hitFeatures);

        //FIXME: how to check if output exists? (How long is waiting time? Rest of code on hold until output pipe available?)
        readOutputArrays();        

        //windowoutputs.push_back(window.getOutput());
    }



    // reconstruct showers using all windows and put them into the event
    //reconstructShowers();


    //FIXME: block below temporarily uses simclusters as a placeholder for OC candidates
    //FIXME: check window treatment

    // making candidate collection
    auto candidates = std::make_unique<reco::PFCandidateCollection>();
    std::cout << "Creating PF candidates " << std::endl;
    //for (unsigned int i=0; i<windows_.size(); i++) {
        //loop over windows
        //std::cout << "Window " << i << std::endl;

        
        // FIXME: need to take care of converting E, px, py, pz <--> XYZT whenever needed
        //float X = -9999., Y=-9999., Z=-9999., E=-9999.;

        // loop over particles reconstructed in the window
        //FIXME: Simclusters size to be revisited if working with windows
        for(size_t it=0;it<simclusters.size();it++) {

            //std::cout << " particle " << it << std::endl;

            //const auto abs_pdg_id = -9999;
            const auto charge = 0; // FIXME!
                       
            //block inspired by calcP4 method in TICL TracksterP4FromEnergySum plugin
            //starts from 'position (X,Y,Z)'
            //math::XYZVector direction(X, Y, Z);
            math::XYZVectorF direction = simclusters.at(it).momentum();
            direction = direction.Unit();
            math::XYZTLorentzVector cartesian(direction.X(), direction.Y(), direction.Z(), simclusters.at(it).energy());
            //// Convert px, py, pz, E vector to CMS standard pt/eta/phi/m vector
            reco::Candidate::LorentzVector p4(cartesian);

            //const auto& four_mom = math::XYZTLorentzVector(X,Y,Z,E);
            const auto& four_mom = p4;
            std::cout << "    PF particle energy " << four_mom.E() << ", Px = " << four_mom.Px() << ", Py = " << four_mom.Py() << ", Pz = " << four_mom.Pz() << std::endl;
            reco::PFCandidate::ParticleType part_type = reco::PFCandidate::X;
            candidates->emplace_back(charge, four_mom, part_type);

        }
    //}



    event.put(std::move(candidates));

    std::cout << "Results produced and put in event" << std::endl;


    // clear all windows
    for (auto& window : windows_) {
        window.clear();
    }

    std::cout << "Windows cleared" << std::endl;
}



void peprCandidateFromHitProducer::fillWindows(const edm::Event& event) {

    if (!windows_.size()) {
        throw cms::Exception("NoWindows") << "no windows initialized";
    }
    
    std::cout << "Number of windows = " << windows_.size() << std::endl;
    
    // copied block from window ntupler code
    // get rechits, get positions and merge collections
    std::vector<HGCRecHitWithPos> allrechits;
    for (auto & token : recHitTokens_) {
       for (const auto& rh : event.get(token)) {
           HGCRecHitWithPos rhp = { const_cast<HGCRecHit*>(&rh), recHitTools_.getPosition(rh.detid()) };
           allrechits.push_back(rhp);
       }
    }
    // sort according to the energy
    std::sort(allrechits.begin(), allrechits.end(), 
        [](const HGCRecHitWithPos& rh1, const HGCRecHitWithPos& rh2) 
            { return rh1.hit->energy() > rh2.hit->energy();});

    
    // fills a vector of the specified size with zeroes (entries will be 0 if rechit is not filled, and 1 if it is filled)
    std::vector<size_t> filledrechits(allrechits.size(),0);


    for (auto & window : windows_) {
        //fill rechits in this window
        for(size_t it=0;it<allrechits.size();it++) {
           if(filledrechits.at(it)>3) continue;
           if(window.maybeAddRecHit(allrechits.at(it)))
               filledrechits.at(it)++;
        }

        window.fillFeatureArrays();
    }


}


void peprCandidateFromHitProducer::writeInputArrays(const std::vector<std::vector<float> >& hitFeatures) {

    std::string inpipeRAM = "/dev/shm/" + inpipeName_; //write in RAM
    //std::ofstream inputArrayStream(inpipeName_.c_str()); //do not write in RAM (local testing)
    std::ofstream inputArrayStream(inpipeRAM.c_str()); 

    for (size_t i=0; i<hitFeatures.size(); i++) {

        for (size_t j=0; j<hitFeatures[i].size(); j++) {
            inputArrayStream << hitFeatures[i][j] << " ";
        }
        inputArrayStream << "\n";
    }

    inputArrayStream.close();
}


void peprCandidateFromHitProducer::readOutputArrays() {

    //FIXME: currently just reading and printing, proper values to be stored
    //FIXME: read until you get a last empty line, because the fifo is being filled not instantly

    std::string outpipeRAM = "/dev/shm/" + outpipeName_; //read from RAM
    //std::ifstream outputArrayStream(outpipeName_.c_str()); //do not read from RAM (local testing)
    std::ifstream outputArrayStream(outpipeRAM.c_str()); 

    std::string line;
    if(outputArrayStream.is_open()) {
        std::cout << "Output array file opened!" << std::endl;
        while ( getline (outputArrayStream,line) )
        {
            std::cout << line << '\n';
        }
        outputArrayStream.close();
    }


}

//remove

DEFINE_FWK_MODULE(peprCandidateFromHitProducer);
