# Running the pepr PF candidate producer for HGCAL

This example demonstrates how to run the particle reconstruction in the HGCAL subdetector via inference of graph neural networks. 

## Setup

Install CMSSW:
```
export SCRAM_ARCH="slc7_amd64_gcc820"
cmsrel CMSSW_11_1_0_pre7
cd CMSSW_11_1_0_pre7/src
cmsenv
scram b -j 8
```

Install custom packages: 
```
git cms-init
git cms-merge-topic gvonsem:pepr_CMSSW_11_1_0_pre7_peprCandDev
git clone --recursive https://github.com/CMS-HGCAL/reco-prodtools.git reco_prodtools
scram b -j 8
```

Create prodtools templates:
```
cd reco_prodtools/templates/python
./produceSkeletons_D49_NoSmear_NoDQMNoHLT_PU_AVE_200_BX_25ns.sh
cd ../../../
scram b python
```

Check out the configuration files to generate events:
```
cd ../../
git clone -b pepr_CMSSW_11_1_0_pre7_peprCandDev https://github.com/gvonsem/production_tests.git
```

## Generate events

First we produce GEN-SIM-DIGI (GSD) events, in this example by shooting particles (e.g. photons) 
in a certain energy range towards the HGCAL subdetector via the `FlatEtaRangeGunProducer`.
```
cd production_tests
cmsRun GSD_GUN.py seed=1 outputFile="file:1_GSD.root" maxEvents=5
```
Once the GSD events are produced, we can run the reconstruction step: 
```
cmsRun RECO.py inputFiles="file://1_GSD.root" outputFile="file:1_RECO.root" outputFileDQM="file:1_DQM.root" maxEvents=5
```
A dedicated EDProducer module, the `peprCandidateFromHitProducer` located 
in the [RecoHGCAL/GraphReco](https://github.com/gvonsem/cmssw/tree/pepr_CMSSW_11_1_0_pre7_peprCandDev/RecoHGCal/GraphReco) package, 
produces PF candidates straight from rechit information, in this example via the [Object Condensation](https://arxiv.org/abs/2002.03605v3) method. 
The inference of trained graph neural network models is done by sending the rechit information per endcap to a custom Triton server, evaluating the model, 
and retrieving the regressed energy and position of clustered particle candidates. 
These candidates are subsequently turned into a PFcandidate collection named `recoPFCandidates_peprCandidateFromHitProducer__RECO`. Particle and charge identification as well as track-cluster matching are work in progress and not included yet. 

**Note:** it may take some time for the event loop in the reconstruction to start, and inference may be slow on a CPU server. The speed of communication with the client and especially the inference will improve drastically once dedicated Triton GPU servers are used. 

The sequence of the producer module is as follows:
* In the constructor of the producer, the Triton client is started.
* The producer checks for open pipes (set up to communicate with the Triton server) and will wait until the pipes are open to send the rechit data to the server.
* In case the pipes are open before the producer reaches the check, the client will wait until the rechit data is passed from the producer.
* The inference itself is done on the Triton server via the trained model that is stored there, and the results are passed back to the module where a collection of reconstructed particle candidates is created.
* The Triton client is automatically closed in the destructor of the producer


