// -*- C++ -*-
//
// Package:    EcalValidation
// Class:      EcalValidation
// Original Author:  Martina Malberti
//
// system include files
#include <memory>
#include "PhysicsTools/Utilities/interface/LumiReWeighting.h"


// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Common/interface/EventBase.h"
#include "DataFormats/Provenance/interface/EventAuxiliary.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "Geometry/CaloTopology/interface/CaloTopology.h"
#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "Geometry/CaloGeometry/interface/CaloSubdetectorGeometry.h"
#include "Geometry/CaloEventSetup/interface/CaloTopologyRecord.h"
#include "Geometry/CaloTopology/interface/CaloSubdetectorTopology.h"
#include "DataFormats/GeometryVector/interface/GlobalPoint.h"
#include "Geometry/EcalAlgo/interface/EcalBarrelGeometry.h"
#include "Geometry/EcalAlgo/interface/EcalEndcapGeometry.h"

#include "DataFormats/EcalDetId/interface/EBDetId.h"
#include "DataFormats/EcalDetId/interface/EEDetId.h"
#include "DataFormats/EcalDetId/interface/ESDetId.h"
#include "DataFormats/EcalDigi/interface/EcalDigiCollections.h"
#include "DataFormats/EcalRecHit/interface/EcalRecHitCollections.h"
#include "DataFormats/ParticleFlowReco/interface/PFRecHit.h"
#include "DataFormats/ParticleFlowReco/interface/PFRecHitFwd.h"
#include "DataFormats/ParticleFlowReco/interface/PFRecHitFraction.h"
#include "DataFormats/EgammaReco/interface/SuperCluster.h"
#include "DataFormats/EgammaReco/interface/SuperClusterFwd.h"
#include "DataFormats/EgammaReco/interface/BasicCluster.h"
#include "DataFormats/EgammaReco/interface/BasicClusterFwd.h"
#include "DataFormats/EgammaReco/interface/PreshowerCluster.h"
#include "DataFormats/EgammaReco/interface/PreshowerClusterFwd.h"
#include "RecoEcal/EgammaCoreTools/interface/EcalClusterTools.h"
#include "RecoEcal/EgammaCoreTools/interface/EcalTools.h"
#include "RecoLocalCalo/EcalRecAlgos/interface/EcalSeverityLevelAlgo.h"
#include "RecoLocalCalo/EcalRecAlgos/interface/EcalCleaningAlgo.h"
#include "RecoEcal/EgammaCoreTools/interface/EcalRecHitLess.h"
#include "CalibCalorimetry/EcalLaserCorrection/interface/EcalLaserDbService.h"
#include "CalibCalorimetry/EcalLaserCorrection/interface/EcalLaserDbRecord.h"

#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"

#include "RecoEgamma/EgammaTools/interface/ECALPositionCalculator.h"

#include "DataFormats/JetReco/interface/CaloJet.h"
#include "DataFormats/JetReco/interface/CaloJetCollection.h"
#include "CondFormats/EcalObjects/interface/EcalPedestals.h"
#include "CondFormats/DataRecord/interface/EcalPedestalsRcd.h"

#include "Validation/EcalValidation/interface/EcalValidationAOD.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"
#include "MagneticField/Engine/interface/MagneticField.h"

#include "TVector3.h"

#include <iostream>
#include <cmath>
#include <fstream>
#include <map>
#include <iterator>

using namespace cms ;
using namespace edm ;
using namespace std ;
using namespace reco;
bool debug(false);
double LumiWeights[200]= {0,	0,	0,	0,	0,	0,	2143.52,	3346.06,	3278.98,	764.818,	282.99,	154.118,	96.0512,58.7164,	34.0694,	24.2334,	16.8311,	11.6877,	8.31968,	6.2278,	4.74167,	3.7003,	2.9184,	2.32909,1.87553,	1.57232,	1.2856,	1.05083,	0.886188,	0.718918,	0.597554,	0.492347,	0.399543,	0.320542,	0.260795,	0.205339,	0.162223,	0.126568,	0.0989134,	0.0773308,	0.0606376,	0.0469714,	0.0364722,	0.0286898,	0.0224534,	0.0179116,	0.0142907,	0.0116515,	0.00949087,	0.00783991,	0.00649521,	0.00533664,	0.00468546,	0.0040227,	0.00371085,	0.00282875,	0.0027305,	0.00246812,	0.0022999,	0.00222884,	0.00217703,	0.00265132,	0.0019689,	0.00263327,	0.00154821,	0.00182158,	0.00349053,	0.00359745,	0.0032268,	0.00213538,	0.00302513,	0.0024201,	0.048402,	0.008067,	0,	0.008067,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	};



//{0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0.117733,	0.138446,	0.535993,	1.05875,	1.04964,	2.02185,	2.64248,	3.48501,	4.75895,	4.50494,	5.62076,	6.53113,	7.77193,	8.04238,	7.71062,	8.71349,	7.91055,	7.08059,	7.26529,	6.22432,	5.39497,	5.09417,	3.99162,	3.79948,	2.72379,	2.54813,	2.03363,	1.57913,	1.5133,	0.99523,	0.698623,	0.763208,	0.458599,	0.420863,	0.308687,	0.269633,	0.141078,	0.163361,	0.0378794,	0.0559184,	0.036519,	0.0266977,	0.0304644,	0.0495151,	0.0523857,	0.0206451,	0.016114,	0.00769526,	0.00294584,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	};
//{0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	121.281,	235.739,	20.9313,	17.8651,	14.0354,	7.12529,	4.786,	3.70699,	2.43642,	1.9902,	1.65699,	1.27054,	1.01939,	0.990199,	0.947373,	0.83607,	0.858608,	0.883404,	0.767167,	0.8434,	0.898385,	0.838823,	0.929265,	0.868659,	0.961145,	1.02493,	1.03648,	1.03746,	0.940245,	1.16427,	1.33602,	1.11229,	1.3686,	1.52682,	1.38993,	2.15667,	2.24537,	1.9635,	2.31615,	2.89896,	3.00573,	3.9015,	2.44967,	3.02912,	4.37512,	5.81201,	9.73745,	21.6388,	32.9734,	66.6246,	72.987,	0,	0,	0,	3375.65,	0,	2025.39,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	};
//
//
// constructors and destructor
//
int ttt(1);
EcalValidationAOD::EcalValidationAOD(const edm::ParameterSet& ps)
{
    cout << "in constructor" << endl;
//LumiReWeighting::LumiReWeighting( std::string generatedFile,
//                                                            std::string dataFile,
//                                                            std::string GenHistName1 = "MCpileup",
//                                                            std::string DataHistName2 = "datapileup" )
//LumiWeights_ = edm::LumiReWeighting("MCInfo.root", "MyDataPileupHistogram_2017.root", "pileup", "rat");
    //now do what ever initialization is needed
    puCollection_              = consumes<vector<PileupSummaryInfo> >    (ps.getParameter<InputTag>("pileupCollection"));
    PV_                        = consumes<reco::VertexCollection>(ps.getParameter<edm::InputTag>("PVTag"));
    //recHitCollection_EB_       = consumes<EcalRecHitCollection>(ps.getParameter<edm::InputTag>("recHitCollection_EB"));
    //recHitCollection_EE_       = consumes<EcalRecHitCollection>(ps.getParameter<edm::InputTag>("recHitCollection_EE"));
    recHitCollection_EB_    = consumes<EcalRecHitCollection>(ps.getParameter<edm::InputTag>("recHitCollection_EB"));
    recHitCollection_EE_    = consumes<EcalRecHitCollection>(ps.getParameter<edm::InputTag>("recHitCollection_EE"));
//    PFrecHitCollection_EB_    = consumes<reco::PFRecHitCollection>(ps.getParameter<edm::InputTag>("PFrecHitCollection_EB"));
//    PFrecHitCollection_EE_    = consumes<reco::PFRecHitCollection>(ps.getParameter<edm::InputTag>("PFrecHitCollection_EE"));
    basicClusterCollection_EB_ = consumes<reco::BasicClusterCollection>(ps.getParameter<edm::InputTag>("basicClusterCollection_EB"));
    basicClusterCollection_EE_ = consumes<reco::BasicClusterCollection>(ps.getParameter<edm::InputTag>("basicClusterCollection_EE"));
    superClusterCollection_EB_ = consumes<reco::SuperClusterCollection>(ps.getParameter<edm::InputTag>("superClusterCollection_EB"));
    superClusterCollection_EE_ = consumes<reco::SuperClusterCollection>(ps.getParameter<edm::InputTag>("superClusterCollection_EE"));
    esClusterCollectionX_      = consumes<PreshowerClusterCollection>(ps.getParameter<edm::InputTag>("ClusterCollectionX_ES"));
    esClusterCollectionY_      = consumes<PreshowerClusterCollection>(ps.getParameter<edm::InputTag>("ClusterCollectionY_ES"));

    tracks_                    = consumes<edm::View<reco::Track> >(ps.getParameter<edm::InputTag>("tracks"));
    beamSpot_                  = consumes<reco::BeamSpot>(ps.getParameter<edm::InputTag>("beamSpot"));
    jets_                      = consumes<reco::CaloJetCollection>(ps.getParameter<edm::InputTag>("jets"));

    ethrEB_                    = ps.getParameter<double>("ethrEB");
    isMC                       = ps.getParameter<bool>("isMC");
    ethrEE_                    = ps.getParameter<double>("ethrEE");
    gainId_                    = ps.getParameter<double>("gainId");

    scEtThrEB_                 = ps.getParameter<double>("scEtThrEB");
    scEtThrEE_                 = ps.getParameter<double>("scEtThrEE");

    ///for Pi0 barrel selection
    isMonEBpi0_ = ps.getUntrackedParameter<bool>("isMonEBpi0",false);

    selePtGamma_ = ps.getParameter<double> ("selePtGamma");
    selePtPi0_ = ps.getParameter<double> ("selePtPi0");
    seleMinvMaxPi0_ = ps.getParameter<double> ("seleMinvMaxPi0");
    seleMinvMinPi0_ = ps.getParameter<double> ("seleMinvMinPi0");
    seleS4S9Gamma_ = ps.getParameter<double> ("seleS4S9Gamma");
    selePi0Iso_ = ps.getParameter<double> ("selePi0Iso");
    ptMinForIsolation_ = ps.getParameter<double> ("ptMinForIsolation");
    selePi0BeltDR_ = ps.getParameter<double> ("selePi0BeltDR");
    selePi0BeltDeta_ = ps.getParameter<double> ("selePi0BeltDeta");


    ///for Pi0 endcap selection
    isMonEEpi0_ = ps.getUntrackedParameter<bool>("isMonEEpi0",false);

    selePtGamma_EE_ = ps.getParameter<double> ("selePtGamma_EE");
    selePtPi0_EE_ = ps.getParameter<double> ("selePtPi0_EE");
    seleS4S9Gamma_EE_ = ps.getParameter<double> ("seleS4S9Gamma_EE");
    seleMinvMaxPi0_EE_ = ps.getParameter<double> ("seleMinvMaxPi0_EE");
    seleMinvMinPi0_EE_ = ps.getParameter<double> ("seleMinvMinPi0_EE");
    ptMinForIsolation_EE_ = ps.getParameter<double> ("ptMinForIsolation_EE");
    selePi0BeltDR_EE_ = ps.getParameter<double> ("selePi0BeltDR_EE");
    selePi0BeltDeta_EE_ = ps.getParameter<double> ("selePi0BeltDeta_EE");
    selePi0Iso_EE_ = ps.getParameter<double> ("selePi0Iso_EE");

    ///Endcap regions for pi0:
    ///try to divide endcap region into 3 parts
    /// eta< 2 ; eta>2 && eta<2.5 ; eta>2.5;
    region1_Pi0_EE_ = ps.getParameter<double> ("region1_Pi0_EE");
    selePtGammaPi0_EE_region1_ = ps.getParameter<double> ("selePtGammaPi0_EE_region1");
    selePtPi0_EE_region1_ = ps.getParameter<double> ("selePtPi0_EE_region1");

    region2_Pi0_EE_ = ps.getParameter<double> ("region2_Pi0_EE");
    selePtGammaPi0_EE_region2_ = ps.getParameter<double> ("selePtGammaPi0_EE_region2");
    selePtPi0_EE_region2_ = ps.getParameter<double> ("selePtPi0_EE_region2");

    selePtGammaPi0_EE_region3_ = ps.getParameter<double> ("selePtGammaPi0_EE_region3");
    selePtPi0_EE_region3_ = ps.getParameter<double> ("selePtPi0_EE_region3");

    isMaskEB_ = ps.getUntrackedParameter<bool>("isMaskEB",false);
    isMaskEE_ = ps.getUntrackedParameter<bool>("isMaskEE",false);

    maskEBFile_=  ps.getUntrackedParameter<string>("maskEBFile","maskEB.txt");
    maskEEFile_=  ps.getUntrackedParameter<string>("maskEEFile","maskEE.txt");

    useRecoFlag_ = ps.getUntrackedParameter<bool>("useRecoFlag",false);
    usePFRecHitFlag_ = ps.getParameter<bool>("usePFRecHitFlag");
    
    clusSeedThr_ = ps.getParameter<double> ("clusSeedThr");
    clusSeedThr_EE_ = ps.getParameter<double> ("clusSeedThr_EE");
    clusEtaSize_ = ps.getParameter<int> ("clusEtaSize");
    clusPhiSize_ = ps.getParameter<int> ("clusPhiSize");

    seleXtalMinEnergy_ = ps.getParameter<double>("seleXtalMinEnergy");
    seleXtalMinEnergy_EE_ = ps.getParameter<double>("seleXtalMinEnergy_EE");

    edm::ParameterSet posCalcParameters = ps.getParameter<edm::ParameterSet>("posCalcParameters");
    posCalculator_ = PositionCalc(posCalcParameters);

    naiveId_ = 0;


    // histos

    edm::Service<TFileService> fs;

    //   TFileDirectory dirEB = fs->mkdir("EB");
    //   TFileDirectory dirEE = fs->mkdir("EE");
    //   TFileDirectory dirES = fs->mkdir("ES");

    TFileDirectory dirPi0 = fs->mkdir("Pi0");
    tree = new TTree("test","");
    h_numberOfEvents = fs->make<TH1D>("h_numberOfEvents","h_numberOfEvents",10,0,10);
    tree->Branch("nPU_",  &nPU_);
    tree->Branch("puTrue_", &puTrue_);
    tree->Branch("puBX_", &puBX_);

    h_PV_n = fs->make<TH1D>("h_PV_n","h_PV_n",150,0.,150.);
    h_PV_nw = fs->make<TH1D>("h_PV_nw","h_PV_nw",150,0.,150.);
    hPU_        = fs->make<TH1F>("hPU",        "number of pileup",      200,  0, 200);
    hPUTrue_    = fs->make<TH1F>("hPUTrue",    "number of true pilepu", 1000, 0, 200);


    // RecHits ----------------------------------------------

    // ... DataBase noise map
    h_DB_noiseMap_EB           = fs->make<TH2D>("h_DB_noiseMap_EB","h_DB_noiseMap_EB",360,1.,361.,172,-86.,86. );
    h_DB_noiseMap_EEP          = fs->make<TH2D>("h_DB_noiseMap_EEP","h_DB_noiseMap_EEP",100,0.,100.,100,0.,100. );
    h_DB_noiseMap_EEM          = fs->make<TH2D>("h_DB_noiseMap_EEM","h_DB_noiseMap_EEM",100,0.,100.,100,0.,100. );


    // ... DataBase LaserCorr map
    h_DB_LaserCorrMap_EB           = fs->make<TH2D>("h_DB_LaserCorrMap_EB","h_DB_LaserCorrMap_EB",360,1.,361.,172,-86.,86. );
    h_DB_LaserCorrMap_EEP          = fs->make<TH2D>("h_DB_LaserCorrMap_EEP","h_DB_LaserCorrMap_EEP",100,0.,100.,100,0.,100. );
    h_DB_LaserCorrMap_EEM          = fs->make<TH2D>("h_DB_LaserCorrMap_EEM","h_DB_LaserCorrMap_EEM",100,0.,100.,100,0.,100. );

    // ... barrel
    h_recHits_EB_size          = fs->make<TH1D>("h_recHits_EB_size", "h_recHitsEB_size", 100000, 0, 100000 );
    h_recHits_EB_energy        = fs->make<TH1D>("h_recHits_EB_energy","h_recHitsEB_energy",11000,-50,500);
    h_recHits_EB_energy_thr    = fs->make<TH1D>("h_recHits_EB_energy_thr","h_recHitsEB_energy_thr",2000,-50,500);
    h_recHits_EB_energyMax     = fs->make<TH1D>("h_recHits_EB_energyMax","h_recHitsEB_energyMax",2000,-50,500);
    h_recHits_EB_sumEt         = fs->make<TH1D>("h_recHits_EB_sumEt","h_recHits_EB_sumEt",40000,0,20000);
    h_recHits_EB_kOutOfTime          = fs->make<TH1D>("h_recHits_EB_kOutOfTime","h_recHits_EB_kOutOfTime",2, 0, 2);
    h_recHits_EB_kWeird          = fs->make<TH1D>("h_recHits_EB_kWeird","h_recHits_EB_kWeird",2, 0, 2);
    h_recHits_EB_time_cut          = fs->make<TH1D>("h_recHits_EB_time_cut","h_recHits_EB_time_cut",400,-100,100);
    h_recHits_EB_Chi2_cut          = fs->make<TH1D>("h_recHits_EB_Chi2_cut","h_recHits_EB_Chi2_cut",1000,0,100);
    h_recHits_EB_time          = fs->make<TH1D>("h_recHits_EB_time","h_recHits_EB_time",400,-100,100);
    h_recHits_EB_Chi2          = fs->make<TH1D>("h_recHits_EB_Chi2","h_recHits_EB_Chi2",1000,0,100);
    h_recHits_EB_OutOfTimeChi2 = fs->make<TH1D>("h_recHits_EB_OutOfTimeChi2","h_recHits_EB_OutOfTimeChi2",1000,0,100);
    h_recHits_EB_E1oE4         = fs->make<TH1D>("h_recHits_EB_E1oE4","h_recHitsEB_E1oE4",150, 0, 1.5);
    h_recHits_EB_iPhiOccupancy = fs->make<TH1D>("h_recHits_EB_iPhiOccupancy","h_recHits_EB_iPhiOccupancy",360,1.,361. );
    h_recHits_EB_iEtaOccupancy = fs->make<TH1D>("h_recHits_EB_iEtaOccupancy","h_recHits_EB_iEtaOccupancy",172,-86.,86.);
    h_recHits_EB_occupancy     = fs->make<TH2D>("h_recHits_EB_occupancy","h_recHits_EB_occupancy",360,1.,361.,172,-86.,86. );
    h_recHits_EB_deviation     = fs->make<TH2D>("h_recHits_EB_deviation","h_recHits_EB_deviation",360,1.,361.,172,-86.,86. );


    h_recHits_EB_energy_spike  = fs->make<TH1D>("h_recHits_EB_energy_spike","h_recHitsEB_energy_spike",2000,0,500);

    // ... barrel (with spike cleaning)
    h_recHits_EB_size_cleaned          = fs->make<TH1D>("h_recHits_EB_size_cleaned", "h_recHitsEB_size_cleaned", 1000, 0, 10000 );
    h_recHits_EB_energy_cleaned        = fs->make<TH1D>("h_recHits_EB_energy_cleaned","h_recHitsEB_energy_cleaned",11000,-50,500);
    h_recHits_EB_energyMax_cleaned     = fs->make<TH1D>("h_recHits_EB_energyMax_cleaned","h_recHitsEB_energyMax_cleaned",2000,-50,500);
    h_recHits_EB_time_cleaned          = fs->make<TH1D>("h_recHits_EB_time_cleaned","h_recHits_EB_time_cleaned",400,-100,100);
    h_recHits_EB_Chi2_cleaned          = fs->make<TH1D>("h_recHits_EB_Chi2_cleaned","h_recHits_EB_Chi2_cleaned",1000,0,100);
    h_recHits_EB_OutOfTimeChi2_cleaned = fs->make<TH1D>("h_recHits_EB_OutOfTimeChi2_cleaned","h_recHits_EB_OutOfTimeChi2_cleaned",1000,0,100);
    h_recHits_EB_recoFlag = fs->make<TH1D>("h_recHits_EB_recoFlag","h_recHits_EB_recoFlag",16,-0.5,15.5);

    // ... endcap
    h_recHits_EE_size           = fs->make<TH1D>("h_recHits_EE_size","h_recHits_EE_size",10000,0,10000);
    h_recHits_EE_recoFlag = fs->make<TH1D>("h_recHits_EE_recoFlag","h_recHits_EE_recoFlag",16,-0.5,15.5);

    h_recHits_EEP_size          = fs->make<TH1D>("h_recHits_EEP_size","h_recHits_EEP_size",1000,0,10000);
    h_recHits_EEP_energy        = fs->make<TH1D>("h_recHits_EEP_energy","h_recHits_EEP_energy",11000,-50,500);
    h_recHits_EEP_energyMax     = fs->make<TH1D>("h_recHits_EEP_energyMax","h_recHitsEEP_energyMax",2000,-50,500);
    h_recHits_EEP_sumEt         = fs->make<TH1D>("h_recHits_EEP_sumEt","h_recHits_EEP_sumEt",2000,0,1000);
    h_recHits_EEP_sumEtCut         = fs->make<TH1D>("h_recHits_EEP_sumEtCut","h_recHits_EEP_sumEtCut",2000,0,1000);
    h_recHits_EEP_time          = fs->make<TH1D>("h_recHits_EEP_time","h_recHits_EEP_time",400,-100,100);
    h_recHits_EEP_Chi2          = fs->make<TH1D>("h_recHits_EEP_Chi2","h_recHits_EEP_Chi2",1000,0,100);
    h_recHits_EEP_OutOfTimeChi2 = fs->make<TH1D>("h_recHits_EEP_OutOfTimeChi2","h_recHits_EEP_OutOfTimeChi2",1000,0,100);
    h_recHits_EEP_E1oE4         = fs->make<TH1D>("h_recHits_EEP_E1oE4","h_recHitsEEP_E1oE4",150, 0, 1.5);
    h_recHits_EEP_iXoccupancy   = fs->make<TH1D>("h_recHits_EEP_iXoccupancy","h_recHits_EEP_iXoccupancy",100,0.,100.);
    h_recHits_EEP_iYoccupancy   = fs->make<TH1D>("h_recHits_EEP_iYoccupancy","h_recHits_EEP_iYoccupancy",100,0.,100.);
    h_recHits_EEP_occupancy     = fs->make<TH2D>("h_recHits_EEP_occupancy","h_recHits_EEP_occupancy",100,0.,100.,100,0.,100. );
    h_recHits_EEP_deviation     = fs->make<TH2D>("h_recHits_EEP_deviation","h_recHits_EEP_deviation",100,0.,100.,100,0.,100. );

    h_recHits_EEM_size          = fs->make<TH1D>("h_recHits_EEM_size","h_recHits_EEM_size",1000,0,10000);
    h_recHits_EEM_energy        = fs->make<TH1D>("h_recHits_EEM_energy","h_recHits_EEM_energy",11000,-50,500);
    h_recHits_EEM_energyMax     = fs->make<TH1D>("h_recHits_EEM_energyMax","h_recHits_EEM_energyMax",2000,-50,500);
    h_recHits_EEM_sumEt         = fs->make<TH1D>("h_recHits_EEM_sumEt","h_recHits_EEM_sumEt",2000,0,1000);
    h_recHits_EEM_sumEtCut      = fs->make<TH1D>("h_recHits_EEM_sumEtCut","h_recHits_EEM_sumEtCut",2000,0,1000);
    h_recHits_EEM_time          = fs->make<TH1D>("h_recHits_EEM_time","h_recHits_EEM_time",400,-100,100);
    h_recHits_EEM_Chi2          = fs->make<TH1D>("h_recHits_EEM_Chi2","h_recHits_EEM_Chi2",1000,0,100);
    h_recHits_EEM_OutOfTimeChi2 = fs->make<TH1D>("h_recHits_EEM_OutOfTimeChi2","h_recHits_EEM_OutOfTimeChi2",1000,0,100);
    h_recHits_EEM_E1oE4         = fs->make<TH1D>("h_recHits_EEM_E1oE4","h_recHitsEEM_E1oE4",150, 0, 1.5);
    h_recHits_EEM_iXoccupancy   = fs->make<TH1D>("h_recHits_EEM_iXoccupancy","h_recHits_EEM_iXoccupancy",100,0.,100.);
    h_recHits_EEM_iYoccupancy   = fs->make<TH1D>("h_recHits_EEM_iYoccupancy","h_recHits_EEM_iYoccupancy",100,0.,100.);
    h_recHits_EEM_occupancy     = fs->make<TH2D>("h_recHits_EEM_occupancy","h_recHits_EEM_occupancy",100,0.,100.,100,0.,100. );
    h_recHits_EEM_deviation     = fs->make<TH2D>("h_recHits_EEM_deviation","h_recHits_EEM_deviation",100,0.,100.,100,0.,100. );

    // eta/phi distributions
    h_recHits_eta          = fs->make<TH1D>("h_recHits_eta","h_recHits_eta",300,-3.,3.);
    h_recHits_EB_eta       = fs->make<TH1D>("h_recHits_EB_eta","h_recHits_EB_eta",170,-85*0.0175,85*0.0175);
    h_recHits_EB_Eta       = fs->make<TH1D>("h_recHits_EB_Eta","h_recHits_EB_Eta",170,-85*0.0175,85*0.0175);
    h_recHits_EB_Phi       = fs->make<TH1D>("h_recHits_EB_Phi","h_recHits_EB_Phi",360,-3.1415927, 3.1415927);

    h_recHits_EEP_eta      = fs->make<TH1D>("h_recHits_EEP_eta","h_recHits_EEP_eta",300,-3.,3.);
    h_recHits_EEM_eta      = fs->make<TH1D>("h_recHits_EEM_eta","h_recHits_EEM_eta",300,-3.,3.);

    h_recHits_EB_phi       = fs->make<TH1D>("h_recHits_EB_phi","h_recHits_EB_phi",360,-3.1415927, 3.1415927);
    h_recHits_EE_phi       = fs->make<TH1D>("h_recHits_EE_phi","h_recHits_EE_phi",360,-3.1415927, 3.1415927);
    h_recHits_EEP_phi      = fs->make<TH1D>("h_recHits_EEP_phi","h_recHits_EEP_phi",360,-3.1415927, 3.1415927);
    h_recHits_EEM_phi      = fs->make<TH1D>("h_recHits_EEM_phi","h_recHits_EEM_phi",360,-3.1415927, 3.1415927);

    h_recHits_eta_MaxEt    = fs->make<TH1D>("h_recHits_eta_MaxEt","h_recHits_eta_MaxEt",150,-3.,3.);
    h_recHits_EB_phi_MaxEt = fs->make<TH1D>("h_recHits_EB_phi_MaxEt","h_recHits_EB_phi_MaxEt",360,-3.1415927, 3.1415927);
    h_recHits_EE_phi_MaxEt = fs->make<TH1D>("h_recHits_EE_phi_MaxEt","h_recHits_EE_phi_MaxEt",360,-3.1415927, 3.1415927);

    // ... all
    h_recHits_recoFlag = fs->make<TH1D>("h_recHits_recoFlag","h_recHits_recoFlag",16,-0.5,15.5);
/*//------------------------PFRechits Histogram defination------------------------------------
    h_PFRecHitEnergy          =fs->make<TH1D>("h_PFRecHitEnergy","h_PFrecHits_energy",10000,-1.,100.);
    h_PFRecHitTime              =fs->make<TH1D>("h_PFRecHitTime","h_PFrecHits_time",400,-100.,100.);
    h_PFRecHitDepth             =fs->make<TH1D>("h_PFRecHitDepth","h_PFrecHits_depth",10000, 0., 1.);

    h_PFRecHitEta             =fs->make<TH1D>("h_PFRecHitEta", "h_PFRecHitEta", 300, -4, 4);
    h_PFRecHitPhi             =fs->make<TH1D>("h_PFRecHitPhi", "h_PFRecHitPhi", 360, -3.14, +3.14);
    h_PFRecHitIEta            =fs->make<TH1D>("h_PFRecHitIEta", "h_PFRecHitIEta",300, -100, 100);

//------------------------....................................................................
*/
    if(usePFRecHitFlag_){
        cout<<"sjhfehihikkg"<<endl;
        PFrecHitCollection_EB_    = consumes<reco::PFRecHitCollection>(ps.getParameter<edm::InputTag>("PFrecHitCollection_EB"));
        PFrecHitCollection_EE_    = consumes<reco::PFRecHitCollection>(ps.getParameter<edm::InputTag>("PFrecHitCollection_EE"));

//------------------------PFRechits Histogram defination------------------------------------
        h_PFRecHitEnergy          =fs->make<TH1D>("h_PFRecHitEnergy","h_PFrecHits_energy",10000,-1.,100.);
        h_PFRecHitTime              =fs->make<TH1D>("h_PFRecHitTime","h_PFrecHits_time",400,-100.,100.);
        h_PFRecHitDepth             =fs->make<TH1D>("h_PFRecHitDepth","h_PFrecHits_depth",10000, 0., 1.);
        h_PFRecHitEta             =fs->make<TH1D>("h_PFRecHitEta", "h_PFRecHitEta", 300, -4, 4);
        h_PFRecHitPhi             =fs->make<TH1D>("h_PFRecHitPhi", "h_PFRecHitPhi", 360, -3.14, +3.14);
        h_PFRecHitIEta            =fs->make<TH1D>("h_PFRecHitIEta", "h_PFRecHitIEta",300, -100, 100);
    }

    // Basic Clusters ----------------------------------------------

    // ... barrel
    h_basicClusters_EB_size    = fs->make<TH1D>("h_basicClusters_EB_size","h_basicClusters_EB_size",200,0.,200.);
    h_basicClusters_EB_nXtals  = fs->make<TH1D>("h_basicClusters_EB_nXtals","h_basicClusters_EB_nXtals",400,0.,400.);
    h_basicClusters_EB_energy  = fs->make<TH1D>("h_basicClusters_EB_energy","h_basicClusters_EB_energy",2000,0.,400.);

    // ... barrel (with spike cleaning)
    h_basicClusters_EB_size_cleaned    = fs->make<TH1D>("h_basicClusters_EB_size_cleaned","h_basicClusters_EB_size_cleaned",200,0.,200.);
    h_basicClusters_EB_nXtals_cleaned  = fs->make<TH1D>("h_basicClusters_EB_nXtals_cleaned","h_basicClusters_EB_nXtals_cleaned",400,0.,400.);
    h_basicClusters_EB_energy_cleaned  = fs->make<TH1D>("h_basicClusters_EB_energy_cleaned","h_basicClusters_EB_energy_cleaned",2000,0.,400.);

    //barrel (with spike cleaning and track match)
    h_basicClusters_EB_dr_cleaned_tkmatched  = fs->make<TH1D>("h_basicClusters_EB_dr_cleaned_tkmatched","h_basicClusters_EB_dr_cleaned_tkmatched",100,0.,0.1);
    h_basicClusters_EB_size_cleaned_tkmatched    = fs->make<TH1D>("h_basicClusters_EB_size_cleaned_tkmatched","h_basicClusters_EB_size_cleaned_tkmatched",200,0.,200.);
    h_basicClusters_EB_nXtals_cleaned_tkmatched  = fs->make<TH1D>("h_basicClusters_EB_nXtals_cleaned_tkmatched","h_basicClusters_EB_nXtals_cleaned_tkmatched",400,0.,400.);
    h_basicClusters_EB_energy_cleaned_tkmatched  = fs->make<TH1D>("h_basicClusters_EB_energy_cleaned_tkmatched","h_basicClusters_EB_energy_cleaned_tkmatched",2000,0.,400.);

    // ... associated barrel rec hits
    h_basicClusters_recHits_EB_recoFlag = fs->make<TH1D>("h_basicClusters_recHits_EB_recoFlag","h_basicClusters_recHits_EB_recoFlag",16,-0.5,15.5);


    // ... endcap
    h_basicClusters_EEP_size   = fs->make<TH1D>("h_basicClusters_EEP_size","h_basicClusters_EEP_size",200,0.,200.);
    h_basicClusters_EEP_nXtals = fs->make<TH1D>("h_basicClusters_EEP_nXtals","h_basicClusters_EEP_nXtals",400,0.,400.);
    h_basicClusters_EEP_energy = fs->make<TH1D>("h_basicClusters_EEP_energy","h_basicClusters_EEP_energy",2000,0.,400.);

    h_basicClusters_EEM_size   = fs->make<TH1D>("h_basicClusters_EEM_size","h_basicClusters_EEM_size",200,0.,200.);
    h_basicClusters_EEM_nXtals = fs->make<TH1D>("h_basicClusters_EEM_nXtals","h_basicClusters_EEM_nXtals",400,0.,400.);
    h_basicClusters_EEM_energy = fs->make<TH1D>("h_basicClusters_EEM_energy","h_basicClusters_EEM_energy",2000,0.,400.);

    h_basicClusters_EEP_dr_tkmatched  = fs->make<TH1D>("h_basicClusters_EEP_dr_tkmatched","h_basicClusters_EEP_dr_tkmatched",100,0.,0.1);
    h_basicClusters_EEP_size_tkmatched   = fs->make<TH1D>("h_basicClusters_EEP_size_tkmatched","h_basicClusters_EEP_size_tkmatched",200,0.,200.);
    h_basicClusters_EEP_nXtals_tkmatched = fs->make<TH1D>("h_basicClusters_EEP_nXtals_tkmatched","h_basicClusters_EEP_nXtals_tkmatched",400,0.,400.);
    h_basicClusters_EEP_energy_tkmatched = fs->make<TH1D>("h_basicClusters_EEP_energy_tkmatched","h_basicClusters_EEP_energy_tkmatched",2000,0.,400.);
    h_basicClusters_EEP_occupancy_esmatched = fs->make<TH2D>
            ("h_basicClusters_EEP_occupancy_esmatched","h_basicClusters_EEP_occupancy_esmatched",150,-3.,3.,360,-3.1415927,3.1415927);
    h_basicClusters_EEP_eta_esmatched = fs->make<TH1D>
                                        ("h_basicClusters_EEP_eta_esmatched","h_basicClusters_EEP_eta_esmatched",150,-3.,3.);
    h_basicClusters_EEP_phi_esmatched = fs->make<TH1D>
                                        ("h_basicClusters_EEP_phi_esmatched","h_basicClusters_EEP_phi_esmatched",360,-3.1415927,3.1415927);

    h_basicClusters_EEM_dr_tkmatched  = fs->make<TH1D>("h_basicClusters_EEM_dr_tkmatched","h_basicClusters_EEM_dr_tkmatched",100,0.,0.1);
    h_basicClusters_EEM_size_tkmatched   = fs->make<TH1D>("h_basicClusters_EEM_size_tkmatched","h_basicClusters_EEM_size_tkmatched",200,0.,200.);
    h_basicClusters_EEM_nXtals_tkmatched = fs->make<TH1D>("h_basicClusters_EEM_nXtals_tkmatched","h_basicClusters_EEM_nXtals_tkmatched",400,0.,400.);
    h_basicClusters_EEM_energy_tkmatched = fs->make<TH1D>("h_basicClusters_EEM_energy_tkmatched","h_basicClusters_EEM_energy_tkmatched",2000,0.,400.);
    h_basicClusters_EEM_occupancy_esmatched = fs->make<TH2D>
            ("h_basicClusters_EEM_occupancy_esmatched","h_basicClusters_EEM_occupancy_esmatched",150,-3.,3.,360,-3.1415927,3.1415927);
    h_basicClusters_EEM_eta_esmatched = fs->make<TH1D>
                                        ("h_basicClusters_EEM_eta_esmatched","h_basicClusters_EEM_eta_esmatched",150,-3.,3.);
    h_basicClusters_EEM_phi_esmatched = fs->make<TH1D>
                                        ("h_basicClusters_EEM_phi_esmatched","h_basicClusters_EEM_phi_esmatched",360,-3.1415927,3.1415927);

    h_basicClusters_eta        = fs->make<TH1D>("h_basicClusters_eta","h_basicClusters_eta",150,-3.,3.);
    h_basicClusters_EB_eta     = fs->make<TH1D>("h_basicClusters_EB_eta","h_basicClusters_EB_eta",150,-3.,3.);
    h_basicClusters_EE_eta     = fs->make<TH1D>("h_basicClusters_EE_eta","h_basicClusters_EE_eta",150,-3.,3.);
    h_basicClusters_EB_phi     = fs->make<TH1D>("h_basicClusters_EB_phi","h_basicClusters_EB_phi",360,-3.1415927,3.1415927);
    h_basicClusters_EE_phi     = fs->make<TH1D>("h_basicClusters_EE_phi","h_basicClusters_EE_phi",360,-3.1415927,3.1415927);

    // ... associated endcap rec hits
    h_basicClusters_recHits_EE_recoFlag = fs->make<TH1D>("h_basicClusters_recHits_EE_recoFlag","h_basicClusters_recHits_EE_recoFlag",16,-0.5,15.5);

    //cleaned+tkmatched
    h_basicClusters_eta_tkmatched        = fs->make<TH1D>("h_basicClusters_eta_tkmatched","h_basicClusters_eta_tkmatched",150,-3.,3.);
    h_basicClusters_EB_eta_tkmatched     = fs->make<TH1D>("h_basicClusters_EB_eta_tkmatched","h_basicClusters_EB_eta_tkmatched",150,-3.,3.);
    h_basicClusters_EE_eta_tkmatched     = fs->make<TH1D>("h_basicClusters_EE_eta_tkmatched","h_basicClusters_EE_eta_tkmatched",150,-3.,3.);
    h_basicClusters_EB_phi_tkmatched     = fs->make<TH1D>("h_basicClusters_EB_phi_tkmatched","h_basicClusters_EB_phi_tkmatched",360,-3.1415927,3.1415927);
    h_basicClusters_EE_phi_tkmatched     = fs->make<TH1D>("h_basicClusters_EE_phi_tkmatched","h_basicClusters_EE_phi_tkmatched",360,-3.1415927,3.1415927);

    // ... associated all rec hits
    h_basicClusters_recHits_recoFlag = fs->make<TH1D>("h_basicClusters_recHits_recoFlag","h_basicClusters_recHits_recoFlag",16,-0.5,15.5);


    // Super Clusters ----------------------------------------------
    // ... barrel
    h_superClusters_EB_size    = fs->make<TH1D>("h_superClusters_EB_size","h_superClusters_EB_size",200,0.,200.);
    h_superClusters_EB_nXtals  = fs->make<TH1D>("h_superClusters_EB_nXtals","h_superClusters_EB_nXtals",400,0.,400.);
    h_superClusters_EB_nBC     = fs->make<TH1D>("h_superClusters_EB_nBC","h_superClusters_EB_nBC",100,0.,100.);
    h_superClusters_EB_energy  = fs->make<TH1D>("h_superClusters_EB_energy","h_superClusters_EB_energy",2000,0.,400.);
    h_superClusters_EB_E1oE4   = fs->make<TH1D>("h_superClusters_EB_E1oE4","h_superClusters_EB_E1oE4",150,0,1.5);
    h_superClusters_EB_occupancy     = fs->make<TH2D>("h_superClusters_EB_occupancy","h_superClusters_EB_occupancy",360,1.,361.,172,-86.,86. );

    // ... barrel (with spike cleaning)
    h_superClusters_EB_size_cleaned    = fs->make<TH1D>("h_superClusters_EB_size_cleaned","h_superClusters_EB_size_cleaned",200,0.,200.);
    h_superClusters_EB_nXtals_cleaned  = fs->make<TH1D>("h_superClusters_EB_nXtals_cleaned","h_superClusters_EB_nXtals_cleaned",400,0.,400.);
    h_superClusters_EB_nBC_cleaned     = fs->make<TH1D>("h_superClusters_EB_nBC_cleaned","h_superClusters_EB_nBC_cleaned",100,0.,100.);
    h_superClusters_EB_energy_cleaned  = fs->make<TH1D>("h_superClusters_EB_energy_cleaned","h_superClusters_EB_energy_cleaned",2000,0.,400.);

    h_superClusters_EB_rawEnergy_cleaned  = fs->make<TH1D>("h_superClusters_EB_rawEnergy_cleaned","h_superClusters_EB_rawEnergy_cleaned",2000,0.,400.);
    h_superClusters_EB_rawEt_cleaned      = fs->make<TH1D>("h_superClusters_EB_rawEt_cleaned","h_superClusters_EB_rawEt_cleaned",2000,0.,400.);

    // ... endcap
    h_superClusters_EEP_size   = fs->make<TH1D>("h_superClusters_EEP_size","h_superClusters_EEP_size",200,0.,200.);
    h_superClusters_EEP_nXtals = fs->make<TH1D>("h_superClusters_EEP_nXtals","h_superClusters_EEP_nXtals",400,0.,400.);
    h_superClusters_EEP_nBC    = fs->make<TH1D>("h_superClusters_EEP_nBC","h_superClusters_EEP_nBC",100,0.,100.);
    h_superClusters_EEP_energy = fs->make<TH1D>("h_superClusters_EEP_energy","h_superClusters_EEP_energy",2000,0.,400.);
    h_superClusters_EEP_rawEnergy = fs->make<TH1D>("h_superClusters_EEP_rawEnergy","h_superClusters_EEP_rawEnergy",2000,0.,400.);
    h_superClusters_EEP_rawEt = fs->make<TH1D>("h_superClusters_EEP_rawEt","h_superClusters_EEP_rawEt",2000,0.,400.);
    h_superClusters_EEP_E1oE4  = fs->make<TH1D>("h_superClusters_EEP_E1oE4","h_superClusters_EEP_E1oE4",150,0,1.5);
    h_superClusters_EEP_occupancy     = fs->make<TH2D>("h_superClusters_EEP_occupancy","h_superClusters_EEP_occupancy",100,0.,100.,100,0.,100.);

    h_superClusters_EEM_size   = fs->make<TH1D>("h_superClusters_EEM_size","h_superClusters_EEM_size",200,0.,200.);
    h_superClusters_EEM_nXtals = fs->make<TH1D>("h_superClusters_EEM_nXtals","h_superClusters_EEM_nXtals",400,0.,400.);
    h_superClusters_EEM_nBC    = fs->make<TH1D>("h_superClusters_EEM_nBC","h_superClusters_EEM_nBC",100,0.,100.);
    h_superClusters_EEM_energy = fs->make<TH1D>("h_superClusters_EEM_energy","h_superClusters_EEM_energy",2000,0.,400.);
    h_superClusters_EEM_rawEnergy = fs->make<TH1D>("h_superClusters_EEM_rawEnergy","h_superClusters_EEM_rawEnergy",2000,0.,400.);
    h_superClusters_EEM_rawEt = fs->make<TH1D>("h_superClusters_EEM_rawEt","h_superClusters_EEM_rawEt",2000,0.,400.);
    h_superClusters_EEM_E1oE4  = fs->make<TH1D>("h_superClusters_EEM_E1oE4","h_superClusters_EEM_E1oE4",150,0,1.5);
    h_superClusters_EEM_occupancy     = fs->make<TH2D>("h_superClusters_EEM_occupancy","h_superClusters_EEM_occupancy",100,0.,100.,100,0.,100.);

    h_superClusters_occupancyPhiEta = fs->make<TH2D>("h_superClusters_occupancyPhiEta","h_superClusters_occupancyPhiEta",360,-3.1415927,3.1415927,150,-3.,3.);

    h_superClusters_eta        = fs->make<TH1D>("h_superClusters_eta","h_superClusters_eta",150,-3.,3.);
    h_superClusters_eta_scCut  = fs->make<TH1D>("h_superClusters_eta_scCut","h_superClusters_eta_scCut",150,-3.,3.);
    h_superClusters_EB_eta     = fs->make<TH1D>("h_superClusters_EB_eta","h_superClusters_EB_eta",150,-3.,3.);
    h_superClusters_EE_eta     = fs->make<TH1D>("h_superClusters_EE_eta","h_superClusters_EE_eta",150,-3.,3.);
    h_superClusters_EB_phi     = fs->make<TH1D>("h_superClusters_EB_phi","h_superClusters_EB_phi",360,-3.1415927,3.1415927);
    h_superClusters_EE_phi     = fs->make<TH1D>("h_superClusters_EE_phi","h_superClusters_EE_phi",360,-3.1415927,3.1415927);

    h2_superClusters_EB_seedTimeVsEnergy = fs->make<TH2D>("h2_superClusters_EB_seedTimeVsEnergy","h2_superClusters_EB_seedTimeVsEnergy",2000,0.,400.,400,-100.,100.);
    h2_superClusters_EE_seedTimeVsEnergy = fs->make<TH2D>("h2_superClusters_EE_seedTimeVsEnergy","h2_superClusters_EE_seedTimeVsEnergy",2000,0.,400.,400,-100.,100.);

    // preshower
    h_esClusters_energy_plane1 = fs->make<TH1D>("h_esClusters_energy_plane1","h_esClusters_energy_plane1",1000,0.,0.01);
    h_esClusters_energy_plane2 = fs->make<TH1D>("h_esClusters_energy_plane2","h_esClusters_energy_plane2",1000,0.,0.01);
    h_esClusters_energy_ratio  = fs->make<TH1D>("h_esClusters_energy_ratio","h_esClusters_energy_ratio",100,0.,10.);

    // Pi0 ----------------------------------------------
    // ... barrel
    h_Pi0_EB_mass   = dirPi0.make<TH1D>("h_Pi0_EB_mass","h_Pi0_EB_mass",100,0.,0.5);
    h_Pi0_EB_pt1    = dirPi0.make<TH1D>("h_Pi0_EB_pt1","h_Pi0_EB_pt1",100,0.,20.);
    h_Pi0_EB_pt2    = dirPi0.make<TH1D>("h_Pi0_EB_pt2","h_Pi0_EB_pt2",100,0.,20.);
    h_Pi0_EB_pt     = dirPi0.make<TH1D>("h_Pi0_EB_pt","h_Pi0_EB_pt",100,0.,20.);
    h_Pi0_EB_eta    = dirPi0.make<TH1D>("h_Pi0_EB_eta","h_Pi0_EB_eta",100,-1.5,1.5);
    h_Pi0_EB_phi    = dirPi0.make<TH1D>("h_Pi0_EB_phi","h_Pi0_EB_phi",360,-3.,3.);

    // ... endcap
    h_Pi0_EE_mass   = dirPi0.make<TH1D>("h_Pi0_EE_mass","h_Pi0_EE_mass",100,0.,0.5);
    h_Pi0_EE_pt1    = dirPi0.make<TH1D>("h_Pi0_EE_pt1","h_Pi0_EE_pt1",100,0.,20.);
    h_Pi0_EE_pt2    = dirPi0.make<TH1D>("h_Pi0_EE_pt2","h_Pi0_EE_pt2",100,0.,20.);
    h_Pi0_EE_pt     = dirPi0.make<TH1D>("h_Pi0_EE_pt","h_Pi0_EE_pt",100,0.,20.);
    h_Pi0_EE_eta    = dirPi0.make<TH1D>("h_Pi0_EE_eta","h_Pi0_EE_eta",100,-3.1415927,3.1415927);
    h_Pi0_EE_phi    = dirPi0.make<TH1D>("h_Pi0_EE_phi","h_Pi0_EE_phi",360,-3.1415927,3.1415927);

    // Jets ----------------------------------------------
    // ... barrel
    h_Jets_EB_emf        =  fs->make<TProfile2D>("h_Jets_EB_emf","h_Jets_EB_emf",150,-3.,3.,360,-3.1415927,3.1415927);
    h_Jets_EB_emf_eta    =  fs->make<TProfile>("h_Jets_EB_emf_eta","h_Jets_EB_emf_eta",150,-3.,3.);
    h_Jets_EB_emf_phi    =  fs->make<TProfile>("h_Jets_EB_emf_phi","h_Jets_EB_emf_phi",360,-3.1415927,3.1415927);

    // ... Endcaps
    h_Jets_EEP_emf        =  fs->make<TProfile2D>("h_Jets_EEP_emf","h_Jets_EEP_emf",150,-3.,3.,360,-3.1415927,3.1415927);
    h_Jets_EEP_emf_eta    =  fs->make<TProfile>("h_Jets_EEP_emf_eta","h_Jets_EEP_emf_eta",150,-3.,3.);
    h_Jets_EEP_emf_phi    =  fs->make<TProfile>("h_Jets_EEP_emf_phi","h_Jets_EEP_emf_phi",360,-3.1415927,3.1415927);

    h_Jets_EEM_emf        =  fs->make<TProfile2D>("h_Jets_EEM_emf","h_Jets_EEM_emf",150,-3.,3.,360,-3.1415927,3.1415927);
    h_Jets_EEM_emf_eta    =  fs->make<TProfile>("h_Jets_EEM_emf_eta","h_Jets_EEM_emf_eta",150,-3.,3.);
    h_Jets_EEM_emf_phi    =  fs->make<TProfile>("h_Jets_EEM_emf_phi","h_Jets_EEM_emf_phi",360,-3.1415927,3.1415927);


    cout << "done constructor" << endl;
}



EcalValidationAOD::~EcalValidationAOD()
{
    // do anything here that needs to be done at desctruction time
    // (e.g. close files, deallocate resources etc.)
}


//
// member functions
//

// ------------ method called to for each event  ------------
void EcalValidationAOD::analyze(const edm::Event& ev, const edm::EventSetup& iSetup)
{
    if(debug) cout << "begin analyze function" << endl;
    nPU_   .clear();
    puTrue_.clear();
    puBX_  .clear();
//    float runNumber = ev.run();   //ps
//    std::cout<<"RUN NUMBER"<<runNumber<<std::endl; //ps
//   float bx      = ev.bunchCrossing();
//   float ls      = ev.luminosityBlock();
//   float orbitNb = ev.orbitNumber();

    //Get the BS position
    edm::Handle<reco::BeamSpot> recoBeamSpotHandle;
    ev.getByToken(beamSpot_,recoBeamSpotHandle);
    const reco::BeamSpot::Point& BSPosition = recoBeamSpotHandle->position();
    // Vertex Collection
    edm::Handle<reco::VertexCollection> vertexes;
    ev.getByToken(PV_,vertexes);
    int nvertices = 0;
    for(reco::VertexCollection::const_iterator it = vertexes->begin(); it !=vertexes->end() ; ++it) {
        if(it->isValid() && it->ndof() > 4. && it->position().Rho() < 2. && fabs(it->position().Z() - BSPosition.Z()) < 24) {
            nvertices++;
        }
    }/*vertex loop*/
//PileUp interactions
//cout<<"nvertices"<<nvertices<<endl;
    if(vertexes->size() ==1 && debug) cout << "cleaned verticies " << nvertices << endl;
//    if(!nvertices) return;
    double MyWeight(0);
    if(isMC) {
        edm::Handle<vector<PileupSummaryInfo> > genPileupHandle;
        ev.getByToken(puCollection_, genPileupHandle);
        if (genPileupHandle.isValid()) {
            for (vector<PileupSummaryInfo>::const_iterator pu = genPileupHandle->begin(); pu != genPileupHandle->end(); ++pu) {
                if (pu->getBunchCrossing() == 0) {
                    hPU_->Fill(pu->getPU_NumInteractions());
                    hPUTrue_->Fill(pu->getTrueNumInteractions());
//cout << "PU " << pu->getTrueNumInteractions() << endl;
//double MyWeight = LumiWeights_.weight( pu->getPU_NumInteractions() );
//cout << MyWeight << "\t" <<  LumiWeights_.weight( pu->getPU_NumInteractions()+20 ) << "\t" <<  LumiWeights_.weight( pu->getPU_NumInteractions()+40 ) << endl;
//cout <<MyWeight << "\t" <<  LumiWeights[pu->getPU_NumInteractions()+20] << "\t" <<  LumiWeights[pu->getPU_NumInteractions()+40] << endl;
                }

                nPU_   .push_back(pu->getPU_NumInteractions());
                puTrue_.push_back(pu->getTrueNumInteractions());
                puBX_  .push_back(pu->getBunchCrossing());

                nPUInfo_++;
            }
        }
        else
            edm::LogWarning("ggNtuplizer") << "no PileupSummaryInfo in event";
        if(nvertices < 200) MyWeight = LumiWeights[nvertices];
        if(nvertices) h_PV_n->Fill(nvertices, MyWeight);
        if(nvertices) h_PV_nw->Fill(nvertices);
//cout << "number of vertices " << nvertices << " " << h_PV_n->Integral() << endl;
    }
    else {
        MyWeight = 1;
        if(nvertices) h_PV_n->Fill(nvertices, MyWeight);
        if(nvertices) h_PV_nw->Fill(nvertices);
//cout << "number of vertices " << nvertices << " " << h_PV_n->Integral() << endl;
    }
    //Get tracks
    edm::Handle<edm::View<reco::Track> > TracksHandle ;
    ev.getByToken (tracks_, TracksHandle) ;
    //Get the magnetic field
    edm::ESHandle<MagneticField> theMagField;
    iSetup.get<IdealMagneticFieldRecord>().get(theMagField);
    //DataBase noise
    edm::ESHandle<EcalPedestals> dbPed   ;
    iSetup.get<EcalPedestalsRcd>().get( dbPed ) ;
    const EcalPedestals* pedestals ( dbPed.product() ) ;
    //DataBase LaserCorrection
    edm::ESHandle<EcalLaserDbService> theLaser;
    iSetup.get<EcalLaserDbRecord>().get(theLaser);

    naiveId_++;
//cout<<"Event no"<<naiveId_<<endl;

    if(ttt%1000 ==0) {
        cout << "MyWeight" << MyWeight <<endl;
        cout<< "clean/all " << nvertices << "\t" << vertexes->size() << endl;
    }
    ++ttt;
    //   ev.getByToken(ebEcalDigiCollection_, SRFlagsEB );
    //
    //   for(EBSrFlagCollection::const_iterator it = SRFlagsEB->begin(); it != SRFlagsEB->end(); ++it)
    //     {
    //       const int flag = it->value();
    //       if( flag != EcalSrFlag::SRF_FULL ) continue;
    //       h_recHits_EB_SRP->Fill(1);
    //     }
    // calo geometry
    edm::ESHandle<CaloGeometry> pGeometry;
    iSetup.get<CaloGeometryRecord>().get(pGeometry);
    const CaloGeometry *geometry = pGeometry.product();

//   const CaloSubdetectorGeometry *geometry_EB = geometry->getSubdetectorGeometry(DetId::Ecal,EcalBarrel);
//   const CaloSubdetectorGeometry *geometry_EE = geometry->getSubdetectorGeometry(DetId::Ecal,EcalEndcap);

    // calo topology
    edm::ESHandle<CaloTopology> pTopology;
    iSetup.get<CaloTopologyRecord>().get(pTopology);
    const CaloTopology *topology = pTopology.product();

    /*// --- REDUCED REC HITS -------------------------------------------------------------------------------------
    edm::Handle<EcalRecHitCollection> redRecHitsEB;
    ev.getByToken( redRecHitCollection_EB_, redRecHitsEB );
    const EcalRecHitCollection* theBarrelEcalredRecHits = redRecHitsEB.product () ;
    if ( ! redRecHitsEB.isValid() ) {
      std::cerr << "EcalValidation::analyze --> redRecHitsEB not found" << std::endl;
    }

    for ( EcalRecHitCollection::const_iterator itr = theBarrelEcalredRecHits->begin () ;
          itr != theBarrelEcalredRecHits->end () ;++itr)
    {

      h_redRecHits_EB_recoFlag      -> Fill( itr -> recoFlag() );
      h_redRecHits_recoFlag         -> Fill( itr -> recoFlag() );

    }

    // ... endcap
    edm::Handle<EcalRecHitCollection> redRecHitsEE;
    ev.getByToken( redRecHitCollection_EE_, redRecHitsEE );
    const EcalRecHitCollection* theEndcapEcalredRecHits = redRecHitsEE.product () ;
    if ( ! redRecHitsEE.isValid() ) {
      std::cerr << "EcalValidation::analyze --> redRecHitsEE not found" << std::endl;
    }

    for ( EcalRecHitCollection::const_iterator itr = theEndcapEcalredRecHits->begin () ;
          itr != theEndcapEcalredRecHits->end () ; ++itr)
    {

        EEDetId eeid( itr -> id() );
        h_redRecHits_EE_recoFlag       -> Fill( itr -> recoFlag() );
        h_redRecHits_recoFlag          -> Fill( itr -> recoFlag() );

    }*/

    // --- REC HITS -------------------------------------------------------------------------------------

    // ... barrel
    edm::Handle<EcalRecHitCollection> recHitsEB;
    ev.getByToken( recHitCollection_EB_, recHitsEB );
    const EcalRecHitCollection* theBarrelEcalRecHits = recHitsEB.product () ;
    if ( ! recHitsEB.isValid() ) {
        std::cerr << "EcalValidation::analyze --> recHitsEB not found" << std::endl;
    }

    float maxRecHitEnergyEB        = -999.;
    float maxRecHitEtaEB           = -999.;
    float maxRecHitPhiEB           = -999.;
    float maxRecHitIEtaEB          = -999.;


    float maxRecHitEnergyEBcleaned = -999.;

    float maxEtRecHitEnergyEB = -999.;
    float maxEtRecHitEtaEB    = -999.;
    float maxEtRecHitPhiEB    = -999.;

    EBDetId ebid_MrecHitEB;

    bool  hasSpike = false;

    int gain = int(gainId_);
    int nhits_passingEB(0);
    int nhits_passingEE(0);
//int nhits_passingEEP(0);
//int nhits_passingEEM(0);
    for ( EcalRecHitCollection::const_iterator itr = theBarrelEcalRecHits->begin () ;
            itr != theBarrelEcalRecHits->end () ; ++itr)
    {

        EBDetId ebid( itr -> id() );
        GlobalPoint mycell = geometry -> getPosition(DetId(itr->id()));
        double et = itr -> energy()*mycell.perp()/mycell.mag();

        if(itr->checkFlag(itr->kOutOfTime)) continue;
        if(itr->checkFlag(itr->kWeird)) continue;
//cut on the Et to  be greater than certain threshold.
//cout << "et  " << et << " enerny " << itr->energy() << endl;
        if (et < ethrEB_ ) continue;
	++nhits_passingEB;
	float R4 = EcalTools::swissCross( ebid, *theBarrelEcalRecHits, 0. );
//std::cout << "et_after cut" << et << std::endl;
//std::cout << "R4=    " << R4 << std::endl;
//std::cout << "mycell.eta()=    " << mycell.eta() << std::endl;
//std::cout << "gain=    " << gain << std::endl;
        if ( itr -> energy() > 3. && abs(ebid.ieta())!=85 )  h_recHits_EB_E1oE4-> Fill( R4 , MyWeight);


        if ( R4 > 0.95 ) {
            h_recHits_EB_energy_spike       -> Fill( itr -> energy(), MyWeight);
        }

        // spike cleaning
        if ( R4 > 0.95 && et > 3.) hasSpike= true;
        if ( R4 > 0.95 && et > 3.) continue;
        if (hasSpike) continue;
        sumRecHitEtEB[naiveId_-1] += et;

        EcalPedestalsMap::const_iterator itped = pedestals->getMap().find( ebid );
        float width = (*itped).rms(gain);
//std::cout << "width=    " << width << std::endl;
        //std::cout << "EB RMS-event: " << width << " - " << ebid.iphi() << " - "  << ebid.ieta() << " - " <<  naiveId_ << " - " << gain << std::endl;

        if(h_DB_noiseMap_EB->GetBinContent(ebid.iphi()-1,ebid.ieta()+86) == 0) h_DB_noiseMap_EB->SetBinContent(ebid.iphi()-1,ebid.ieta()+86,width);

        float seedLaserCorrection = theLaser->getLaserCorrection(ebid, ev.time());
//std::cout << "seedLaserCorrection=    " << seedLaserCorrection << std::endl;
//std::cout << "ev.time()"<< ev.time() <<endl;

        if(h_DB_LaserCorrMap_EB->GetBinContent(ebid.iphi()-1,ebid.ieta()+86) == 0) h_DB_LaserCorrMap_EB->SetBinContent(ebid.iphi()-1,ebid.ieta()+86,seedLaserCorrection);
        h_recHits_EB_energy        -> Fill( itr -> energy() , MyWeight);
        h_recHits_EB_recoFlag      -> Fill( itr -> recoFlag() , MyWeight);
        h_recHits_recoFlag         -> Fill( itr -> recoFlag() , MyWeight);
        // max E rec hit
        if (itr -> energy() > maxRecHitEnergyEB ) {
            maxRecHitEnergyEB = itr -> energy() ;
            maxRecHitEtaEB = mycell.eta();
            maxRecHitPhiEB = mycell.phi();
            maxRecHitIEtaEB = ebid.ieta();
        }
        if ( et > maxEtRecHitEnergyEB ) {
            maxEtRecHitEnergyEB = et ;
            maxEtRecHitEtaEB    = mycell.eta();
            maxEtRecHitPhiEB    = mycell.phi();
            ebid_MrecHitEB      = ebid;
        }
//     cout << "looking for out of time rechit " << itr->checkFlag(itr->kOutOfTime) << endl;
//if(itr->checkFlag(itr->kOutOfTime)) cout <<"energy & time " << itr->time() << "\t" << itr->energy() << "\t" << itr->chi2() << endl;
//cout << "  kGood                " << itr->checkFlag(itr->kGood)          			<< endl;
////cout << "  kGood                " << itr->kGood          			<< endl;
//cout << "  kPoorReco            " << itr->checkFlag(itr->kPoorReco           	)	 << endl;
//cout << "  kOutOfTime           " << itr->checkFlag(itr->kOutOfTime           	)	 << endl;
//cout << "  kFaultyHardware      " << itr->checkFlag(itr->kFaultyHardware       	)	 << endl;
//cout << "  kNoisy               " << itr->checkFlag(itr->kNoisy               	) 		 << endl;
//cout << "  kPoorCalib           " << itr->checkFlag(itr->kPoorCalib           	)	 		 << endl;
//cout << "  kSaturated           " << itr->checkFlag(itr->kSaturated           	)	 		 << endl;
//cout << "  kLeadingEdgeRecovered" << itr->checkFlag(itr->kLeadingEdgeRecovered	)	 		 << endl;
//cout << "  kNeighboursRecovered " << itr->checkFlag(itr->kNeighboursRecovered 	)	 		 << endl;
//cout << "  kTowerRecovered      " << itr->checkFlag(itr->kTowerRecovered      	)	 		 << endl;
//cout << "  kDead                " << itr->checkFlag(itr->kDead                	)	 		 << endl;
//cout << "  kKilled              " << itr->checkFlag(itr->kKilled              	)	 		 << endl;
//cout << "  kTPSaturated         " << itr->checkFlag(itr->kTPSaturated         	)	 		 << endl;
//cout << "  kL1SpikeFlag         " << itr->checkFlag(itr->kL1SpikeFlag         	)	 		 << endl;
//cout << "  kWeird               " << itr->checkFlag(itr->kWeird               	)	 		 << endl;
//cout << "  kDiWeird             " << itr->checkFlag(itr->kDiWeird             	)	 		 << endl;
//cout << "  kHasSwitchToGain6    " << itr->checkFlag(itr->kHasSwitchToGain6    	)	 		 << endl;
//cout << "  kHasSwitchToGain1    " << itr->checkFlag(itr->kHasSwitchToGain1    	)	 		 << endl;
//cout << "  kUnknown             "  <<itr->checkFlag(itr-> kUnknown              )		 		 << endl;

        if ( et > ethrEB_ ) {
            h_recHits_EB_time          -> Fill( itr -> time(), MyWeight);
            h_recHits_EB_Eta->Fill(mycell.eta(), MyWeight);
            h_recHits_EB_Phi->Fill(mycell.phi(), MyWeight);
            h_recHits_EB_Chi2          -> Fill( itr -> chi2(), MyWeight);
            h_recHits_EB_kOutOfTime->Fill(itr->checkFlag(itr->kOutOfTime));
            h_recHits_EB_kWeird->Fill(itr->checkFlag(itr->kWeird));
            h_recHits_EB_energy_thr    -> Fill( itr -> energy(), MyWeight);
            if( !(itr->checkFlag(itr->kOutOfTime))   &&  !(itr->checkFlag(itr->kWeird)))
                h_recHits_EB_time_cut          -> Fill( itr -> time(), MyWeight);
            if( !(itr->checkFlag(itr->kOutOfTime))   &&  !(itr->checkFlag(itr->kWeird)))
                h_recHits_EB_Chi2_cut          -> Fill( itr -> chi2(), MyWeight);
            //h_recHits_EB_OutOfTimeChi2 -> Fill( itr -> outOfTimeChi2(), MyWeight);
            h_recHits_EB_occupancy     -> Fill( ebid.iphi() , ebid.ieta(), MyWeight);
            h_recHits_EB_deviation     -> Fill( ebid.iphi() , ebid.ieta(), MyWeight);
            h_recHits_EB_iPhiOccupancy -> Fill( ebid.iphi(), MyWeight);
            h_recHits_EB_iEtaOccupancy -> Fill( ebid.ieta(), MyWeight);
        }


        h_recHits_EB_energy_cleaned -> Fill( itr -> energy(), MyWeight);
        if ( (et > ethrEB_)  && !hasSpike) {
            h_recHits_EB_time_cleaned          -> Fill( itr -> time(), MyWeight);
            h_recHits_EB_Chi2_cleaned          -> Fill( itr -> chi2(), MyWeight);
            //h_recHits_EB_OutOfTimeChi2_cleaned -> Fill( itr -> outOfTimeChi2(), MyWeight);
        }

        // max E rec hit - cleaned
        if (itr -> energy() > maxRecHitEnergyEBcleaned ) {
            maxRecHitEnergyEBcleaned = itr -> energy() ;
        }
    }
//   cout << "total hits passing EB thresholds " << nhits_passingEB << endl;
    if(maxRecHitEnergyEB >=8.0) {
        h_recHits_EB_size           -> Fill( nhits_passingEB, MyWeight);
        if (!hasSpike) h_recHits_EB_size_cleaned -> Fill( nhits_passingEB, MyWeight);

        h_recHits_EB_energyMax         -> Fill( maxRecHitEnergyEB , MyWeight);
        h_recHits_EB_energyMax_cleaned -> Fill( maxRecHitEnergyEBcleaned , MyWeight);
    }
    if(debug)  cout << "done with EB RecHits " << endl;
    // ... endcap
    edm::Handle<EcalRecHitCollection> recHitsEE;
    ev.getByToken( recHitCollection_EE_, recHitsEE );
    const EcalRecHitCollection* theEndcapEcalRecHits = recHitsEE.product () ;
    if ( ! recHitsEE.isValid() ) {
        std::cerr << "EcalValidation::analyze --> recHitsEE not found" << std::endl;
    }

    int nHitsEE = 0;
    int nHitsEEP = 0;
    int nHitsEEM = 0;

    float maxRecHitEnergyEEP        = -999.;
    float maxRecHitEnergyEEM        = -999.;
    float maxRecHitEtaEEP           = -999.;
    float maxRecHitEtaEEM           = -999.;
    float maxRecHitPhiEEP           = -999.;
    float maxRecHitPhiEEM           = -999.;

    float maxEtRecHitEnergyEEP = -999.;
    float maxEtRecHitEnergyEEM = -999.;
    float maxEtRecHitEtaEEP    = -999.;
    float maxEtRecHitEtaEEM    = -999.;
    float maxEtRecHitPhiEEP    = -999.;
    float maxEtRecHitPhiEEM    = -999.;

    EEDetId eeid_MrecHitEEM;
    EEDetId eeid_MrecHitEEP;

    for ( EcalRecHitCollection::const_iterator itr = theEndcapEcalRecHits->begin () ;
            itr != theEndcapEcalRecHits->end () ; ++itr)
    {

        EEDetId eeid( itr -> id() );
        GlobalPoint mycell = geometry->getPosition(itr->detid());

        EcalPedestalsMap::const_iterator itped = pedestals->getMap().find( eeid );
        float width = (*itped).rms(gain);

        //if(width > 3) std::cout << "EE RMS-event: " << width << " - " << eeid.ix()-0.5 << " - " << eeid.iy()-0.5 << " - " <<  naiveId_ << " - " << gain << std::endl;


        if ( eeid.zside() > 0 ) {
            if(h_DB_noiseMap_EEP->GetBinContent(eeid.ix()-0.5,eeid.iy()-0.5) == 0) h_DB_noiseMap_EEP->SetBinContent(eeid.ix()-0.5,eeid.iy()-0.5,width);
        }

        if ( eeid.zside() < 0 ) {
            if(h_DB_noiseMap_EEM->GetBinContent(eeid.ix()-0.5,eeid.iy()-0.5) == 0) h_DB_noiseMap_EEM->SetBinContent(eeid.ix()-0.5,eeid.iy()-0.5,width);
        }

        float seedLaserCorrection = theLaser->getLaserCorrection(eeid, ev.time());
        if(h_DB_LaserCorrMap_EEP->GetBinContent(eeid.ix()-0.5,eeid.iy()-0.5) == 0 && eeid.zside() > 0) h_DB_LaserCorrMap_EEP->SetBinContent(eeid.ix()-0.5,eeid.iy()-0.5,seedLaserCorrection);
        if(h_DB_LaserCorrMap_EEM->GetBinContent(eeid.ix()-0.5,eeid.iy()-0.5) == 0 && eeid.zside() < 0) h_DB_LaserCorrMap_EEM->SetBinContent(eeid.ix()-0.5,eeid.iy()-0.5,seedLaserCorrection);

        double et = itr -> energy()*mycell.perp()/mycell.mag();
        if(et < 8 ) continue;
        if(debug) cout << "2et of EE hits " << et << endl;
        if(itr->checkFlag(itr->kOutOfTime)) continue;
        if(itr->checkFlag(itr->kWeird)) continue;
        if(debug) cout << "ee hits " << itr->checkFlag(itr->kOutOfTime) << "\t" << itr->checkFlag(itr->kWeird) << endl;
        ++nHitsEE;
        ++nhits_passingEE;
        if ( eeid.zside() > 0 ) {
            sumRecHitEtEEP[naiveId_-1] += et;
            if(fabs(mycell.eta()) < 2.5) sumRecHitEtCutEEP[naiveId_-1] += et;
        }
        if ( eeid.zside() < 0 ) {
            sumRecHitEtEEM[naiveId_-1] += et;
            if(fabs(mycell.eta()) < 2.5) sumRecHitEtCutEEM[naiveId_-1] += et;
        }

        // EE+
        if ( eeid.zside() > 0 ) {

            nHitsEEP++;
            h_recHits_EEP_energy        -> Fill( itr -> energy(), MyWeight);
            h_recHits_EE_recoFlag       -> Fill( itr -> recoFlag(), MyWeight);
            h_recHits_recoFlag          -> Fill( itr -> recoFlag(), MyWeight);

            // max E rec hit

            // max Et rec hit
            if (et > maxEtRecHitEnergyEEP &&
                    !(eeid.ix()>=41 && eeid.ix()<=60 && eeid.iy()>=41 && eeid.iy()<=60) ) {
                maxEtRecHitEnergyEEP = et ;
                maxEtRecHitEtaEEP    = mycell.eta() ;
                maxEtRecHitPhiEEP    = mycell.phi() ;
                eeid_MrecHitEEP      = eeid ;
            }

            // only channels above noise
            //if (  itr -> energy() > ethrEE_ ){
            h_recHits_EEP_time          -> Fill( itr -> time(), MyWeight);
            h_recHits_EEP_Chi2          -> Fill( itr -> chi2(), MyWeight);
            //h_recHits_EEP_OutOfTimeChi2 -> Fill( itr -> outOfTimeChi2(), MyWeight);
            h_recHits_EEP_occupancy     -> Fill( eeid.ix() - 0.5, eeid.iy() - 0.5, MyWeight);
            h_recHits_EEP_deviation     -> Fill( eeid.ix() - 0.5, eeid.iy() - 0.5, MyWeight);
            h_recHits_EEP_iXoccupancy   -> Fill( eeid.ix() - 0.5, MyWeight);
            h_recHits_EEP_iYoccupancy   -> Fill( eeid.iy() - 0.5, MyWeight);
            //}
        }

        // EE-
        if ( eeid.zside() < 0 ) {

            nHitsEEM++;
            h_recHits_EEM_energy        -> Fill( itr -> energy(), MyWeight);
            h_recHits_EE_recoFlag       -> Fill( itr -> recoFlag(), MyWeight);
            h_recHits_recoFlag          -> Fill( itr -> recoFlag(), MyWeight);

            // max E rec hit
            if (itr -> energy() > maxRecHitEnergyEEM &&
                    !(eeid.ix()>=41 && eeid.ix()<=60 && eeid.iy()>=41 && eeid.iy()<=60) ) {
                maxRecHitEnergyEEM = itr -> energy() ;
                maxRecHitEtaEEM    = mycell.eta() ;
                maxRecHitPhiEEM    = mycell.phi() ;
            }
            if (et > maxEtRecHitEnergyEEM &&
                    !(eeid.ix()>=41 && eeid.ix()<=60 && eeid.iy()>=41 && eeid.iy()<=60)) {
                maxEtRecHitEnergyEEM = et ;
                maxEtRecHitEtaEEM    = mycell.eta() ;
                maxEtRecHitPhiEEM    = mycell.phi() ;
                eeid_MrecHitEEM      = eeid ;
            }

            // only channels above noise
            if (  itr -> energy() > ethrEE_ ) {
                h_recHits_EEM_time          -> Fill( itr -> time(), MyWeight);
                h_recHits_EEM_Chi2          -> Fill( itr -> chi2(), MyWeight);
                //h_recHits_EEM_OutOfTimeChi2 -> Fill( itr -> outOfTimeChi2(), MyWeight);
                h_recHits_EEM_occupancy     -> Fill( eeid.ix()- 0.5, eeid.iy() - 0.5, MyWeight);
                h_recHits_EEM_deviation     -> Fill( eeid.ix()- 0.5, eeid.iy() - 0.5, MyWeight);
                h_recHits_EEM_iXoccupancy   -> Fill( eeid.ix() - 0.5, MyWeight);
                h_recHits_EEM_iYoccupancy   -> Fill( eeid.iy() - 0.5, MyWeight);

            }

       }

    } // end loop over EE rec hits

    if(debug)  cout << "done with EE RecHits " << endl;

    // size
    if(debug)   cout << "total hits passing EE thresholds " << nhits_passingEE << endl;
    h_recHits_EE_size    -> Fill( nhits_passingEE, MyWeight);
    if(recHitsEE->size() ==0)
//cout << "event number "  << ev.id().event()<<"\t" <<ev.id().luminosityBlock() << "\t";
        h_recHits_EEP_size   -> Fill( nHitsEEP , MyWeight);
    h_recHits_EEM_size   -> Fill( nHitsEEM , MyWeight);

    // energy
    h_recHits_EEP_energyMax -> Fill( maxRecHitEnergyEEP , MyWeight);
    h_recHits_EEM_energyMax -> Fill( maxRecHitEnergyEEM , MyWeight);


    // eta/phi for max E rec hit

    h_recHits_eta     -> Fill( maxRecHitEtaEB  , MyWeight);
    h_recHits_eta     -> Fill( maxRecHitEtaEEP , MyWeight);
    h_recHits_eta     -> Fill( maxRecHitEtaEEM , MyWeight);

    if (maxRecHitIEtaEB>0)   h_recHits_EB_eta -> Fill((maxRecHitIEtaEB -0.5)*0.0175, MyWeight);
    if (maxRecHitIEtaEB<0)   h_recHits_EB_eta -> Fill((maxRecHitIEtaEB +0.5)*0.0175, MyWeight);

    h_recHits_EEP_eta -> Fill( maxRecHitEtaEEP , MyWeight);
    h_recHits_EEM_eta -> Fill( maxRecHitEtaEEM , MyWeight);

    h_recHits_EB_phi  -> Fill( maxRecHitPhiEB  , MyWeight);
    h_recHits_EE_phi  -> Fill( maxRecHitPhiEEP , MyWeight);
    h_recHits_EE_phi  -> Fill( maxRecHitPhiEEM , MyWeight);

    h_recHits_EEP_phi -> Fill( maxRecHitPhiEEP , MyWeight);
    h_recHits_EEM_phi -> Fill( maxRecHitPhiEEM , MyWeight);



    // histos for max Et rec hit
    float noiseThr = 1.;

    if ( maxEtRecHitEnergyEB > noiseThr )  {
        h_recHits_eta_MaxEt    -> Fill( maxEtRecHitEtaEB, MyWeight) ;
        h_recHits_EB_phi_MaxEt -> Fill( maxEtRecHitPhiEB, MyWeight);
    }


    if ( maxEtRecHitEnergyEEP > noiseThr )  {
        h_recHits_eta_MaxEt    -> Fill( maxEtRecHitEtaEEP , MyWeight) ;
        h_recHits_EE_phi_MaxEt -> Fill( maxEtRecHitPhiEEP , MyWeight);
    }

    if ( maxEtRecHitEnergyEEM > noiseThr )  {
        h_recHits_eta_MaxEt     -> Fill( maxEtRecHitEtaEEM , MyWeight) ;
        h_recHits_EE_phi_MaxEt  -> Fill( maxEtRecHitPhiEEM , MyWeight);
    }


//----------------------------------PFrechits-------------------------------
    if(usePFRecHitFlag_){
    edm::Handle<reco::PFRecHitCollection> PFrecHitsEB;
    ev.getByToken( PFrecHitCollection_EB_, PFrecHitsEB );
    const PFRecHitCollection* theBarrelEcalPFRecHits = PFrecHitsEB.product ();

//cout<<"I am inside usePFRecHitFlag_ condition"<<endl;
   
    for(reco::PFRecHitCollection::const_iterator PFitr = theBarrelEcalPFRecHits->begin() ; PFitr != theBarrelEcalPFRecHits->end() ; ++PFitr ){
    EBDetId ebid( PFitr -> detId() );
    GlobalPoint mycell = geometry -> getPosition(DetId(PFitr->detId()));
    double et = PFitr -> energy()*mycell.perp()/mycell.mag();

/*    if (PFitr -> energy() > maxPFRecHitEnergyEB ) {
            maxPFRecHitEnergyEB = PFitr -> energy() ;
            maxPFRecHitEtaEB = mycell.eta();
            maxPFRecHitPhiEB = mycell.phi();
            maxPFRecHitIEtaEB = ebid.ieta();
        }
*/
//    cout<<"et of PF"<<et<<endl;
//    if( et < ethrEB_ ) continue;
    h_PFRecHitEnergy->Fill(PFitr -> energy(), MyWeight);
    h_PFRecHitDepth->Fill(PFitr -> depth(), MyWeight);
    h_PFRecHitTime->Fill(PFitr -> time(), MyWeight);
//    h_PFRecHitpt2->Fill(PFitr -> pt2(), MyWeight);    

    h_PFRecHitEta->Fill(mycell.eta(), MyWeight);
    h_PFRecHitPhi->Fill(mycell.phi(), MyWeight);
    h_PFRecHitIEta->Fill(ebid.ieta(), MyWeight);

       }
    }

//-------------------------------end of PFrechits---------------------------


    //--- BASIC CLUSTERS --------------------------------------------------------------

    // ... barrel
    edm::Handle<reco::BasicClusterCollection> basicClusters_EB_h;
    ev.getByToken( basicClusterCollection_EB_, basicClusters_EB_h );
    const reco::BasicClusterCollection* theBarrelBasicClusters = basicClusters_EB_h.product () ;
    if ( ! basicClusters_EB_h.isValid() ) {
        std::cerr << "EcalValidation::analyze --> basicClusters_EB_h not found" << std::endl;
    }



    int nBCcleaned = 0;
    int nBCcleanedTkmatched = 0;

    for (reco::BasicClusterCollection::const_iterator itBC = theBarrelBasicClusters->begin();
            itBC != theBarrelBasicClusters->end(); ++itBC ) {

        //Get the associated RecHits
        const std::vector<std::pair<DetId,float> > & hits= itBC->hitsAndFractions();
        int nRechits_BC(0);
        for (std::vector<std::pair<DetId,float> > ::const_iterator rh = hits.begin(); rh!=hits.end(); ++rh) {

            EBRecHitCollection::const_iterator itrechit = theBarrelEcalRecHits->find((*rh).first);
            EBDetId ebid( itrechit -> id() );
            GlobalPoint mycell = geometry -> getPosition(DetId(itrechit->id()));
            double et = itrechit -> energy()*mycell.perp()/mycell.mag();
            if(et<4) continue;
            if (itrechit==theBarrelEcalRecHits->end()) continue;
            if(itrechit -> checkFlag(itrechit ->kOutOfTime)) continue;
            if(itrechit -> checkFlag(itrechit ->kWeird)) continue;
            if(debug) cout << "EB BC et (after all cuts )" << et << endl;
            h_basicClusters_recHits_EB_recoFlag -> Fill ( itrechit -> recoFlag(), MyWeight);
            h_basicClusters_recHits_recoFlag    -> Fill ( itrechit -> recoFlag(), MyWeight);
            ++nRechits_BC;
        }
        if(!nRechits_BC) continue;
        if(itBC->energy() < 6) continue;
        if(debug) cout << "num of rechits in BC " << hits.size() << "\t" << nRechits_BC << endl;
        //h_basicClusters_EB_nXtals -> Fill( (*itBC).hitsAndFractions().size(), MyWeight);
        h_basicClusters_EB_nXtals -> Fill( nRechits_BC, MyWeight);
        h_basicClusters_EB_energy -> Fill( itBC->energy(), MyWeight);
        h_basicClusters_eta       -> Fill( itBC->eta(), MyWeight);
        h_basicClusters_EB_eta    -> Fill( itBC->eta(), MyWeight);
        h_basicClusters_EB_phi    -> Fill( itBC->phi(), MyWeight);

        float E1 = EcalClusterTools::eMax   ( *itBC, theBarrelEcalRecHits);
//     float E9 = EcalClusterTools::e3x3   ( *itBC, theBarrelEcalRecHits, topology );
        float E4 = EcalClusterTools::eTop   ( *itBC, theBarrelEcalRecHits, topology )+
                   EcalClusterTools::eRight ( *itBC, theBarrelEcalRecHits, topology )+
                   EcalClusterTools::eBottom( *itBC, theBarrelEcalRecHits, topology )+
                   EcalClusterTools::eLeft  ( *itBC, theBarrelEcalRecHits, topology );

        if ( (1.-E4/E1) > 0.95 && E1 > 3.) continue;
        //h_basicClusters_EB_nXtals_cleaned -> Fill( (*itBC).hitsAndFractions().size(), MyWeight);
        h_basicClusters_EB_nXtals_cleaned -> Fill( nRechits_BC, MyWeight);
        h_basicClusters_EB_energy_cleaned -> Fill( itBC->energy(), MyWeight);
        nBCcleaned++;

        ///Do the Track-cluster matching for the cleaned ones

        float theBestDr = 99999.;
        for (edm::View<reco::Track>::const_iterator tkIt = TracksHandle->begin (); tkIt != TracksHandle->end (); ++tkIt )
        {
            ECALPositionCalculator posCalc;
            const math::XYZPoint vertex(BSPosition.x(),BSPosition.y(),tkIt->vz());
            const math::XYZVector trackMom =  tkIt->momentum();

            float phi= posCalc.ecalPhi(theMagField.product(),trackMom,vertex,tkIt -> charge());
            float deltaphi=fabs( phi - itBC->phi() );
            if(deltaphi>6.283185308) deltaphi -= 6.283185308;
            if(deltaphi>3.141592654) deltaphi = 6.283185308-deltaphi;

            float eta= posCalc.ecalEta(trackMom,vertex);
            float deltaeta=fabs( eta - itBC->eta() );
            if(deltaeta>6.283185308) deltaeta -= 6.283185308;
            if(deltaeta>3.141592654) deltaeta = 6.283185308-deltaeta;

            //compute dR squared
            float thisDr = deltaeta*deltaeta + deltaphi*deltaphi;
            if ( thisDr < theBestDr ) theBestDr = thisDr;
        }

        h_basicClusters_EB_dr_cleaned_tkmatched -> Fill ( theBestDr , MyWeight);
        //Matching CL-Tk if the dr < 0.2
        if ( theBestDr > 0.04 ) continue;

        h_basicClusters_EB_nXtals_cleaned_tkmatched -> Fill( nRechits_BC, MyWeight);
        h_basicClusters_EB_energy_cleaned_tkmatched -> Fill( itBC->energy(), MyWeight);
        nBCcleanedTkmatched++;

        h_basicClusters_EB_eta_tkmatched    -> Fill( itBC->eta(), MyWeight);
        h_basicClusters_EB_phi_tkmatched    -> Fill( itBC->phi(), MyWeight);
    }

    if(debug) cout << "nBcs " << nBCcleaned << "\t" << nBCcleanedTkmatched++  << endl;

    //h_basicClusters_EB_size         -> Fill( basicClusters_EB_h->size(), MyWeight);
    h_basicClusters_EB_size         -> Fill( nBCcleaned, MyWeight);
    h_basicClusters_EB_size_cleaned -> Fill( nBCcleaned , MyWeight);
    h_basicClusters_EB_size_cleaned_tkmatched -> Fill( nBCcleanedTkmatched , MyWeight);


    if(debug)  cout << "done with EB BC " << endl;
    // ... endcap
    edm::Handle<reco::BasicClusterCollection> basicClusters_EE_h;
    ev.getByToken( basicClusterCollection_EE_, basicClusters_EE_h );
    if ( ! basicClusters_EE_h.isValid() ) {
        std::cerr << "EcalValidation::analyze --> basicClusters_EE_h not found" << std::endl;
    }


    int nBasicClustersEEP = 0;
    int nBasicClustersEEM = 0;
    int nBasicClustersEEPTkmatched = 0;
    int nBasicClustersEEMTkmatched = 0;
    int nCrystal_BC_EE(0);
    for (unsigned int icl = 0; icl < basicClusters_EE_h->size(); ++icl) {

        if( (*basicClusters_EE_h)[icl].energy() < 8) continue;
        //Get the associated RecHits
        const std::vector<std::pair<DetId,float> > & hits= (*basicClusters_EE_h)[icl].hitsAndFractions();
        for (std::vector<std::pair<DetId,float> > ::const_iterator rh = hits.begin(); rh!=hits.end(); ++rh) {

            EcalRecHitCollection::const_iterator itrechit = theEndcapEcalRecHits->find((*rh).first);
            EEDetId eeid( itrechit -> id() );
            GlobalPoint mycell = geometry->getPosition(itrechit->detid());
            double et = itrechit -> energy()*mycell.perp()/mycell.mag();
            if(et < 8) continue;
            if(itrechit->checkFlag(itrechit->kOutOfTime)) continue;
            if(itrechit->checkFlag(itrechit->kWeird)) continue;
//if(debug) cout << "ee hits " << itrechit->checkFlag(itrechit->kOutOfTime) << "\t" << itrechit->checkFlag(itrechit->kWeird) << endl;
            if(debug) cout << "et of EE hits in BC " << et << endl;
            if (itrechit==theEndcapEcalRecHits->end()) continue;
            h_basicClusters_recHits_EE_recoFlag -> Fill ( itrechit -> recoFlag(), MyWeight);
            h_basicClusters_recHits_recoFlag    -> Fill ( itrechit -> recoFlag(), MyWeight);
            ++nCrystal_BC_EE;
        }
        if(!nCrystal_BC_EE) continue;

        h_basicClusters_eta       -> Fill( (*basicClusters_EE_h)[icl].eta(), MyWeight);
        h_basicClusters_EE_eta    -> Fill( (*basicClusters_EE_h)[icl].eta(), MyWeight);
        h_basicClusters_EE_phi    -> Fill( (*basicClusters_EE_h)[icl].phi(), MyWeight);

        if ((*basicClusters_EE_h)[icl].z() > 0) {
            //h_basicClusters_EEP_nXtals -> Fill( (*basicClusters_EE_h)[icl].hitsAndFractions().size(), MyWeight);
            h_basicClusters_EEP_nXtals -> Fill( nCrystal_BC_EE, MyWeight);
            h_basicClusters_EEP_energy -> Fill( (*basicClusters_EE_h)[icl].energy(), MyWeight);
            nBasicClustersEEP++;


            float theBestDr = 99999.;
            for (edm::View<reco::Track>::const_iterator tkIt = TracksHandle->begin (); tkIt != TracksHandle->end (); ++tkIt )
            {
                ECALPositionCalculator posCalc;
                const math::XYZPoint vertex(BSPosition.x(),BSPosition.y(),tkIt->vz());
                const math::XYZVector trackMom =  tkIt->momentum();

                float phi= posCalc.ecalPhi(theMagField.product(),trackMom,vertex,tkIt -> charge());
                float deltaphi=fabs( phi - (*basicClusters_EE_h)[icl].phi());
                if(deltaphi>6.283185308) deltaphi -= 6.283185308;
                if(deltaphi>3.141592654) deltaphi = 6.283185308-deltaphi;

                float eta= posCalc.ecalEta(trackMom,vertex);
                float deltaeta=fabs( eta - (*basicClusters_EE_h)[icl].eta());
                if(deltaeta>6.283185308) deltaeta -= 6.283185308;
                if(deltaeta>3.141592654) deltaeta = 6.283185308-deltaeta;

                //compute dR squared
                float thisDr = deltaeta*deltaeta + deltaphi*deltaphi;
                if ( thisDr < theBestDr ) theBestDr = thisDr;
            }
            h_basicClusters_EEP_dr_tkmatched     -> Fill( theBestDr , MyWeight);
            //Matching CL-Tk if the dr < 0.2
            if ( theBestDr > 0.04 ) continue;

            h_basicClusters_EEP_nXtals_tkmatched -> Fill( (*basicClusters_EE_h)[icl].hitsAndFractions().size(), MyWeight);
            h_basicClusters_EEP_energy_tkmatched -> Fill( (*basicClusters_EE_h)[icl].energy(), MyWeight);
            nBasicClustersEEPTkmatched++;

            h_basicClusters_eta_tkmatched       -> Fill(  (*basicClusters_EE_h)[icl].eta(), MyWeight);
            h_basicClusters_EE_eta_tkmatched    -> Fill(  (*basicClusters_EE_h)[icl].eta(), MyWeight);
            h_basicClusters_EE_phi_tkmatched    -> Fill(  (*basicClusters_EE_h)[icl].phi(), MyWeight);

        }

        if ((*basicClusters_EE_h)[icl].z() < 0) {
//      h_basicClusters_EEM_nXtals -> Fill( (*basicClusters_EE_h)[icl].hitsAndFractions().size(), MyWeight);
            h_basicClusters_EEM_nXtals -> Fill( nCrystal_BC_EE, MyWeight);
            h_basicClusters_EEM_energy -> Fill( (*basicClusters_EE_h)[icl].energy(), MyWeight);
            nBasicClustersEEM++;

            float theBestDr = 99999.;
            for (edm::View<reco::Track>::const_iterator tkIt = TracksHandle->begin (); tkIt != TracksHandle->end (); ++tkIt )
            {
                ECALPositionCalculator posCalc;
                const math::XYZPoint vertex(BSPosition.x(),BSPosition.y(),tkIt->vz());
                const math::XYZVector trackMom =  tkIt->momentum();

                float phi= posCalc.ecalPhi(theMagField.product(),trackMom,vertex,tkIt -> charge());
                float deltaphi=fabs( phi - (*basicClusters_EE_h)[icl].phi());
                if(deltaphi>6.283185308) deltaphi -= 6.283185308;
                if(deltaphi>3.141592654) deltaphi = 6.283185308-deltaphi;

                float eta= posCalc.ecalEta(trackMom,vertex);
                float deltaeta=fabs( eta - (*basicClusters_EE_h)[icl].eta());
                if(deltaeta>6.283185308) deltaeta -= 6.283185308;
                if(deltaeta>3.141592654) deltaeta = 6.283185308-deltaeta;

                //compute dR squared
                float thisDr = deltaeta*deltaeta + deltaphi*deltaphi;
                if ( thisDr < theBestDr ) theBestDr = thisDr;
            }
            h_basicClusters_EEM_dr_tkmatched     -> Fill( theBestDr , MyWeight);
            //Matching CL-Tk if the dr < 0.2
            if ( theBestDr > 0.04 ) continue;

            h_basicClusters_EEM_nXtals_tkmatched -> Fill( (*basicClusters_EE_h)[icl].hitsAndFractions().size(), MyWeight);
            h_basicClusters_EEM_energy_tkmatched -> Fill( (*basicClusters_EE_h)[icl].energy(), MyWeight);
            nBasicClustersEEMTkmatched++;

            h_basicClusters_eta_tkmatched       -> Fill(  (*basicClusters_EE_h)[icl].eta(), MyWeight);
            h_basicClusters_EE_eta_tkmatched    -> Fill(  (*basicClusters_EE_h)[icl].eta(), MyWeight);
            h_basicClusters_EE_phi_tkmatched    -> Fill(  (*basicClusters_EE_h)[icl].phi(), MyWeight);

        }
    }

    h_basicClusters_EEP_size->Fill( nBasicClustersEEP , MyWeight);
    h_basicClusters_EEM_size->Fill( nBasicClustersEEM , MyWeight);
    h_basicClusters_EEP_size_tkmatched->Fill( nBasicClustersEEPTkmatched , MyWeight);
    h_basicClusters_EEM_size_tkmatched->Fill( nBasicClustersEEMTkmatched , MyWeight);

    if(debug)  cout << "done with EE BC " << endl;
    // Super Clusters
    // ... barrel
    edm::Handle<reco::SuperClusterCollection> superClusters_EB_h;
    ev.getByToken( superClusterCollection_EB_, superClusters_EB_h );
    const reco::SuperClusterCollection* theBarrelSuperClusters = superClusters_EB_h.product () ;
    if ( ! superClusters_EB_h.isValid() ) {
        std::cerr << "EcalValidation::analyze --> superClusters_EB_h not found" << std::endl;
    }


    int nSCcleaned = 0;
    int nSC = 0;

    for (reco::SuperClusterCollection::const_iterator itSC = theBarrelSuperClusters->begin();
            itSC != theBarrelSuperClusters->end(); ++itSC ) {
        double scEt = itSC -> energy() * sin(2.*atan( exp(- itSC->position().eta() )));
        double scRawEt = itSC -> rawEnergy() * sin(2.*atan( exp(- itSC->position().eta() )));
        if (scEt < scEtThrEB_ ) continue;
        int nBC_SC(0);
        ///get the BC within the SC
        for(reco::CaloCluster_iterator bcIt = (*itSC).clustersBegin(); bcIt!=(*itSC).clustersEnd(); bcIt++) {
//cout << "basic cluster energy " <<(*bcIt)->energy()  << endl;
            if((*bcIt)->energy() < 6)continue;
            ++nBC_SC;
        }
//cout << "# BC in SC " << nBC_SC << "\t initial # " <<  (*itSC).clustersSize() << endl;


        //nSC++;

        float phi = itSC -> phi();
        float eta = itSC -> eta();

        if ((phi >= -3.01942 && phi <= -3.00197 && eta >= -1.36  && eta <= -1.32) ||
                (phi >= -3.00197 && phi <= -2.98451 && eta >= -1.32  && eta <= -1.28) ||
                (phi >= 0.994838 && phi <= 1.01229  && eta >= 0.16   && eta <= 0.2  ) ||
                (phi >= 1.09956  && phi <= 1.11701  && eta >= 1.16   && eta <= 1.2  )    ) {}
        else h_superClusters_eta_scCut -> Fill( itSC -> eta(), MyWeight );


        float E1 = EcalClusterTools::eMax   ( *itSC, theBarrelEcalRecHits);
        float E4 = EcalClusterTools::eTop   ( *itSC, theBarrelEcalRecHits, topology )+
                   EcalClusterTools::eRight ( *itSC, theBarrelEcalRecHits, topology )+
                   EcalClusterTools::eBottom( *itSC, theBarrelEcalRecHits, topology )+
                   EcalClusterTools::eLeft  ( *itSC, theBarrelEcalRecHits, topology );

//loop over all the rechits

        int nHits_SCEB(0);
        const std::vector<std::pair<DetId,float> > & hits= (*itSC).hitsAndFractions();
//if(debug) cout << "SC # of rechits " << hits.size() << endl;
        for (std::vector<std::pair<DetId,float> > ::const_iterator rh = hits.begin(); rh!=hits.end(); ++rh) {

            EBRecHitCollection::const_iterator itrechit = theBarrelEcalRecHits->find((*rh).first);
            EBDetId eeid( itrechit -> id() );
            GlobalPoint mycell = geometry->getPosition(itrechit->detid());
            double et = itrechit -> energy()*mycell.perp()/mycell.mag();
            if(et < 4) continue;
            if(itrechit->checkFlag(itrechit->kOutOfTime)) continue;
            if(itrechit->checkFlag(itrechit->kWeird)) continue;
            if(debug) cout << "et of EB hits " << et << endl;
//if(debug) cout << "eb hits " << itrechit->checkFlag(itrechit->kOutOfTime) << "\t" << itrechit->checkFlag(itrechit->kWeird) << endl;
            if (itrechit==theBarrelEcalRecHits->end()) continue;
            ++nHits_SCEB;
        }
        if(!nHits_SCEB) continue;
        if(debug) cout << "nhits in SC "  << "\t" << nHits_SCEB << endl;



        //Now get the seed:
        EBDetId theSeedIdEB = EcalClusterTools::getMaximum( (*itSC).hitsAndFractions(), theBarrelEcalRecHits ).first;
        EcalRecHitCollection::const_iterator theSeedEB = theBarrelEcalRecHits->find (theSeedIdEB) ;
        GlobalPoint mycell = geometry->getPosition(theSeedEB->detid());
        double et = theSeedEB -> energy()*mycell.perp()/mycell.mag();
        if(et < 4.0 ) continue;
        if(theSeedEB->checkFlag(theSeedEB->kOutOfTime)) continue;
        if(theSeedEB->checkFlag(theSeedEB->kWeird)) continue;
        if(debug) cout << "SC et of EB hits " << et << endl;
//if(debug) cout << "SC ee hits " << theSeedEB->checkFlag(theSeedEB->kOutOfTime) << "\t" << theSeedEB->checkFlag(theSeedEB->kWeird) << endl;
        if ( (1.- E4/E1) > 0.95  && E1 > 3. ) continue;

        nSC++;
        //h_superClusters_EB_nXtals -> Fill( (*itSC).hitsAndFractions().size(), MyWeight);
        h_superClusters_EB_nXtals -> Fill( nHits_SCEB, MyWeight);
        //h_superClusters_EB_nBC    -> Fill( itSC -> clustersSize(), MyWeight);
        h_superClusters_EB_nBC    -> Fill( nBC_SC, MyWeight);
        h_superClusters_EB_energy -> Fill( itSC -> energy(), MyWeight);
        h_superClusters_eta       -> Fill( itSC -> eta(), MyWeight);
        h_superClusters_EB_eta    -> Fill( itSC -> eta(), MyWeight);
        h_superClusters_EB_phi    -> Fill( itSC -> phi(), MyWeight);
        h_superClusters_occupancyPhiEta -> Fill(itSC -> phi(),itSC -> eta(), MyWeight);
        h_superClusters_EB_occupancy ->Fill(theSeedIdEB.iphi(),theSeedIdEB.ieta(), MyWeight);

        EcalPedestalsMap::const_iterator itped = pedestals->getMap().find( theSeedIdEB );
        float width = (*itped).rms(gain);

        if ( E1 > 3. ) h_superClusters_EB_E1oE4  -> Fill( 1.- E4/E1, MyWeight);

        //Fill plots for SCL seed
        h2_superClusters_EB_seedTimeVsEnergy -> Fill ( theSeedEB->energy(), theSeedEB->time(), MyWeight);

        if ( (1.- E4/E1) > 0.95  && E1 > 3. ) continue;

        h_superClusters_EB_nXtals_cleaned -> Fill( (*itSC).hitsAndFractions().size(), MyWeight);
        //h_superClusters_EB_nBC_cleaned    -> Fill( itSC -> clustersSize(), MyWeight);
        h_superClusters_EB_nBC_cleaned    -> Fill( nBC_SC, MyWeight);
        h_superClusters_EB_energy_cleaned -> Fill( itSC -> energy(), MyWeight);
        h_superClusters_EB_rawEnergy_cleaned -> Fill( itSC -> rawEnergy(), MyWeight);
        h_superClusters_EB_rawEt_cleaned -> Fill( scRawEt , MyWeight);

        nSCcleaned++;

    }
    if(debug)
        cout << "# of SC ucleaned/cleaned EB " << nSC << "\t" << nSCcleaned << endl;
    h_superClusters_EB_size         -> Fill( nSC , MyWeight);
    h_superClusters_EB_size_cleaned -> Fill( nSCcleaned , MyWeight);

    if(debug)  cout << "done with SC barrel " << endl;


    // ... endcap
    edm::Handle<reco::SuperClusterCollection> superClusters_EE_h;
    ev.getByToken( superClusterCollection_EE_, superClusters_EE_h );
    const reco::SuperClusterCollection* theEndcapSuperClusters = superClusters_EE_h.product () ;
    if ( ! superClusters_EE_h.isValid() ) {
        std::cerr << "EcalValidation::analyze --> superClusters_EE_h not found" << std::endl;
    }

    int nSuperClustersEEP = 0;
    int nSuperClustersEEM = 0;
    int nHits_SCEE(0);
    for (reco::SuperClusterCollection::const_iterator itSC = theEndcapSuperClusters->begin();
            itSC != theEndcapSuperClusters->end(); ++itSC ) {

        double scEt = itSC -> energy() * sin(2.*atan( exp(- itSC->position().eta() )));
        double scRawEt = itSC -> rawEnergy() * sin(2.*atan( exp(- itSC->position().eta() )));

        if (scEt < scEtThrEE_ ) continue;
        float phi = itSC -> phi();
        float eta = itSC -> eta();

        if ((phi >= -3.01942 && phi <= -3.00197 && eta >= -1.36  && eta <= -1.32) ||
                (phi >= -3.00197 && phi <= -2.98451 && eta >= -1.32  && eta <= -1.28) ||
                (phi >= 0.994838 && phi <= 1.01229  && eta >= 0.16   && eta <= 0.2  ) ||
                (phi >= 1.09956  && phi <= 1.11701  && eta >= 1.16   && eta <= 1.2  )    ) {}
        else h_superClusters_eta_scCut -> Fill( itSC -> eta(), MyWeight);

        //Now get the seed:
        EEDetId theSeedIdEE = EcalClusterTools::getMaximum( (*itSC).hitsAndFractions(), theEndcapEcalRecHits ).first;
        EcalRecHitCollection::const_iterator theSeedEE = theEndcapEcalRecHits->find (theSeedIdEE) ;

//we should get a loop over the number of hits to get the selected hits.

        const std::vector<std::pair<DetId,float> > & hits= (*itSC).hitsAndFractions();
        for (std::vector<std::pair<DetId,float> > ::const_iterator rh = hits.begin(); rh!=hits.end(); ++rh) {

            EcalRecHitCollection::const_iterator itrechit = theEndcapEcalRecHits->find((*rh).first);
            EEDetId eeid( itrechit -> id() );
            GlobalPoint mycell = geometry->getPosition(itrechit->detid());
            double et = itrechit -> energy()*mycell.perp()/mycell.mag();
            if(et < 8) continue;
            if(itrechit->checkFlag(itrechit->kOutOfTime)) continue;
            if(itrechit->checkFlag(itrechit->kWeird)) continue;
            if(debug) cout << "et of EE hits " << et << endl;
//if(debug) cout << "ee hits " << itrechit->checkFlag(itrechit->kOutOfTime) << "\t" << itrechit->checkFlag(itrechit->kWeird) << endl;
            if (itrechit==theEndcapEcalRecHits->end()) continue;
            ++nHits_SCEE;
        }
        if(!nHits_SCEE) continue;
        if(debug) cout << "nhits in SC in EE "  << "\t" << nHits_SCEE << endl;

        int nBC_SC(0);
        ///get the BC within the SC
        for(reco::CaloCluster_iterator bcIt = (*itSC).clustersBegin(); bcIt!=(*itSC).clustersEnd(); bcIt++) {
//cout << "basic cluster energy " <<(*bcIt)->energy()  << endl;
            if((*bcIt)->energy() < 6)continue;
            ++nBC_SC;
        }
        if(!nBC_SC) continue;
        if(theSeedIdEE.zside() > 0) h_superClusters_EEP_occupancy ->Fill(theSeedIdEE.ix()-0.5,theSeedIdEE.iy()-0.5, MyWeight);
        if(theSeedIdEE.zside() < 0) h_superClusters_EEM_occupancy ->Fill(theSeedIdEE.ix()-0.5,theSeedIdEE.iy()-0.5, MyWeight);

        EcalPedestalsMap::const_iterator itped = pedestals->getMap().find( theSeedIdEE );
        float width = (*itped).rms(gain);

        GlobalPoint mycell = geometry->getPosition(theSeedEE->detid());
        double et = theSeedEE -> energy()*mycell.perp()/mycell.mag();
        if(et < 8.0 ) continue;
        if(theSeedEE->checkFlag(theSeedEE->kOutOfTime)) continue;
        if(theSeedEE->checkFlag(theSeedEE->kWeird)) continue;
        if(debug) cout << "SC et of Seed hits " << et << endl;
//if(debug) cout << "SC ee hits " << theSeedEE->checkFlag(theSeedEE->kOutOfTime) << "\t" << theSeedEE->checkFlag(theSeedEE->kWeird) << endl;

        h_superClusters_eta       -> Fill( itSC -> eta(), MyWeight);
        h_superClusters_EE_eta    -> Fill( itSC -> eta(), MyWeight);
        h_superClusters_EE_phi    -> Fill( itSC -> phi(), MyWeight);
        h_superClusters_occupancyPhiEta -> Fill(itSC -> phi(),itSC -> eta(), MyWeight);


        float E1 = EcalClusterTools::eMax( *itSC, theEndcapEcalRecHits);
        float E4 = EcalClusterTools::eTop( *itSC, theEndcapEcalRecHits, topology)+
                   EcalClusterTools::eRight( *itSC, theEndcapEcalRecHits, topology)+
                   EcalClusterTools::eBottom( *itSC, theEndcapEcalRecHits, topology)+
                   EcalClusterTools::eLeft( *itSC, theEndcapEcalRecHits, topology);

        //Fill plots for SCL seed
        h2_superClusters_EE_seedTimeVsEnergy -> Fill ( theSeedEE->energy(), theSeedEE->time(), MyWeight);


        if  ( itSC -> z() > 0 ) {
            //h_superClusters_EEP_nXtals -> Fill( (*itSC).hitsAndFractions().size(), MyWeight);
            h_superClusters_EEP_nXtals -> Fill( nHits_SCEE, MyWeight);
            //h_superClusters_EEP_nBC    -> Fill( itSC -> clustersSize(), MyWeight);
            h_superClusters_EEP_nBC    -> Fill( nBC_SC, MyWeight);
            h_superClusters_EEP_energy -> Fill( itSC -> energy(), MyWeight);
            h_superClusters_EEP_rawEnergy -> Fill( itSC -> rawEnergy(), MyWeight);
            h_superClusters_EEP_rawEt -> Fill( scRawEt , MyWeight);
            if ( E1 > 3. ) h_superClusters_EEP_E1oE4  -> Fill( 1.- E4/E1, MyWeight);
            nSuperClustersEEP++;
        }

        if  ( itSC -> z() < 0 ) {
            //h_superClusters_EEM_nXtals -> Fill( (*itSC).hitsAndFractions().size(), MyWeight);
            h_superClusters_EEM_nXtals -> Fill( nHits_SCEE, MyWeight);
            //h_superClusters_EEM_nBC    -> Fill( itSC -> clustersSize(), MyWeight);
            h_superClusters_EEM_nBC    -> Fill( nBC_SC, MyWeight);
            h_superClusters_EEM_energy -> Fill( itSC -> energy(), MyWeight);
            h_superClusters_EEM_rawEnergy -> Fill( itSC -> rawEnergy(), MyWeight);
            h_superClusters_EEM_rawEt -> Fill( scRawEt , MyWeight);
            if ( E1 > 3. ) h_superClusters_EEM_E1oE4  -> Fill( 1.- E4/E1 , MyWeight);
            nSuperClustersEEM++;
        }
    }
    if(debug) {
        cout << "no. of SCs EEP/EEM and nhits in SC "  << nSuperClustersEEP << "\t" << nSuperClustersEEM << "\t" << nHits_SCEE << endl;
    }
    if(debug)  cout << "done with SC EC " << endl;

    h_superClusters_EEP_size->Fill( nSuperClustersEEP , MyWeight);
    h_superClusters_EEM_size->Fill( nSuperClustersEEM , MyWeight);

    //--------------------------------------------------------

    // ES clusters in X plane
    Handle<PreshowerClusterCollection> esClustersX;
    ev.getByToken( esClusterCollectionX_, esClustersX);
    const PreshowerClusterCollection *ESclustersX = esClustersX.product();

    // ES clusters in Y plane
    Handle<PreshowerClusterCollection> esClustersY;
    ev.getByToken( esClusterCollectionY_, esClustersY);
    const PreshowerClusterCollection *ESclustersY = esClustersY.product();

    // Do the ES-BasicCluster matching
    for (unsigned int icl = 0; icl < basicClusters_EE_h->size(); ++icl) {

        const CaloCluster ecalBasicCluster = (*basicClusters_EE_h)[icl];

        bool isMatchedESx = false;
        bool isMatchedESy = false;
        for (PreshowerClusterCollection::const_iterator iESClus = ESclustersX->begin(); iESClus != ESclustersX->end();
                ++iESClus) {
            const CaloClusterPtr preshBasicCluster = iESClus->basicCluster();
            if (*preshBasicCluster == ecalBasicCluster) isMatchedESx = true;
        }

        for (PreshowerClusterCollection::const_iterator iESClus = ESclustersY->begin(); iESClus != ESclustersY->end();
                ++iESClus) {
            const CaloClusterPtr preshBasicCluster = iESClus->basicCluster();
            if (*preshBasicCluster == ecalBasicCluster) isMatchedESy = true;
        }

        //Fot the matching require at least a matching with 1 ES plane
        if ( !isMatchedESx && !isMatchedESy ) continue;

        if ((*basicClusters_EE_h)[icl].z() > 0) {

            h_basicClusters_EEP_occupancy_esmatched -> Fill (  ecalBasicCluster.eta(),ecalBasicCluster.phi(), MyWeight);
            h_basicClusters_EEP_eta_esmatched -> Fill (  ecalBasicCluster.eta(), MyWeight);
            h_basicClusters_EEP_phi_esmatched -> Fill (  ecalBasicCluster.phi(), MyWeight);

        }

        if ((*basicClusters_EE_h)[icl].z() < 0) {

            h_basicClusters_EEM_occupancy_esmatched -> Fill (  ecalBasicCluster.eta(),ecalBasicCluster.phi(), MyWeight);
            h_basicClusters_EEM_eta_esmatched -> Fill (  ecalBasicCluster.eta(), MyWeight);
            h_basicClusters_EEM_phi_esmatched -> Fill (  ecalBasicCluster.phi(), MyWeight);

        }


    }


    // loop over all super clusters
    for (reco::SuperClusterCollection::const_iterator itSC = theEndcapSuperClusters->begin();
            itSC != theEndcapSuperClusters->end(); ++itSC ) {

        if ( fabs(itSC->eta()) < 1.65 || fabs(itSC->eta()) > 2.6 ) continue;

        // Loop over all ECAL Basic clusters in the supercluster
        for (CaloCluster_iterator ecalBasicCluster = itSC->clustersBegin(); ecalBasicCluster!= itSC->clustersEnd();
                ecalBasicCluster++) {
            const CaloClusterPtr ecalBasicClusterPtr = *(ecalBasicCluster);

            float ESenergyPlane1 = -999.;
            float ESenergyPlane2 = -999.;

            for (PreshowerClusterCollection::const_iterator iESClus = ESclustersX->begin(); iESClus != ESclustersX->end();
                    ++iESClus) {
                const CaloClusterPtr preshBasicCluster = iESClus->basicCluster();
                const PreshowerCluster *esCluster = &*iESClus;
                if (preshBasicCluster == ecalBasicClusterPtr) {
                    ESenergyPlane1 = esCluster->energy();
                    h_esClusters_energy_plane1 ->Fill(esCluster->energy(), MyWeight);
                }
            }  // end of x loop

            for (PreshowerClusterCollection::const_iterator iESClus = ESclustersY->begin(); iESClus != ESclustersY->end();
                    ++iESClus) {
                const CaloClusterPtr preshBasicCluster = iESClus->basicCluster();
                const PreshowerCluster *esCluster = &*iESClus;
                if (preshBasicCluster == ecalBasicClusterPtr) {
                    ESenergyPlane2 = esCluster->energy();
                    h_esClusters_energy_plane2 -> Fill(esCluster->energy(), MyWeight);
                }
            } // end of y loop

            if ( ESenergyPlane1 != -999. && ESenergyPlane2 != -999. )
                h_esClusters_energy_ratio -> Fill(ESenergyPlane1/ESenergyPlane2, MyWeight);


        } // end loop over all basic clusters in the supercluster
    }// end loop over superclusters


    // ---------- Do histos for pi0 peak

    doPi0Barrel(geometry, topology, recHitsEB);
    doPi0Endcap(geometry, topology, recHitsEE);

    /// Get Jets and draw EMF maps
    edm::Handle<reco::CaloJetCollection> JetHandle ;
    ev.getByToken (jets_,JetHandle);

    for(unsigned int i=0; i<JetHandle->size(); ++i)
    {


        // ... Barrel
        if ( (*JetHandle)[i].p4().eta() < 1.4442 ) {

            h_Jets_EB_emf        -> Fill ( (*JetHandle)[i].p4().eta(), (*JetHandle)[i].p4().phi(), (*JetHandle)[i].emEnergyFraction(), MyWeight);
            h_Jets_EB_emf_eta    -> Fill ( (*JetHandle)[i].p4().eta(), (*JetHandle)[i].emEnergyFraction(), MyWeight);
            h_Jets_EB_emf_phi    -> Fill ( (*JetHandle)[i].p4().phi(), (*JetHandle)[i].emEnergyFraction(), MyWeight);

        }

        // ... Endcaps
        if ( (*JetHandle)[i].p4().eta() > 1.566 ) {

            h_Jets_EEP_emf        -> Fill ( (*JetHandle)[i].p4().eta(), (*JetHandle)[i].p4().phi(), (*JetHandle)[i].emEnergyFraction(), MyWeight);
            h_Jets_EEP_emf_eta    -> Fill ( (*JetHandle)[i].p4().eta(), (*JetHandle)[i].emEnergyFraction(), MyWeight);
            h_Jets_EEP_emf_phi    -> Fill ( (*JetHandle)[i].p4().phi(), (*JetHandle)[i].emEnergyFraction(), MyWeight);

        }

        if ( (*JetHandle)[i].p4().eta() < -1.566 ) {

            h_Jets_EEM_emf        -> Fill ( (*JetHandle)[i].p4().eta(), (*JetHandle)[i].p4().phi(), (*JetHandle)[i].emEnergyFraction(), MyWeight);
            h_Jets_EEM_emf_eta    -> Fill ( (*JetHandle)[i].p4().eta(), (*JetHandle)[i].emEnergyFraction(), MyWeight);
            h_Jets_EEM_emf_phi    -> Fill ( (*JetHandle)[i].p4().phi(), (*JetHandle)[i].emEnergyFraction(), MyWeight);

        }

    }// loop on jets
    tree->Fill();
}


// ------------ method called once each job just before starting event loop  ------------
void
EcalValidationAOD::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void
EcalValidationAOD::endJob() {


    h_numberOfEvents ->Fill(0.,naiveId_);

    /// ---------- Compute and Fill RecHits occupancy deviations

    //Barrel
    std::cout << "\nNUMBER OF EVENTS: " << naiveId_ << std::endl;
    for(int ii = 0; ii < naiveId_; ii++)
        h_recHits_EB_sumEt->Fill(sumRecHitEtEB[ii]);

    //Intialize RMS, Mean and #AliveChannels for each eta slice.
    float recHits_EB_mean[173];
    float recHits_EB_RMS [173];
    int alive_channels_EB[173];
    for ( int i = 1; i < 173; i++ )
    {
        recHits_EB_mean[i] = 0;
        recHits_EB_RMS[i] = 0;
        alive_channels_EB[i] = 0;
    }

    for ( int iPhi = 1; iPhi < 361; iPhi++ )
        for ( int iEta = 1; iEta < 173; iEta++ )
        {
            if ( h_recHits_EB_occupancy -> GetBinContent(iPhi, iEta) == 0 ) continue ;
            alive_channels_EB[iEta] ++ ;
            recHits_EB_mean[iEta] += h_recHits_EB_occupancy -> GetBinContent(iPhi, iEta) ;
            recHits_EB_RMS[iEta]  += (h_recHits_EB_occupancy -> GetBinContent(iPhi, iEta)) * (h_recHits_EB_occupancy -> GetBinContent(iPhi, iEta)) ;
        }

    //Now compute RMS using the formula RMS = sqrt (E[x^2] - E[x]^2)
    for ( int i = 1; i < 173; i++ )
    {
        recHits_EB_mean[i] = recHits_EB_mean[i]/alive_channels_EB[i] ;
        recHits_EB_RMS[i]  = recHits_EB_RMS[i]/alive_channels_EB[i] ;
        recHits_EB_RMS[i]  = sqrt(recHits_EB_RMS[i] - recHits_EB_mean[i]*recHits_EB_mean[i]) ;
    }

    for ( int iPhi = 1; iPhi < 361; iPhi++ )
        for ( int iEta = 1; iEta < 173; iEta++ )
        {
            if ( h_recHits_EB_occupancy -> GetBinContent(iPhi, iEta) == 0 || recHits_EB_RMS[iEta] == 0 ) continue ;

            h_recHits_EB_deviation -> SetBinContent(iPhi, iEta,
                                                    ( fabs(h_recHits_EB_occupancy -> GetBinContent(iPhi, iEta) - recHits_EB_mean[iEta]))/recHits_EB_RMS[iEta]
                                                   );

        }


    //EndCaps

    //EE+

    for(int ii = 0; ii < naiveId_; ii++) {
        h_recHits_EEP_sumEt->Fill(sumRecHitEtEEP[ii]);
        h_recHits_EEP_sumEtCut->Fill(sumRecHitEtCutEEP[ii]);
    }

    //Intialize RMS, Mean and #AliveChannels for each eta slice.
    //Eta slices are not corresponding to DeltaEta = cost, but are defined as Delta[sqrt(x^2 + y^2)] = 1
    float recHits_EEP_mean[100];
    float recHits_EEP_RMS [100];
    int alive_channels_EEP[100];
    for ( int i = 0; i < 100; i++ )
    {
        recHits_EEP_mean[i] = 0;
        recHits_EEP_RMS[i] = 0;
        alive_channels_EEP[i] = 0;
    }

    for ( int iX = 1; iX < 101; iX++ )
        for ( int iY = 1; iY < 101; iY++ )
        {
            if ( h_recHits_EEP_occupancy -> GetBinContent(iX, iY) == 0 ) continue ;
            int iRing = (int) sqrt( (float) (iX-51)*(iX-51) + (float)(iY-51)*(iY-51) ) ;
            alive_channels_EEP[iRing] ++ ;
            recHits_EEP_mean[iRing] += h_recHits_EEP_occupancy -> GetBinContent(iX, iY) ;
            recHits_EEP_RMS[iRing]  += (h_recHits_EEP_occupancy -> GetBinContent(iX, iY)) * (h_recHits_EEP_occupancy -> GetBinContent(iX, iY)) ;
        }

    //Now compute RMS using the formula RMS = sqrt (E[x^2] - E[x]^2)
    for ( int i = 0; i < 100; i++ )
    {
        if ( alive_channels_EEP[i] == 0 ) continue ;
        recHits_EEP_mean[i] = recHits_EEP_mean[i]/alive_channels_EEP[i] ;
        recHits_EEP_RMS[i]  = recHits_EEP_RMS[i]/alive_channels_EEP[i] ;
        recHits_EEP_RMS[i]  = sqrt(recHits_EEP_RMS[i] - recHits_EEP_mean[i]*recHits_EEP_mean[i]) ;
    }

    for ( int iX = 1; iX < 101; iX++ )
        for ( int iY = 1; iY < 101; iY++ )
        {
            int iRing = (int) sqrt( (float) (iX-51)*(iX-51) + (float)(iY-51)*(iY-51) ) ;
            if ( h_recHits_EEP_occupancy -> GetBinContent(iX, iY) == 0 || recHits_EEP_RMS[iRing] == 0 ) continue ;

            h_recHits_EEP_deviation -> SetBinContent(iX, iY,
                    ( fabs(h_recHits_EEP_occupancy -> GetBinContent(iX, iY) - recHits_EEP_mean[iRing]))/recHits_EEP_RMS[iRing]
                                                    );

        }

    //EE-

    for(int ii = 0; ii < naiveId_; ii++) {
        h_recHits_EEM_sumEt->Fill(sumRecHitEtEEM[ii]);
        h_recHits_EEM_sumEtCut->Fill(sumRecHitEtCutEEM[ii]);
    }

    //Intialize RMS, Mean and #AliveChannels for each eta slice.
    //Eta slices are not corresponding to DeltaEta = cost, but are defined as Delta[sqrt(x^2 + y^2)] = 1
    float recHits_EEM_mean[100];
    float recHits_EEM_RMS [100];
    int alive_channels_EEM[100];
    for ( int i = 0; i < 100; i++ )
    {
        recHits_EEM_mean[i] = 0;
        recHits_EEM_RMS[i] = 0;
        alive_channels_EEM[i] = 0;
    }

    for ( int iX = 1; iX < 101; iX++ )
        for ( int iY = 1; iY < 101; iY++ )
        {
            if ( h_recHits_EEM_occupancy -> GetBinContent(iX, iY) == 0 ) continue ;
            int iRing = (int) sqrt( (float) (iX-51)*(iX-51) + (float)(iY-51)*(iY-51) ) ;
            alive_channels_EEM[iRing] ++ ;
            recHits_EEM_mean[iRing] += h_recHits_EEM_occupancy -> GetBinContent(iX, iY) ;
            recHits_EEM_RMS[iRing]  += (h_recHits_EEM_occupancy -> GetBinContent(iX, iY)) * (h_recHits_EEM_occupancy -> GetBinContent(iX, iY)) ;
        }

    //Now compute RMS using the formula RMS = sqrt (E[x^2] - E[x]^2)
    for ( int i = 0; i < 100; i++ )
    {
        if ( alive_channels_EEM[i] == 0 ) continue ;
        recHits_EEM_mean[i] = recHits_EEM_mean[i]/alive_channels_EEM[i] ;
        recHits_EEM_RMS[i]  = recHits_EEM_RMS[i]/alive_channels_EEM[i] ;
        recHits_EEM_RMS[i]  = sqrt(recHits_EEM_RMS[i] - recHits_EEM_mean[i]*recHits_EEM_mean[i]) ;
    }

    for ( int iX = 1; iX < 101; iX++ )
        for ( int iY = 1; iY < 101; iY++ )
        {
            int iRing = (int) sqrt( (float) (iX-51)*(iX-51) + (float)(iY-51)*(iY-51) ) ;
            if ( h_recHits_EEM_occupancy -> GetBinContent(iX, iY) == 0 || recHits_EEM_RMS[iRing] == 0 ) continue ;

            h_recHits_EEM_deviation -> SetBinContent(iX, iY,
                    ( fabs(h_recHits_EEM_occupancy -> GetBinContent(iX, iY) - recHits_EEM_mean[iRing]))/recHits_EEM_RMS[iRing]
                                                    );

        }

}


// ------------ Pi0 peak  ECAL BARREL ------------
void EcalValidationAOD::doPi0Barrel ( const CaloGeometry *geometry,
                                      const CaloTopology *topology,
                                      edm::Handle<EcalRecHitCollection> recHitsEB )
{
    const CaloSubdetectorGeometry *geometry_EB = geometry->getSubdetectorGeometry(DetId::Ecal,EcalBarrel);
    const CaloSubdetectorGeometry *geometry_ES = geometry->getSubdetectorGeometry(DetId::Ecal,EcalPreshower);

    const CaloSubdetectorTopology *topology_EB = topology->getSubdetectorTopology(DetId::Ecal,EcalBarrel);

    const EcalRecHitCollection* theBarrelEcalRecHits = recHitsEB.product () ;

    //---------------------Barrel--------------------------------//
    if (isMonEBpi0_ )
    {
        if (recHitsEB.isValid() && (recHitsEB->size() > 0))
        {

            vector<EBDetId> maskedIds;
            maskedIds.clear();

            if(isMaskEB_)
            {
                //Read the mask
                std::ifstream mask_file(maskEBFile_.c_str());
                if (mask_file.is_open())
                {

                    char buffer[100];
                    int ieta, iphi, ret;
                    while( ! mask_file.getline(buffer,100).eof() ) {
                        ret=sscanf(buffer,"%d %d ", &ieta, &iphi);
                        if ((ret!=2))
                        {
                            //cout<< " Mask file "<<maskEBFile_.c_str()<<" contains wrong info - must be two digits per line "<<endl;
                            continue;
                        }
                        // cout <<" EB Mask: ieta iphi "<<ieta<<" "<<iphi<<endl;

                        EBDetId maskId(ieta,iphi);
                        maskedIds.push_back(maskId);

                    } // end buffer reading

                    mask_file.close();
                } else {
                    //cout<< " EB Mask File "<<maskEBFile_.c_str()<<" was not opened properly. Status: "<<mask_file.is_open()<<endl;
                }//mask_file.isopen

                //cout<< " EB Pi0: # of masked Xtals: "<<maskedIds.size()<<endl;

            }// isMaskEB



            std::vector<EcalRecHit> seeds;
            seeds.clear();

            vector<EBDetId> usedXtals;
            usedXtals.clear();

            std::vector<EBDetId> detIdEBRecHits;
            detIdEBRecHits.clear();

            std::vector<EcalRecHit> EBRecHits;
            EBRecHits.clear();

            float etot =0;
            EcalRecHitCollection::const_iterator itb;

            // Build Seeds
            for(itb=recHitsEB->begin(); itb!=recHitsEB->end(); ++itb) {

                if ( (useRecoFlag_) && (itb->recoFlag() != 0) ) continue;

                EBDetId id(itb->id());
                double energy = itb->energy();
                if( energy < seleXtalMinEnergy_) continue;

                EBDetId det = itb->id();

                if(isMaskEB_)
                {

                    detIdEBRecHits.push_back(det);
                    EBRecHits.push_back(*itb);

                    bool  maskedId=false;
                    std::vector<EBDetId>::const_iterator mIds;
                    for(mIds=maskedIds.begin(); mIds!=maskedIds.end(); mIds++) {
                        if(*mIds==id) {
                            //cout<< " EB Seed was masked : E ieta iphi recoFlag "<<itb->energy()<<" "<<id.ieta()<<" "<<id.iphi()<<" "<<itb->recoFlag()<<endl;
                            maskedId = true;
                            break;
                        }
                    }

                    if ( (energy > clusSeedThr_) && (!(maskedId)) ) seeds.push_back(*itb);

                    if(!(maskedId))
                    {
                        etot+= itb->energy();
                    }

                }
                else
                {

                    detIdEBRecHits.push_back(det);
                    EBRecHits.push_back(*itb);
                    if ( (energy > clusSeedThr_) ) seeds.push_back(*itb);

                    etot+= itb->energy();

                }

            } // Eb rechits


            int nClus;
            vector<float> eClus;
            vector<float> etClus;
            vector<float> etaClus;
            vector<float> thetaClus;
            vector<float> phiClus;
            vector<EBDetId> max_hit;

            vector< vector<EcalRecHit> > RecHitsCluster;
            vector< vector<EcalRecHit> > RecHitsCluster5x5;
            vector<float> s4s9Clus;
            vector<float> s9s25Clus;


            nClus=0;

            // Make own simple clusters (3x3, 5x5 or clusPhiSize_ x clusEtaSize_)
            sort(seeds.begin(), seeds.end(), ecalRecHitLess());

            //Cycle on seeds
            for (std::vector<EcalRecHit>::iterator itseed=seeds.begin(); itseed!=seeds.end(); itseed++) {
                EBDetId seed_id = itseed->id();
                std::vector<EBDetId>::const_iterator usedIds;

                bool seedAlreadyUsed=false;
                for(usedIds=usedXtals.begin(); usedIds!=usedXtals.end(); usedIds++) {
                    if(*usedIds==seed_id) {
                        seedAlreadyUsed=true;
                        //cout<< " Seed with energy "<<itseed->energy()<<" was used !"<<endl;
                        break;
                    }
                }
                if(seedAlreadyUsed)continue;

                std::vector<DetId> clus_v = topology_EB->getWindow(seed_id,clusEtaSize_,clusPhiSize_);
                std::vector<std::pair<DetId,float> > clus_used;

                vector<EcalRecHit> RecHitsInWindow;
                vector<EcalRecHit> RecHitsInWindow5x5;

                double simple_energy = 0;

                for (std::vector<DetId>::iterator det=clus_v.begin(); det!=clus_v.end(); det++) {
                    EBDetId EBdet = *det;
                    //      cout<<" det "<< EBdet<<" ieta "<<EBdet.ieta()<<" iphi "<<EBdet.iphi()<<endl;
                    bool  HitAlreadyUsed=false;
                    for(usedIds=usedXtals.begin(); usedIds!=usedXtals.end(); usedIds++) {
                        if(*usedIds==*det) {
                            HitAlreadyUsed=true;
                            break;
                        }
                    }
                    if(HitAlreadyUsed)continue;


                    std::vector<EBDetId>::iterator itdet = find( detIdEBRecHits.begin(),detIdEBRecHits.end(),EBdet);
                    if(itdet == detIdEBRecHits.end()) continue;

                    int nn = int(itdet - detIdEBRecHits.begin());
                    usedXtals.push_back(*det);
                    RecHitsInWindow.push_back(EBRecHits[nn]);
                    clus_used.push_back(std::make_pair(*det,1));
                    simple_energy = simple_energy + EBRecHits[nn].energy();

                }

                if(simple_energy <= 0) continue;

                math::XYZPoint clus_pos = posCalculator_.Calculate_Location(clus_used,theBarrelEcalRecHits,geometry_EB,geometry_ES);

                float theta_s = 2. * atan(exp(-clus_pos.eta()));
                float et_s = simple_energy * sin(theta_s);

                //Compute S4/S9 variable
                //We are not sure to have 9 RecHits so need to check eta and phi:
                ///check s4s9
                float s4s9_tmp[4];
                for(int i=0; i<4; i++)s4s9_tmp[i]= 0;

                int seed_ieta = seed_id.ieta();
                int seed_iphi = seed_id.iphi();

                convxtalid( seed_iphi,seed_ieta);

                float e3x3 = 0;
                float e5x5 = 0;

                for(unsigned int j=0; j<RecHitsInWindow.size(); j++) {
                    EBDetId det = (EBDetId)RecHitsInWindow[j].id();

                    int ieta = det.ieta();
                    int iphi = det.iphi();

                    convxtalid(iphi,ieta);

                    float en = RecHitsInWindow[j].energy();

                    int dx = diff_neta_s(seed_ieta,ieta);
                    int dy = diff_nphi_s(seed_iphi,iphi);

                    if(abs(dx)<=1 && abs(dy)<=1) {
                        e3x3 += en;
                        if(dx <= 0 && dy <=0) s4s9_tmp[0] += en;
                        if(dx >= 0 && dy <=0) s4s9_tmp[1] += en;
                        if(dx <= 0 && dy >=0) s4s9_tmp[2] += en;
                        if(dx >= 0 && dy >=0) s4s9_tmp[3] += en;
                    }


                }

                if(e3x3 <= 0) continue;

                float s4s9_max = *max_element( s4s9_tmp,s4s9_tmp+4)/e3x3;

                ///calculate e5x5
                std::vector<DetId> clus_v5x5 = topology_EB->getWindow(seed_id,5,5);
                for( std::vector<DetId>::const_iterator idItr = clus_v5x5.begin(); idItr != clus_v5x5.end(); idItr++) {
                    EBDetId det = *idItr;

                    //inside collections
                    std::vector<EBDetId>::iterator itdet = find( detIdEBRecHits.begin(),detIdEBRecHits.end(),det);
                    if(itdet == detIdEBRecHits.end()) continue;

                    int nn = int(itdet - detIdEBRecHits.begin());

                    RecHitsInWindow5x5.push_back(EBRecHits[nn]);
                    e5x5 += EBRecHits[nn].energy();

                }



                if(e5x5 <= 0) continue;

                eClus.push_back(simple_energy);
                etClus.push_back(et_s);
                etaClus.push_back(clus_pos.eta());
                thetaClus.push_back(theta_s);
                phiClus.push_back(clus_pos.phi());
                s4s9Clus.push_back(s4s9_max);
                s9s25Clus.push_back(e3x3/e5x5);
                RecHitsCluster.push_back(RecHitsInWindow);
                RecHitsCluster5x5.push_back(RecHitsInWindow5x5);

                nClus++;

            }

            // Selection, based on Simple clustering
            //pi0 candidates
            int npi0_s=0;

            for(Int_t i=0 ; i<nClus ; i++) {
                for(Int_t j=i+1 ; j<nClus ; j++) {
                    //cout<<" i "<<i<<"  etClus[i] "<<etClus[i]<<" j "<<j<<"  etClus[j] "<<etClus[j]<<endl;
                    if( etClus[i]>selePtGamma_ && etClus[j]>selePtGamma_ && s4s9Clus[i]>seleS4S9Gamma_ && s4s9Clus[j]>seleS4S9Gamma_) {


                        float p0x = etClus[i] * cos(phiClus[i]);
                        float p1x = etClus[j] * cos(phiClus[j]);
                        float p0y = etClus[i] * sin(phiClus[i]);
                        float p1y = etClus[j] * sin(phiClus[j]);
                        float p0z = eClus[i] * cos(thetaClus[i]);
                        float p1z = eClus[j] * cos(thetaClus[j]);


                        float pt_pair = sqrt( (p0x+p1x)*(p0x+p1x) + (p0y+p1y)*(p0y+p1y));

                        if (pt_pair < selePtPi0_)continue;

                        float m_inv = sqrt ( (eClus[i] + eClus[j])*(eClus[i] + eClus[j]) - (p0x+p1x)*(p0x+p1x) - (p0y+p1y)*(p0y+p1y) - (p0z+p1z)*(p0z+p1z) );
                        if ( (m_inv<seleMinvMaxPi0_) && (m_inv>seleMinvMinPi0_) ) {

                            //New Loop on cluster to measure isolation:
                            vector<int> IsoClus;
                            IsoClus.clear();
                            float Iso = 0;
                            TVector3 pairVect = TVector3((p0x+p1x), (p0y+p1y), (p0z+p1z));
                            for(Int_t k=0 ; k<nClus ; k++) {


                                if(etClus[k] < ptMinForIsolation_) continue;

                                if(k==i || k==j)continue;
                                TVector3 ClusVect = TVector3(etClus[k] *cos(phiClus[k]), etClus[k] * sin(phiClus[k]) , eClus[k] * cos(thetaClus[k]));

                                float dretacl = fabs(etaClus[k] - pairVect.Eta());
                                float drcl = ClusVect.DeltaR(pairVect);
                                //cout<< "   Iso: k, E, drclpi0, detaclpi0, dphiclpi0 "<<k<<" "<<eClus[k]<<" "<<drcl<<" "<<dretacl<<endl;
                                if((drcl<selePi0BeltDR_) && (dretacl<selePi0BeltDeta_) ) {
                                    //cout<< "   ... good iso cluster #: "<<k<<" etClus[k] "<<etClus[k] <<endl;
                                    Iso = Iso + etClus[k];
                                    IsoClus.push_back(k);
                                }
                            }

                            if(Iso/pt_pair<selePi0Iso_) {

                                h_Pi0_EB_eta->Fill(pairVect.Eta()/*, MyWeight*/);
                                h_Pi0_EB_phi->Fill(pairVect.Phi()/*, MyWeight*/);

                                int order[2] = {0,0};

                                if(etClus[i] > etClus[j]) {
                                    order[0] = i;
                                    order[1] = j;
                                }
                                if(etClus[j] > etClus[i]) {
                                    order[0] = j;
                                    order[1] = i;
                                }

                                h_Pi0_EB_mass->           Fill(m_inv/*, MyWeight*/);
                                h_Pi0_EB_pt1 ->Fill(etClus[order[0]]/*, MyWeight*/);
                                h_Pi0_EB_pt2 ->Fill(etClus[order[1]]/*, MyWeight*/);
                                h_Pi0_EB_pt  ->         Fill(pt_pair/*, MyWeight*/);

                                npi0_s++;

                            }
                        }
                    }
                } // End of the "j" loop over Simple Clusters
            } // End of the "i" loop over Simple Clusters
        }

    }// End isMonEBpi0_

    //------------------ End of pi0 in EB --------------------------//
}


// ------------ Pi0 peak  ECAL BARREL ------------

void EcalValidationAOD::doPi0Endcap ( const CaloGeometry *geometry,
                                      const CaloTopology *topology,
                                      edm::Handle<EcalRecHitCollection> recHitsEE )
{
    const CaloSubdetectorGeometry *geometry_EE = geometry->getSubdetectorGeometry(DetId::Ecal,EcalEndcap);
    const CaloSubdetectorGeometry *geometry_ES = geometry->getSubdetectorGeometry(DetId::Ecal,EcalPreshower);

    const CaloSubdetectorTopology *topology_EE = topology->getSubdetectorTopology(DetId::Ecal,EcalEndcap);

    const EcalRecHitCollection* theEndcapEcalRecHits = recHitsEE.product () ;

    if (isMonEEpi0_ )
    {
        if (recHitsEE.isValid() && (recHitsEE->size() > 0))
        {
            vector<EEDetId> maskedIds;
            maskedIds.clear();

            if(isMaskEE_)
            {
                //Read the mask
                std::ifstream mask_file(maskEEFile_.c_str());
                if (mask_file.is_open())
                {
                    char buffer[100];
                    int ix, iy, iz, ret;
                    while( ! mask_file.getline(buffer,100).eof() ) {
                        ret=sscanf(buffer,"%d %d %d", &iz, &ix, &iy);
                        if ((ret!=3))
                        {
                            //cout<< " Mask file "<<maskEEFile_.c_str()<<" contains wrong info - must be three digits per line "<<endl;
                            continue;
                        }
                        //cout <<" EE Mask: iz ix iy "<<iz<<" "<<ix<<" "<<iy<<endl;

                        EEDetId maskId(ix,iy,iz);
                        maskedIds.push_back(maskId);

                    } // end buffer reading

                    mask_file.close();
                } else {
                    //cout<< " EE Mask File "<<maskEEFile_.c_str()<<" was not opened properly. Status: "<<mask_file.is_open()<<endl;
                }//mask_file.isopen

                //cout<< " EE Pi0: # of masked Xtals: "<<maskedIds.size()<<endl;

            }// isMaskEE


            std::vector<EcalRecHit> seeds_EE;
            seeds_EE.clear();

            vector<EEDetId> usedXtals_EE;
            usedXtals_EE.clear();

            std::vector<EEDetId> detIdEERecHits;
            detIdEERecHits.clear();

            std::vector<EcalRecHit> EERecHits;
            EERecHits.clear();

            float etot_EE =0;
            EERecHitCollection::const_iterator ite;

            // Build Seeds
            for(ite=recHitsEE->begin(); ite!=recHitsEE->end(); ++ite) {

                if ( (useRecoFlag_) && (ite->recoFlag() != 0) ) continue;

                double energy = ite->energy();
                if( energy < seleXtalMinEnergy_EE_) continue;

                EEDetId det = ite->id();
                EEDetId id(ite->id());

                if(isMaskEE_)
                {

                    detIdEERecHits.push_back(det);
                    EERecHits.push_back(*ite);

                    bool  maskedId=false;
                    std::vector<EEDetId>::const_iterator mIds;
                    for(mIds=maskedIds.begin(); mIds!=maskedIds.end(); mIds++) {
                        if(*mIds==id) {
                            //cout<< " EE Seed was masked : E ix iy iz "<<ite->energy()<<" "<<id.ix()<<" "<<id.iy()<<" "<<id.zside()<<endl;
                            maskedId = true;
                            break;
                        }
                    }

                    if ( (energy > clusSeedThr_EE_) && (!(maskedId)) ) seeds_EE.push_back(*ite);

                    if(!(maskedId))
                    {
                        etot_EE+= ite->energy();
                    }


                }
                else
                {
                    detIdEERecHits.push_back(det);
                    EERecHits.push_back(*ite);

                    if ( (energy > clusSeedThr_EE_) ) seeds_EE.push_back(*ite);

                    etot_EE+= ite->energy();

                } //ismaskEE
            } // Ee rechits


            int nClus_EE;
            vector<float> eClus_EE;
            vector<float> etClus_EE;
            vector<float> etaClus_EE;
            vector<float> thetaClus_EE;
            vector<float> phiClus_EE;
            vector< vector<EcalRecHit> > RecHitsCluster_EE;
            vector< vector<EcalRecHit> > RecHitsCluster_EE_5x5;
            vector<float> s4s9Clus_EE;
            vector<float> s9s25Clus_EE;
            vector<float> xClusEndCap;
            vector<float> yClusEndCap;
            vector<float> zClusEndCap;



            nClus_EE=0;

            // Make own simple clusters (3x3, 5x5 or clusPhiSize_ x clusEtaSize_)
            sort(seeds_EE.begin(), seeds_EE.end(), ecalRecHitLess());

            //Cycle on seeds_EE
            for (std::vector<EcalRecHit>::iterator itseed=seeds_EE.begin(); itseed!=seeds_EE.end(); itseed++) {
                EEDetId seed_id = itseed->id();
                std::vector<EEDetId>::const_iterator usedIds;

                bool seedAlreadyUsed=false;
                for(usedIds=usedXtals_EE.begin(); usedIds!=usedXtals_EE.end(); usedIds++) {
                    if(*usedIds==seed_id) {
                        seedAlreadyUsed=true;
                        //cout<< " Seed with energy "<<itseed->energy()<<" was used !"<<endl;
                        break;
                    }
                }
                if(seedAlreadyUsed)continue;

                std::vector<DetId> clus_v = topology_EE->getWindow(seed_id,clusEtaSize_,clusPhiSize_);
                std::vector<std::pair<DetId,float> > clus_used;

                vector<EcalRecHit> RecHitsInWindow;
                vector<EcalRecHit> RecHitsInWindow5x5;

                double simple_energy = 0;

                for (std::vector<DetId>::iterator det=clus_v.begin(); det!=clus_v.end(); det++) {
                    EEDetId EEdet = *det;
                    //      cout<<" det "<< EBdet<<" ieta "<<EBdet.ieta()<<" iphi "<<EBdet.iphi()<<endl;
                    bool  HitAlreadyUsed=false;
                    for(usedIds=usedXtals_EE.begin(); usedIds!=usedXtals_EE.end(); usedIds++) {
                        if(*usedIds==*det) {
                            HitAlreadyUsed=true;
                            break;
                        }
                    }
                    if(HitAlreadyUsed)continue;


                    std::vector<EEDetId>::iterator itdet = find( detIdEERecHits.begin(),detIdEERecHits.end(),EEdet);
                    if(itdet == detIdEERecHits.end()) continue;

                    int nn = int(itdet - detIdEERecHits.begin());
                    usedXtals_EE.push_back(*det);
                    RecHitsInWindow.push_back(EERecHits[nn]);
                    clus_used.push_back(std::make_pair(*det,1));
                    simple_energy = simple_energy + EERecHits[nn].energy();

                }

                if(simple_energy <= 0) continue;

                math::XYZPoint clus_pos = posCalculator_.Calculate_Location(clus_used, theEndcapEcalRecHits, geometry_EE, geometry_ES);

                float theta_s = 2. * atan(exp(-clus_pos.eta()));
                float et_s = simple_energy * sin(theta_s);

                //Compute S4/S9 variable
                //We are not sure to have 9 RecHits so need to check eta and phi:
                ///check s4s9
                float s4s9_tmp[4];
                for(int i=0; i<4; i++)s4s9_tmp[i]= 0;

                int ixSeed = seed_id.ix();
                int iySeed = seed_id.iy();
                float e3x3 = 0;
                float e5x5 = 0;

                for(unsigned int j=0; j<RecHitsInWindow.size(); j++) {
                    EEDetId det_this = (EEDetId)RecHitsInWindow[j].id();

                    int dx = ixSeed - det_this.ix();
                    int dy = iySeed - det_this.iy();

                    float en = RecHitsInWindow[j].energy();

                    if(abs(dx)<=1 && abs(dy)<=1) {
                        e3x3 += en;
                        if(dx <= 0 && dy <=0) s4s9_tmp[0] += en;
                        if(dx >= 0 && dy <=0) s4s9_tmp[1] += en;
                        if(dx <= 0 && dy >=0) s4s9_tmp[2] += en;
                        if(dx >= 0 && dy >=0) s4s9_tmp[3] += en;
                    }


                }

                if(e3x3 <= 0) continue;

                float s4s9_max = *max_element( s4s9_tmp,s4s9_tmp+4)/e3x3;

                ///calculate e5x5
                std::vector<DetId> clus_v5x5 = topology_EE->getWindow(seed_id,5,5);
                for( std::vector<DetId>::const_iterator idItr = clus_v5x5.begin(); idItr != clus_v5x5.end(); idItr++) {
                    EEDetId det = *idItr;

                    //inside collections
                    std::vector<EEDetId>::iterator itdet = find( detIdEERecHits.begin(),detIdEERecHits.end(),det);
                    if(itdet == detIdEERecHits.end()) continue;

                    int nn = int(itdet - detIdEERecHits.begin());

                    RecHitsInWindow5x5.push_back(EERecHits[nn]);
                    e5x5 += EERecHits[nn].energy();

                }



                if(e5x5 <= 0) continue;



                xClusEndCap.push_back(clus_pos.x());
                yClusEndCap.push_back(clus_pos.y());
                zClusEndCap.push_back(clus_pos.z());

                eClus_EE.push_back(simple_energy);
                etClus_EE.push_back(et_s);
                etaClus_EE.push_back(clus_pos.eta());
                thetaClus_EE.push_back(theta_s);
                phiClus_EE.push_back(clus_pos.phi());
                s4s9Clus_EE.push_back(s4s9_max);
                s9s25Clus_EE.push_back(e3x3/e5x5);
                RecHitsCluster_EE.push_back(RecHitsInWindow);
                RecHitsCluster_EE_5x5.push_back(RecHitsInWindow5x5);

                nClus_EE++;

            }

            // Selection, based on Simple clustering
            //pi0 candidates
            int npi0_s=0;

            for(Int_t i=0 ; i<nClus_EE ; i++) {
                for(Int_t j=i+1 ; j<nClus_EE ; j++) {
                    //cout<<" i "<<i<<"  etClus_EE[i] "<<etClus_EE[i]<<" j "<<j<<"  etClus_EE[j] "<<etClus_EE[j]<<endl;
                    if( s4s9Clus_EE[i] < seleS4S9Gamma_EE_ || s4s9Clus_EE[j] < seleS4S9Gamma_EE_) continue ;

                    float p0x = etClus_EE[i] * cos(phiClus_EE[i]);
                    float p1x = etClus_EE[j] * cos(phiClus_EE[j]);
                    float p0y = etClus_EE[i] * sin(phiClus_EE[i]);
                    float p1y = etClus_EE[j] * sin(phiClus_EE[j]);
                    float p0z = eClus_EE[i] * cos(thetaClus_EE[i]);
                    float p1z = eClus_EE[j] * cos(thetaClus_EE[j]);


                    float pt_pair = sqrt( (p0x+p1x)*(p0x+p1x) + (p0y+p1y)*(p0y+p1y));

                    if (pt_pair < selePtPi0_)continue;

                    float m_inv = sqrt ( (eClus_EE[i] + eClus_EE[j])*(eClus_EE[i] + eClus_EE[j]) - (p0x+p1x)*(p0x+p1x) - (p0y+p1y)*(p0y+p1y) - (p0z+p1z)*(p0z+p1z) );

                    if ( m_inv > seleMinvMaxPi0_EE_ || m_inv < seleMinvMinPi0_EE_) continue;

                    ////try different cut for different regions
                    TVector3 pairVect = TVector3((p0x+p1x), (p0y+p1y), (p0z+p1z));
                    float etapair = fabs(pairVect.Eta());
                    float ptmin = etClus_EE[i] < etClus_EE[j] ?  etClus_EE[i] : etClus_EE[j];

                    if(etapair <= region1_Pi0_EE_) {
                        if(ptmin < selePtGammaPi0_EE_region1_ || pt_pair < selePtPi0_EE_region1_) continue;
                    } else if( etapair <= region2_Pi0_EE_) {
                        if(ptmin < selePtGammaPi0_EE_region2_ || pt_pair < selePtPi0_EE_region2_) continue;
                    } else {
                        if(ptmin < selePtGammaPi0_EE_region3_ || pt_pair < selePtPi0_EE_region3_) continue;
                    }


                    //New Loop on cluster to measure isolation:
                    vector<int> IsoClus_EE;
                    IsoClus_EE.clear();
                    float Iso = 0;
                    for(Int_t k=0 ; k<nClus_EE ; k++) {


                        if(etClus_EE[k] < ptMinForIsolation_EE_) continue;

                        if(k==i || k==j)continue;
                        TVector3 ClusVect = TVector3(etClus_EE[k] *cos(phiClus_EE[k]), etClus_EE[k] * sin(phiClus_EE[k]) , eClus_EE[k] * cos(thetaClus_EE[k]));

                        float dretacl = fabs(etaClus_EE[k] - pairVect.Eta());
                        float drcl = ClusVect.DeltaR(pairVect);
                        //cout<< "   Iso: k, E, drclpi0, detaclpi0, dphiclpi0 "<<k<<" "<<eClus_EE[k]<<" "<<drcl<<" "<<dretacl<<endl;
                        if((drcl<selePi0BeltDR_EE_) && (dretacl<selePi0BeltDeta_EE_) ) {
                            //cout<< "   ... good iso cluster #: "<<k<<" etClus_EE[k] "<<etClus_EE[k] <<endl;
                            Iso = Iso + etClus_EE[k];
                            IsoClus_EE.push_back(k);
                        }
                    }

                    if(Iso/pt_pair<selePi0Iso_EE_) {

                        h_Pi0_EE_eta->Fill(pairVect.Eta()/*, MyWeight*/);
                        h_Pi0_EE_phi->Fill(pairVect.Phi()/*, MyWeight*/);

                        int order[2] = {0,0};

                        if(etClus_EE[i] > etClus_EE[j]) {
                            order[0] = i;
                            order[1] = j;
                        }
                        if(etClus_EE[j] > etClus_EE[i]) {
                            order[0] = j;
                            order[1] = i;
                        }

                        h_Pi0_EE_mass->              Fill(m_inv/*, MyWeight*/);
                        h_Pi0_EE_pt1 ->Fill(etClus_EE[order[0]]/*, MyWeight*/);
                        h_Pi0_EE_pt2 ->Fill(etClus_EE[order[1]]/*, MyWeight*/);
                        h_Pi0_EE_pt  ->            Fill(pt_pair/*, MyWeight*/);

                        npi0_s++;

                    }
                } // End of the "j" loop over Simple Clusters
            } // End of the "i" loop over Simple Clusters

        }//------------------ End of pi0 in EE --------------------------//

    }//End isMonEEpi0_

}



// ----------additional functions-------------------

void EcalValidationAOD::convxtalid(Int_t &nphi,Int_t &neta)
{
    // Barrel only
    // Output nphi 0...359; neta 0...84; nside=+1 (for eta>0), or 0 (for eta<0).
    // neta will be [-85,-1] , or [0,84], the minus sign indicates the z<0 side.

    if(neta > 0) neta -= 1;
    if(nphi > 359) nphi=nphi-360;

} //end of convxtalid

int EcalValidationAOD::diff_neta_s(Int_t neta1, Int_t neta2) {
    Int_t mdiff;
    mdiff=(neta1-neta2);
    return mdiff;
}

// Calculate the distance in xtals taking into account the periodicity of the Barrel
int EcalValidationAOD::diff_nphi_s(Int_t nphi1,Int_t nphi2) {
    Int_t mdiff;
    if(abs(nphi1-nphi2) < (360-abs(nphi1-nphi2))) {
        mdiff=nphi1-nphi2;
    }
    else {
        mdiff=360-abs(nphi1-nphi2);
        if(nphi1>nphi2) mdiff=-mdiff;
    }
    return mdiff;
}

//define this as a plug-in
DEFINE_FWK_MODULE(EcalValidationAOD);
