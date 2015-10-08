#ifndef CONST_H
#define CONST_H

static const TString sBRf = "5e2";
static float BRfactor     = 5e2;
static float BRbelle      = 2.1e-8;

static const float luminosity = 20.28;
static const TString slumi = "20.3 fb^{-1}";

static const float fullluminosity = 20.28;

static const float MeV2GeV = 1.e-3;
static const float GeV2MeV = 1.e+3;

static const float mTrainingMinGlob = 750; // 850
static const float mTrainingMaxGlob = 2500; // 2710
static const float mBlindMinGlob = 1690; // since the bin size is 30 MeV in the m3body histo
static const float mBlindMaxGlob = 1870; // since the bin size is 30 MeV in the m3body histo
static const float mSideBandLeftLowerMeVGlob  = 1450;
static const float mSideBandRightUpperMeVGlob = 2110;
static const float mSRminMeVGlob = 1713;
static const float mSRmaxMeVGlob = 1841;
static const float mBinSize = 30; // MeV

static const float optBDTcut = 0.933;
static const float minBDTcut = -0.9;
static const int   nBDTbins  = 38;

static const float SigmaRhoOmegaPhiMeV = 30;
static const float SigmaRhoOmegaPhiOS1MeV = 30;
static const float SigmaRhoOmegaPhiOS2MeV = 27;

static const double mb2fb = 1.e12;
static const double nb2fb = 1.e+6;
static const double pb2fb = 1.e+3;
static const double nb2mb = 1.e-6;

static const float muonMass    = 0.105658367; // GeV
static const float muonMassMeV = 105.658367; // MeV
static const float tauMassGeV  = 1.77682; // GeV
static const float tauMassMeV  = 1776.82; // MeV
static const float elecMassGeV = 0.00051099891; // GeV
static const float elecMassMeV = 0.51099891; // MeV
static const float kaonMassMeV = 493.677; // MeV
static const float pionMassMeV = 139.57018; // MeV

#endif
