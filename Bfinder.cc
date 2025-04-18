// vim:set ts=4 sw=4 fdm=marker et:
// Ntuplt creator for B meson related analysis.
// Maintain and contact: ta-wei wang
// Email: "tawei@mit.edu" or "ta-wei.wang@cern.ch"
#include "Bfinder/Bfinder/interface/format.h"
#include "Bfinder/Bfinder/interface/Bntuple.h"
#include "Bfinder/Bfinder/interface/utilities.h"
//
// class declaration
//

class Bfinder : public edm::one::EDAnalyzer<edm::one::WatchRuns>
{//{{{
public:
  explicit Bfinder(const edm::ParameterSet&);
  ~Bfinder();
  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
 
private:
  virtual void beginJob() ;
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;
 
  virtual void beginRun(edm::Run const&, edm::EventSetup const&);
  virtual void endRun(edm::Run const&, edm::EventSetup const&);
  virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);
  virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);
        
  virtual void BranchOut2MuTk(
                              BInfoBranches &BInfo,
                              // std::vector<pat::PackedCandidate> input_tracks,
                              std::vector<const reco::Track*> input_tracks,
                              reco::Vertex thePrimaryV,
                              std::vector<bool> isNeededTrack,
                              TLorentzVector v4_mu1,
                              TLorentzVector v4_mu2,
                              reco::TransientTrack muonPTT,
                              reco::TransientTrack muonMTT,
                              std::vector<int> &B_counter,
                              float *mass_window,
                              float MuMu_MASS,
                              float Tk_MASS,
                              int channel_number
                              );
  virtual void BranchOut2MuX_XtoTkTk(
                                     BInfoBranches &BInfo,
                                     // std::vector<pat::PackedCandidate> input_tracks,
                                     std::vector<const reco::Track*> input_tracks,
                                     reco::Vertex thePrimaryV,
                                     std::vector<bool> isNeededTrack,
                                     TLorentzVector v4_mu1,
                                     TLorentzVector v4_mu2,
                                     reco::TransientTrack muonPTT,
                                     reco::TransientTrack muonMTT,
                                     std::vector<int> &B_counter,
                                     float *mass_window,
                                     float MuMu_MASS,
                                     float TkTk_MASS,
                                     float TkTk_window,
                                     float Tk1_MASS,
                                     float Tk2_MASS,     
                                     int channel_number,
                                     int fit_option
                                     );

  static const reco::Track* getFromPC(const reco::Candidate& rc) {
    const pat::PackedCandidate* pp = dynamic_cast<const pat::PackedCandidate*>(&rc);
    if (pp == nullptr)
      return nullptr;
    try {
      const reco::Track* trk = &pp->pseudoTrack();
      return trk;
    } catch (edm::Exception const&) {
    }
    return nullptr;
  }
  
  // ----------member data ---------------------------
  edm::ESHandle<MagneticField> bField;
  edm::ESGetToken<MagneticField, IdealMagneticFieldRecord> idealMagneticFieldRecordToken_;
  edm::ParameterSet theConfig;

  bool detailMode_;
  bool dropUnusedTracks_;
  //std::vector<std::string> TriggersForMatching_;
  std::vector<int> Bchannel_;
  std::vector<std::string> MuonTriggerMatchingPath_;
  std::vector<std::string> MuonTriggerMatchingFilter_;
  //edm::InputTag hltLabel_;
  edm::EDGetTokenT< reco::GenParticleCollection > genLabel_;
  edm::EDGetTokenT< edm::View<pat::Muon> > muonLabel_;
  edm::EDGetTokenT< edm::View<pat::PackedCandidate> > trackLabel_;
  edm::EDGetTokenT< edm::ValueMap<float> > chi2Map_;
  edm::EDGetTokenT< edm::ValueMap<reco::DeDxData> > dedxMap_;
  edm::EDGetTokenT< reco::BeamSpot > bsLabel_;
  edm::EDGetTokenT< reco::VertexCollection > pvLabel_;

  double tkPtCut_;
  double tkEtaCut_;
  double jpsiPtCut_;
  double uj_VtxChiProbCut_;
  std::vector<double> bPtCut_;
  std::vector<double> bEtaCut_;
  std::vector<double> VtxChiProbCut_;
  std::vector<double> svpvDistanceCut_;
  std::vector<double> MaxDocaCut_;
  std::vector<double> alphaCut_;
  bool doTkPreCut_;
  bool doMuPreCut_;
  bool makeBntuple_;
  bool doBntupleSkim_;
  bool printInfo_;
  bool readDedx_;
  // edm::EDGetTokenT<edm::ValueMap<float> > MVAMapLabel_;
  // edm::EDGetTokenT< std::vector<float> > MVAMapLabelpA_;
  // edm::InputTag MVAMapLabelInputTag_;

  edm::Service<TFileService> fs;
  TTree *root;
  EvtInfoBranches     EvtInfo;
  VtxInfoBranches     VtxInfo;
  MuonInfoBranches    MuonInfo;
  TrackInfoBranches   TrackInfo;
  BInfoBranches       BInfo;
  GenInfoBranches     GenInfo;
  CommonFuncts        Functs;
  BntupleBranches     *Bntuple = new BntupleBranches;
  TTree* nt0;
  TTree* nt1;
  TTree* nt2;
  TTree* nt3;
  TTree* nt5;
  TTree* nt6;
  TTree* nt7;
  TTree* ntGen;

  //histograms
  TH1F *MuonCutLevel;
  TH1F *TrackCutLevel;
  TH1F *XbujCutLevel;
  //How many channel
  static int const Nchannel = 20;
  std::vector<TH1F*> XbMassCutLevel;

};//}}}

void Bfinder::beginJob()
{//{{{
  root = fs->make<TTree>("root","root");
  nt0   = fs->make<TTree>("ntKp","");     Bntuple->buildBranch(nt0);
  nt1   = fs->make<TTree>("ntpi","");     Bntuple->buildBranch(nt1);
  nt2   = fs->make<TTree>("ntKs","");     Bntuple->buildBranch(nt2);
  nt3   = fs->make<TTree>("ntKstar","");  Bntuple->buildBranch(nt3);
  nt5   = fs->make<TTree>("ntphi","");    Bntuple->buildBranch(nt5);
  nt6   = fs->make<TTree>("ntmix","");    Bntuple->buildBranch(nt6);
  nt7   = fs->make<TTree>("ntJpsi","");   Bntuple->buildBranch(nt7,true);
  ntGen = fs->make<TTree>("ntGen","");    Bntuple->buildGenBranch(ntGen);
  EvtInfo.regTree(root);
  VtxInfo.regTree(root);
  MuonInfo.regTree(root, detailMode_, MuonTriggerMatchingPath_.size(), MuonTriggerMatchingFilter_.size());
  TrackInfo.regTree(root, detailMode_);
  BInfo.regTree(root, detailMode_);
  GenInfo.regTree(root);
}//}}}

Bfinder::Bfinder(const edm::ParameterSet& iConfig):theConfig(iConfig)
{//{{{
  //now do what ever initialization is needed
  detailMode_ = iConfig.getParameter<bool>("detailMode");
  dropUnusedTracks_ = iConfig.getParameter<bool>("dropUnusedTracks");

  idealMagneticFieldRecordToken_ = esConsumes();

  //TriggersForMatching_= iConfig.getUntrackedParameter<std::vector<std::string> >("TriggersForMatching");
  //hltLabel_           = iConfig.getParameter<edm::InputTag>("HLTLabel");
  Bchannel_= iConfig.getParameter<std::vector<int> >("Bchannel");
  MuonTriggerMatchingPath_ = iConfig.getParameter<std::vector<std::string> >("MuonTriggerMatchingPath");
  MuonTriggerMatchingFilter_ = iConfig.getParameter<std::vector<std::string> >("MuonTriggerMatchingFilter");
  genLabel_           = consumes< reco::GenParticleCollection >(iConfig.getParameter<edm::InputTag>("GenLabel"));
  trackLabel_         = consumes< edm::View<pat::PackedCandidate> >(iConfig.getParameter<edm::InputTag>("TrackLabel"));
  chi2Map_ = consumes< edm::ValueMap<float> >(iConfig.getParameter< edm::InputTag >("TrackChi2Label"));
  dedxMap_ = consumes< edm::ValueMap<reco::DeDxData> >(iConfig.getParameter<edm::InputTag>("TrackDedxLabel"));
  muonLabel_          = consumes< edm::View<pat::Muon> >(iConfig.getParameter<edm::InputTag>("MuonLabel"));
  bsLabel_        = consumes< reco::BeamSpot >(iConfig.getParameter<edm::InputTag>("BSLabel"));
  pvLabel_        = consumes< reco::VertexCollection >(iConfig.getParameter<edm::InputTag>("PVLabel"));

  tkPtCut_ = iConfig.getParameter<double>("tkPtCut");
  tkEtaCut_ = iConfig.getParameter<double>("tkEtaCut");
  jpsiPtCut_ = iConfig.getParameter<double>("jpsiPtCut");
  uj_VtxChiProbCut_ = iConfig.getParameter<double>("uj_VtxChiProbCut");
  bPtCut_ = iConfig.getParameter<std::vector<double> >("bPtCut");
  bEtaCut_ = iConfig.getParameter<std::vector<double> >("bEtaCut");
  VtxChiProbCut_ = iConfig.getParameter<std::vector<double> >("VtxChiProbCut");
  svpvDistanceCut_ = iConfig.getParameter<std::vector<double> >("svpvDistanceCut");
  MaxDocaCut_ = iConfig.getParameter<std::vector<double> >("MaxDocaCut");
  alphaCut_ = iConfig.getParameter<std::vector<double> >("alphaCut");
  doTkPreCut_ = iConfig.getParameter<bool>("doTkPreCut");
  doMuPreCut_ = iConfig.getParameter<bool>("doMuPreCut");
  makeBntuple_ = iConfig.getParameter<bool>("makeBntuple");
  doBntupleSkim_ = iConfig.getParameter<bool>("doBntupleSkim");
  printInfo_ = iConfig.getParameter<bool>("printInfo");
  readDedx_ = iConfig.getParameter<bool>("readDedx");
  // MVAMapLabelInputTag_ = iConfig.getParameter<edm::InputTag>("MVAMapLabel");
  // MVAMapLabel_ = consumes<edm::ValueMap<float> >(iConfig.getParameter<edm::InputTag>("MVAMapLabel"));
  // MVAMapLabelpA_ = consumes< std::vector<float> >(iConfig.getParameter<edm::InputTag>("MVAMapLabel"));
  // Dedx_Token1_ = consumes<edm::ValueMap<reco::DeDxData> >(iConfig.getParameter<edm::InputTag>("Dedx_Token1"));
  // Dedx_Token2_ = consumes<edm::ValueMap<reco::DeDxData> >(iConfig.getParameter<edm::InputTag>("Dedx_Token2"));

  MuonCutLevel        = fs->make<TH1F>("MuonCutLevel"     , "MuonCutLevel"    , 10, 0, 10);
  TrackCutLevel       = fs->make<TH1F>("TrackCutLevel"    , "TrackCutLevel"   , 10, 0, 10);
  XbujCutLevel        = fs->make<TH1F>("XbujCutLevel"     , "XbujCutLevel"    , 10, 0, 10);
  for(unsigned int i = 0; i < Bchannel_.size(); i++){
    TH1F* XbMassCutLevel_temp      = fs->make<TH1F>(TString::Format("XbMassCutLevel_i")   ,TString::Format("XbMassCutLevel_i")  , 10, 0, 10);
    XbMassCutLevel.push_back(XbMassCutLevel_temp);
  }
}//}}}

Bfinder::~Bfinder()
{//{{{
  // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)
}//}}}

//
// member functions
//

// ------------ method called for each event  ------------
void Bfinder::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  //checking input parameter size
  if( (Bchannel_.size() != bPtCut_.size()) || (bPtCut_.size() != bEtaCut_.size()) || (bEtaCut_.size() != VtxChiProbCut_.size()) || (VtxChiProbCut_.size() != svpvDistanceCut_.size()) || (svpvDistanceCut_.size() != MaxDocaCut_.size()) || (MaxDocaCut_.size() != alphaCut_.size())){
    std::cout<<"Unmatched input parameter vector size, EXIT"<<std::endl;
    return;
  }
  //checking Trigger matching size
  if(MuonTriggerMatchingPath_.size() > MAX_TRIGGER || MuonTriggerMatchingFilter_.size() > MAX_TRIGGER){
    std::cout<<"Number of trigger matching pathes exceeded limit, EXIT"<<std::endl;
    return;
  }

  //std::cout << "*************************\nReconstructing event number: " << iEvent.id() << "\n";
  using namespace edm;
  using namespace reco;
  //ESHandle<MagneticField> bField;
  bField = iSetup.getHandle(idealMagneticFieldRecordToken_);

  // Change used muon and track collections
  auto muons = iEvent.getHandle(muonLabel_); // edm::Handle< std::vector<pat::Muon> >
  auto tks = iEvent.getHandle( trackLabel_ ); // edm::Handle< std::vector<pat::PackedCandidate> >
  // auto losttks = iEvent.getHandle(losttrackLabel_);
  auto chi2Handle = iEvent.getHandle(chi2Map_); // edm::Handle<edm::ValueMap<float>>
  auto dedxHandle = iEvent.getHandle(dedxMap_); // edm::Handle<edm::ValueMap<reco::DeDxData>>

  edm::Handle<std::vector<reco::GenParticle>> gens;
  if (!iEvent.isRealData()) {
    gens = iEvent.getHandle(genLabel_);
  }
  
  // edm::Handle< std::vector<reco::Track> > etracks;
  // iEvent.getByToken(trackLabelReco_, etracks);
  // if(etracks->size() != tks->size()) 
  //   { 
  //     fprintf(stderr,"ERROR: number of tracks in pat::GenericParticle is different from reco::Track.\n"); 
  //     exit(0);
  //   }

  //CLEAN all memory
  memset(&EvtInfo     ,0x00,sizeof(EvtInfo)   );
  memset(&VtxInfo     ,0x00,sizeof(VtxInfo)   );
  memset(&MuonInfo    ,0x00,sizeof(MuonInfo)  );
  memset(&TrackInfo   ,0x00,sizeof(TrackInfo) );
  memset(&BInfo       ,0x00,sizeof(BInfo)    );
  memset(&GenInfo     ,0x00,sizeof(GenInfo)   );

  // EvtInfo section{{{
  EvtInfo.RunNo   = iEvent.id().run();
  EvtInfo.EvtNo   = iEvent.id().event();
  EvtInfo.BxNo    = iEvent.bunchCrossing();
  EvtInfo.LumiNo  = iEvent.luminosityBlock();
  EvtInfo.Orbit   = iEvent.orbitNumber();
  EvtInfo.McFlag  = !iEvent.isRealData();
  //EvtInfo.hltnames->clear();
  //EvtInfo.nTrgBook= N_TRIGGER_BOOKINGS;

  // Handle primary vertex properties
  Vertex thePrimaryV; //, thePrimaryVmaxPt, thePrimaryVmaxMult;
  math::XYZPoint RefVtx; //, RefVtxmaxPt, RefVtxmaxMult;
  //get beamspot information
  Vertex theBeamSpotV;
  reco::BeamSpot beamSpot;
  edm::Handle<reco::BeamSpot> beamSpotHandle;
  iEvent.getByToken(bsLabel_, beamSpotHandle);
  if (beamSpotHandle.isValid()){
    beamSpot = *beamSpotHandle;
    theBeamSpotV = Vertex(beamSpot.position(), beamSpot.covariance3D());
    EvtInfo.BSx             = beamSpot.x0();
    EvtInfo.BSy             = beamSpot.y0();
    EvtInfo.BSz             = beamSpot.z0();
    EvtInfo.BSxErr          = beamSpot.x0Error();
    EvtInfo.BSyErr          = beamSpot.y0Error();
    EvtInfo.BSzErr          = beamSpot.z0Error();
    EvtInfo.BSdxdz          = beamSpot.dxdz();
    EvtInfo.BSdydz          = beamSpot.dydz();
    EvtInfo.BSdxdzErr       = beamSpot.dxdzError();
    EvtInfo.BSdydzErr       = beamSpot.dydzError();
    EvtInfo.BSWidthX        = beamSpot.BeamWidthX();
    EvtInfo.BSWidthXErr     = beamSpot.BeamWidthXError();
    EvtInfo.BSWidthY        = beamSpot.BeamWidthY();
    EvtInfo.BSWidthYErr     = beamSpot.BeamWidthYError();
  }else{
    std::cout<< "No beam spot available from EventSetup \n";
  }

  //get vertex informationa
  edm::Handle<reco::VertexCollection> VertexHandle;
  iEvent.getByToken(pvLabel_, VertexHandle);

  /*  
      if (!VertexHandle.failedToGet() && VertexHandle->size()>0){
      //int nVtxTrks = 0;//outdated PV definition
      double max_tkSt = 0;
      for(std::vector<reco::Vertex>::const_iterator it_vtx = VertexHandle->begin(); it_vtx != VertexHandle->end(); it_vtx++){
      if (!it_vtx->isValid()) continue;
      //find primary vertex with largest St
      double tkSt = 0;
      for(std::vector<reco::TrackBaseRef>::const_iterator it_tk = it_vtx->tracks_begin();
      it_tk != it_vtx->tracks_end(); it_tk++){
      tkSt += it_tk->get()->pt();
      }
      if (tkSt > max_tkSt){
      max_tkSt = tkSt;
      thePrimaryV = Vertex(*it_vtx);
      }
      }
      }else{ 
      thePrimaryV = Vertex(beamSpot.position(), beamSpot.covariance3D());
      }
      RefVtx = thePrimaryV.position();
  */

  double PVBS_Pt_Max = -100.;
  reco::Vertex PVtx_BS;
  if( VertexHandle.isValid() && !VertexHandle.failedToGet() && VertexHandle->size() > 0) {
    //const vector<reco::Vertex> VerticesBS = *VertexHandle;
    for(std::vector<reco::Vertex>::const_iterator it_vtx = VertexHandle->begin();it_vtx != VertexHandle->end(); it_vtx++ ) {
      if (VtxInfo.Size>=MAX_Vertices) {
        std::cout << "PVBS " << VtxInfo.Size << std::endl;
        fprintf(stderr,"ERROR: number of  Vertices exceeds the size of array.\n");
        break;//exit(0);
      }
      VtxInfo.isValid[VtxInfo.Size] = it_vtx->isValid();
      VtxInfo.isFake[VtxInfo.Size] = it_vtx->isFake();
      VtxInfo.Ndof[VtxInfo.Size] = it_vtx->ndof();
      VtxInfo.NormalizedChi2[VtxInfo.Size] = it_vtx->normalizedChi2();
      VtxInfo.x[VtxInfo.Size] = it_vtx->x(); 
      VtxInfo.y[VtxInfo.Size] = it_vtx->y();
      VtxInfo.z[VtxInfo.Size] = it_vtx->z();
      VtxInfo.Pt_Sum[VtxInfo.Size] = 0.;
      VtxInfo.Pt_Sum2[VtxInfo.Size] = 0.;
      //if its hiSelectedVertex, then there will be only one vertex and will have no associated tracks
      if(int(VertexHandle->end()-VertexHandle->begin())==1){
        thePrimaryV = *it_vtx;
        VtxInfo.Size++;
        break;
      }

      for (reco::Vertex::trackRef_iterator it = it_vtx->tracks_begin(); it != it_vtx->tracks_end(); it++) {
        VtxInfo.Pt_Sum[VtxInfo.Size] += (*it)->pt();
        VtxInfo.Pt_Sum2[VtxInfo.Size] += ((*it)->pt() * (*it)->pt());
      }
      if( VtxInfo.Pt_Sum[VtxInfo.Size] >= PVBS_Pt_Max ){
        PVBS_Pt_Max = VtxInfo.Pt_Sum[VtxInfo.Size];
        thePrimaryV = *it_vtx;
      }            
      VtxInfo.Size++;
    }
  }else{ 
    thePrimaryV = Vertex(beamSpot.position(), beamSpot.covariance3D());
  }
  RefVtx = thePrimaryV.position();

  EvtInfo.PVx     = thePrimaryV.position().x();
  EvtInfo.PVy     = thePrimaryV.position().y();
  EvtInfo.PVz     = thePrimaryV.position().z();
  EvtInfo.PVxE    = thePrimaryV.xError();
  EvtInfo.PVyE    = thePrimaryV.yError();
  EvtInfo.PVzE    = thePrimaryV.zError();
  EvtInfo.PVnchi2 = thePrimaryV.normalizedChi2();
  EvtInfo.PVchi2  = thePrimaryV.chi2();

  //printf("-----*****DEBUG:End of EvtInfo.\n");

  // Double check size=0.
  MuonInfo.size   = 0;
  TrackInfo.size  = 0;
  BInfo.uj_size   = 0;
  BInfo.size      = 0;
  GenInfo.size    = 0;
    
  std::vector<int> B_counter;
  for(unsigned int i = 0; i < Bchannel_.size(); i++){
    B_counter.push_back(0);
  }

  // std::vector<pat::Muon>              input_muons;
  auto input_muons = *muons;
  // std::vector<pat::PackedCandidate>   input_tracks;
  // input_tracks = *tks;
  auto input_pc_tracks = *tks; // std::vector<pat::PackedCandidate>
  std::vector<const reco::Track*> input_tracks;
  for(auto tk_it=tks->begin(); tk_it != tks->end(); tk_it++){
      // static const reco::Track* getFromPC(const reco::Candidate& rc) {
    auto rt = getFromPC((*tk_it));
    input_tracks.push_back(rt);
  }

  try{
    const reco::GenParticle* genMuonPtr[MAX_MUON];
    // memset(genMuonPtr,0x00,MAX_MUON);
    memset(genMuonPtr,0x00,MAX_MUON*sizeof(genMuonPtr[0]));
    // const reco::GenParticle* genTrackPtr[MAX_TRACK];
    // memset(genTrackPtr,0x00,MAX_GEN);
    // memset(genTrackPtr,0x00,MAX_TRACK*sizeof(genTrackPtr[0]));
    //standard check for validity of input data
    int genTrackPtr[MAX_TRACK];
    
    if (input_muons.size() == 0){
      if (printInfo_) std::cout << "There's no muon : " << iEvent.id() << std::endl;
    }else{
      if (printInfo_) std::cout << "Got " << input_muons.size() << " muons / ";
      if (input_tracks.size() == 0){
        if (printInfo_) std::cout << "There's no track: " << iEvent.id() << std::endl;
      }else{
        if (printInfo_) std::cout << "Got " << input_tracks.size() << " tracks" << std::endl;
        if (input_tracks.size() > 0 && input_muons.size() > 1){

          //MuonInfo section{{{
          int PassedMuon = 0;
          int mu_hindex = -1;
          for(auto mu_it=input_muons.begin();
              mu_it != input_muons.end() ; mu_it++){
            mu_hindex = int(mu_it - input_muons.begin());
            if(MuonInfo.size >= MAX_MUON){
              fprintf(stderr,"ERROR: number of muons exceeds the size of array.\n");
              break;//exit(0);
            }

            //Muon cut level
            MuonCutLevel->Fill(0);
            if (!(mu_it->isTrackerMuon() || mu_it->isGlobalMuon())) ;
            else {
              MuonCutLevel->Fill(1);
              if (!muon::isGoodMuon(*mu_it,muon::TMOneStationTight)) ;
              else {
                MuonCutLevel->Fill(2);
                if(!mu_it->innerTrack().isNonnull()) ;
                else {
                  MuonCutLevel->Fill(3);
                  if (  fabs(mu_it->innerTrack()->dxy(RefVtx)) >= 3.        || 
                        fabs(mu_it->innerTrack()->dz(RefVtx))  >= 30.       
                        ) ;
                  else {
                    MuonCutLevel->Fill(4);
                    if (mu_it->innerTrack()->hitPattern().pixelLayersWithMeasurement()<1    ||
                        mu_it->innerTrack()->normalizedChi2()>1.8
                        ) ;
                    else {
                      MuonCutLevel->Fill(5);
                      if (mu_it->innerTrack()->hitPattern().trackerLayersWithMeasurement()<6) ;
                      else 
                        MuonCutLevel->Fill(6);
                    }   
                  }
                }
              }
            }

            //Muon Id flag
            // *Bfinder* 
            MuonInfo.BfinderMuID[MuonInfo.size] = false;
            if(mu_it->innerTrack().isNonnull()){
              if( (mu_it->isTrackerMuon() || mu_it->isGlobalMuon()) 
                  // && (muon::isGoodMuon(*mu_it,muon::TMOneStationTight)) 
                  && fabs(mu_it->innerTrack()->dxy(RefVtx)) < 4.
                  && fabs(mu_it->innerTrack()->dz(RefVtx))  < 35.
                  && mu_it->innerTrack()->hitPattern().pixelLayersWithMeasurement() > 0 
                  && mu_it->innerTrack()->hitPattern().trackerLayersWithMeasurement() > 5
                  && mu_it->innerTrack()->normalizedChi2() <= 1.8
                  )
                MuonInfo.BfinderMuID[MuonInfo.size] = true;
            }
            // *SoftMu*
            MuonInfo.SoftMuID[MuonInfo.size] = false;
            if(mu_it->innerTrack().isNonnull()){
              if( (mu_it->isTrackerMuon() && mu_it->isGlobalMuon()) 
                  // && muon::isGoodMuon(*mu_it,muon::TMOneStationTight) 
                  && mu_it->innerTrack()->hitPattern().trackerLayersWithMeasurement() > 5
                  && mu_it->innerTrack()->hitPattern().pixelLayersWithMeasurement() > 0 
                  // && mu_it->innerTrack()->quality(reco::TrackBase::highPurity)
                  && fabs(mu_it->innerTrack()->dxy(RefVtx)) < 0.3
                  && fabs(mu_it->innerTrack()->dz(RefVtx))  < 20.
                  )
                MuonInfo.SoftMuID[MuonInfo.size] = true;
            }

            //outdated selections
            //if(!(mu_it->innerTrack().isNonnull()*mu_it->globalTrack().isNonnull())) {continue;}
            //if (!(mu_it->isGlobalMuon()*mu_it->track().isAvailable()*mu_it->globalTrack().isAvailable())) continue;
            //if (mu_it->p()>200 || mu_it->pt()>200)                  continue;
            //if (mu_it->innerTrack()->hitPattern().trackerLayersWithMeasurement()<6  &&
            //    mu_it->innerTrack()->hitPattern().numberOfValidStripHits()<11
            //   ) continue;

            //Get Muon HLT Trigger matching
            MuonInfo.MuTrgMatchPathSize = MuonTriggerMatchingPath_.size();
            for(int _m = 0; _m < MuonInfo.MuTrgMatchPathSize; _m++){
              pat::TriggerObjectStandAloneCollection match = mu_it->triggerObjectMatchesByPath(MuonTriggerMatchingPath_[_m].c_str());
              if (match.empty()) {
                MuonInfo.MuTrgMatchPathTrgObjE[MuonInfo.size*MuonInfo.MuTrgMatchPathSize+_m] = -999.;
                MuonInfo.MuTrgMatchPathTrgObjPt[MuonInfo.size*MuonInfo.MuTrgMatchPathSize+_m] = -999.;
                MuonInfo.MuTrgMatchPathTrgObjEta[MuonInfo.size*MuonInfo.MuTrgMatchPathSize+_m] = -999.;
                MuonInfo.MuTrgMatchPathTrgObjPhi[MuonInfo.size*MuonInfo.MuTrgMatchPathSize+_m] = -999.;
                //std::cout << "Muon didn't match Trigger Object" << std::endl;
              } else {
                MuonInfo.MuTrgMatchPathTrgObjE[MuonInfo.size*MuonInfo.MuTrgMatchPathSize+_m] = match[0].energy();
                MuonInfo.MuTrgMatchPathTrgObjPt[MuonInfo.size*MuonInfo.MuTrgMatchPathSize+_m] = match[0].pt();
                MuonInfo.MuTrgMatchPathTrgObjEta[MuonInfo.size*MuonInfo.MuTrgMatchPathSize+_m] = match[0].eta();
                MuonInfo.MuTrgMatchPathTrgObjPhi[MuonInfo.size*MuonInfo.MuTrgMatchPathSize+_m] = match[0].phi();
                //std::cout << "Propagation succeeeded; eta = " << match[0].eta() << ", phi = " << match[0].phi() << std::endl;
              }
            }
            MuonInfo.MuTrgMatchFilterSize = MuonTriggerMatchingFilter_.size();
            MuonInfo.isTriggered[MuonInfo.size] = false;
            for(int _m = 0; _m < MuonInfo.MuTrgMatchFilterSize; _m++){
              pat::TriggerObjectStandAloneCollection match = mu_it->triggerObjectMatchesByFilter(MuonTriggerMatchingFilter_[_m].c_str());
              if (match.empty()) {
                MuonInfo.MuTrgMatchFilterTrgObjE[MuonInfo.size*MuonInfo.MuTrgMatchFilterSize+_m] = -999.;
                MuonInfo.MuTrgMatchFilterTrgObjPt[MuonInfo.size*MuonInfo.MuTrgMatchFilterSize+_m] = -999.;
                MuonInfo.MuTrgMatchFilterTrgObjEta[MuonInfo.size*MuonInfo.MuTrgMatchFilterSize+_m] = -999.;
                MuonInfo.MuTrgMatchFilterTrgObjPhi[MuonInfo.size*MuonInfo.MuTrgMatchFilterSize+_m] = -999.;
                //std::cout << "Muon didn't match Trigger Object" << std::endl;
              } else {
                MuonInfo.MuTrgMatchFilterTrgObjE[MuonInfo.size*MuonInfo.MuTrgMatchFilterSize+_m] = match[0].energy();
                MuonInfo.MuTrgMatchFilterTrgObjPt[MuonInfo.size*MuonInfo.MuTrgMatchFilterSize+_m] = match[0].pt();
                MuonInfo.MuTrgMatchFilterTrgObjEta[MuonInfo.size*MuonInfo.MuTrgMatchFilterSize+_m] = match[0].eta();
                MuonInfo.MuTrgMatchFilterTrgObjPhi[MuonInfo.size*MuonInfo.MuTrgMatchFilterSize+_m] = match[0].phi();
                //std::cout << "Propagation succeeeded; eta = " << match[0].eta() << ", phi = " << match[0].phi() << std::endl;
                MuonInfo.isTriggered[MuonInfo.size] = true;
              }
            }
                        
            //Muon general info.
            MuonInfo.index          [MuonInfo.size] = MuonInfo.size;
            MuonInfo.handle_index   [MuonInfo.size] = mu_hindex;
            MuonInfo.charge         [MuonInfo.size] = mu_it->charge();
            MuonInfo.pt             [MuonInfo.size] = mu_it->pt();
            MuonInfo.eta            [MuonInfo.size] = mu_it->eta();
            MuonInfo.phi            [MuonInfo.size] = mu_it->phi();
            MuonInfo.isTrackerMuon  [MuonInfo.size] = mu_it->isTrackerMuon();
            MuonInfo.isGlobalMuon   [MuonInfo.size] = mu_it->isGlobalMuon();
            MuonInfo.iso_trk        [MuonInfo.size] = mu_it->trackIso();//R<0.3
            MuonInfo.iso_ecal       [MuonInfo.size] = mu_it->ecalIso();
            MuonInfo.iso_hcal       [MuonInfo.size] = mu_it->hcalIso();
            MuonInfo.type           [MuonInfo.size] = mu_it->type();//CaloMuon = 1<<4  GlobalMuon = 1<<1  PFMuon = 1<<5  StandAloneMuon = 1<<3  TrackerMuon = 1<<2
            MuonInfo.n_matches      [MuonInfo.size] = mu_it->numberOfMatches();//only in chamber
            MuonInfo.geninfo_index  [MuonInfo.size] = -1;//initialize for later use
            MuonInfo.TMOneStationTight[MuonInfo.size] = muon::isGoodMuon(*mu_it,muon::TMOneStationTight);//For Muon ID for convenience
            MuonInfo.TrackerMuonArbitrated[MuonInfo.size] = muon::isGoodMuon(*mu_it,muon::TrackerMuonArbitrated);//For Muon ID for convenience
            MuonInfo.isSoftMuon[MuonInfo.size] = muon::isSoftMuon(*mu_it, thePrimaryV);
            genMuonPtr              [MuonInfo.size] = 0;
            if (!iEvent.isRealData()) genMuonPtr [MuonInfo.size] = mu_it->genParticle();

            //Muon standalone info.
            MuonInfo.isStandAloneMuon[MuonInfo.size] = false;
            if(mu_it->isStandAloneMuon()){
              MuonInfo.isStandAloneMuon[MuonInfo.size] = true;
              reco::TrackRef tkref;
              tkref = mu_it->standAloneMuon();
              const reco::Track &trk = *tkref;
              MuonInfo.StandAloneMuon_charge         [MuonInfo.size] = trk.charge();
              MuonInfo.StandAloneMuon_pt             [MuonInfo.size] = trk.pt();
              MuonInfo.StandAloneMuon_eta            [MuonInfo.size] = trk.eta();
              MuonInfo.StandAloneMuon_phi            [MuonInfo.size] = trk.phi();
              MuonInfo.StandAloneMuon_d0             [MuonInfo.size] = trk.d0();
              MuonInfo.StandAloneMuon_dz             [MuonInfo.size] = trk.dz();
              MuonInfo.StandAloneMuon_dzPV           [MuonInfo.size] = trk.dz(RefVtx);
              MuonInfo.StandAloneMuon_dxyPV          [MuonInfo.size] = trk.dxy(RefVtx);
            }

            MuonInfo.outerTrackisNonnull[MuonInfo.size] = mu_it->outerTrack().isNonnull();
            MuonInfo.innerTrackisNonnull[MuonInfo.size] = mu_it->innerTrack().isNonnull();
            //Muon inner track info.
            if(mu_it->innerTrack().isNonnull()){
              //Muon inner track track quality
              //enum TrackQuality { undefQuality = -1, loose = 0, tight = 1, highPurity = 2, confirmed = 3, goodIterative = 4, looseSetWithPV = 5, highPuritySetWithPV = 6, qualitySize = 7}
              for(int tq = 0; tq < reco::TrackBase::qualitySize; tq++){
                if (mu_it->innerTrack()->quality(static_cast<reco::TrackBase::TrackQuality>(tq))) MuonInfo.innerTrackQuality[MuonInfo.size] += 1 << (tq);
                //std::cout<<"type: "<<mu_it->innerTrack()->quality(static_cast<reco::TrackBase::TrackQuality>(tq))<<std::endl;
              }
              MuonInfo.highPurity              [MuonInfo.size] = mu_it->innerTrack()->quality(reco::TrackBase::highPurity);
              MuonInfo.normchi2                [MuonInfo.size] = mu_it->innerTrack()->normalizedChi2();
              MuonInfo.i_striphit              [MuonInfo.size] = mu_it->innerTrack()->hitPattern().numberOfValidStripHits();
              MuonInfo.i_pixelhit              [MuonInfo.size] = mu_it->innerTrack()->hitPattern().numberOfValidPixelHits();
              MuonInfo.i_nStripLayer           [MuonInfo.size] = mu_it->innerTrack()->hitPattern().stripLayersWithMeasurement();
              MuonInfo.i_nPixelLayer           [MuonInfo.size] = mu_it->innerTrack()->hitPattern().pixelLayersWithMeasurement();
              MuonInfo.i_chi2                  [MuonInfo.size] = mu_it->innerTrack()->chi2();
              MuonInfo.i_ndf                   [MuonInfo.size] = mu_it->innerTrack()->ndof();
              //MuonInfo.fpbarrelhit             [MuonInfo.size] = mu_it->innerTrack()->hitPattern().hasValidHitInFirstPixelBarrel();
              //MuonInfo.fpendcaphit             [MuonInfo.size] = mu_it->innerTrack()->hitPattern().hasValidHitInFirstPixelEndcap();
              //https://github.com/cms-sw/cmssw/blob/CMSSW_9_2_3/DataFormats/TrackReco/src/HitPattern.cc#L321
              //https://github.com/cms-sw/cmssw/blob/CMSSW_9_2_3/DataFormats/SiPixelDetId/interface/PixelSubdetector.h#L11
              MuonInfo.fpbarrelhit             [MuonInfo.size] = mu_it->innerTrack()->hitPattern().hasValidHitInPixelLayer(PixelSubdetector::PixelBarrel,1);
              MuonInfo.fpendcaphit             [MuonInfo.size] = mu_it->innerTrack()->hitPattern().hasValidHitInPixelLayer(PixelSubdetector::PixelEndcap,1);
              MuonInfo.ptErr                   [MuonInfo.size] = mu_it->track()->ptError();
              MuonInfo.etaErr                  [MuonInfo.size] = mu_it->track()->etaError();
              MuonInfo.phiErr                  [MuonInfo.size] = mu_it->track()->phiError();
              MuonInfo.d0                      [MuonInfo.size] = mu_it->track()->d0();
              MuonInfo.dz                      [MuonInfo.size] = mu_it->track()->dz();
              MuonInfo.dzPV                    [MuonInfo.size] = mu_it->track()->dz(RefVtx);//==mu_it->innerTrack()->dxy(thePrimaryV.position());
              MuonInfo.dxyPV                   [MuonInfo.size] = mu_it->track()->dxy(RefVtx);//==mu_it->innerTrack()->dz(thePrimaryV.position());
              //mu_it->innerTrack()->hitPattern().trackerLayersWithMeasurement() == MuonInfo.i_nStripLayer + MuonInfo.i_nPixelLayer
            }
            //Muon global track info.
            MuonInfo.globalTrackisNonnull[MuonInfo.size] = mu_it->globalTrack().isNonnull();
            if(mu_it->isGlobalMuon()){
              MuonInfo.g_striphit [MuonInfo.size] = mu_it->globalTrack()->hitPattern().numberOfValidStripHits();
              MuonInfo.g_pixelhit [MuonInfo.size] = mu_it->globalTrack()->hitPattern().numberOfValidPixelHits();
              MuonInfo.g_chi2     [MuonInfo.size] = mu_it->globalTrack()->chi2();
              MuonInfo.g_ndf      [MuonInfo.size] = mu_it->globalTrack()->ndof();
              MuonInfo.nmuhit     [MuonInfo.size] = mu_it->globalTrack()->hitPattern().numberOfValidMuonHits();
            }else{
              MuonInfo.g_striphit [MuonInfo.size] = -1;
              MuonInfo.g_pixelhit [MuonInfo.size] = -1;
              MuonInfo.g_chi2     [MuonInfo.size] = -1;
              MuonInfo.g_ndf      [MuonInfo.size] = -1;
              MuonInfo.nmuhit     [MuonInfo.size] = -1;
            }
            //Muon quality
            //https://twiki.cern.ch/twiki/bin/view/CMSPublic/WorkBookMuonAnalysis
            int qm = 0;
            for(int qi=1; qi!= 24; ++qi){
              if (muon::isGoodMuon(*mu_it, muon::SelectionType(qi))){
                qm += 1 << qi;
              }
            }
            MuonInfo.muqual         [MuonInfo.size] = qm;   

            // Basic Muon selection for fitting
            // Muon in acceptance
            //bool MuInAcc = false;
            //if(fabs(mu_it->eta())<1.3 && mu_it->pt()>3.3) MuInAcc = true;
            //if(fabs(mu_it->eta())>1.3 && fabs(mu_it->eta())<2.2 && mu_it->p()>2.9) MuInAcc = true;
            //if(fabs(mu_it->eta())>2.2 && fabs(mu_it->eta())<2.4 && mu_it->pt()>0.8) MuInAcc = true;

            //Can not be just CaloMuon or empty type
            if((MuonInfo.type[MuonInfo.size]|(1<<4))==(1<<4)){
              MuonInfo.isNeededMuon[MuonInfo.size] = false;
            }
            else if(doMuPreCut_ &&
                    (  // !muon::isGoodMuon(*mu_it,muon::TMOneStationTight)
                     //|| !MuInAcc
                     //||  some other cut
                     !(mu_it->isTrackerMuon() || mu_it->isGlobalMuon())
                       )
                    ){
              MuonInfo.isNeededMuon[MuonInfo.size] = false;
            }
            else {
              PassedMuon ++;
              MuonInfo.isNeededMuon[MuonInfo.size] = true;
            }

            MuonInfo.size++;
          }//end of MuonInfo}}}
          //std::cout<<"PassedMuon: "<<PassedMuon<<std::endl;
          // printf("-----*****DEBUG:End of MuonInfo.\n");

          //Preselect tracks{{{
          std::vector<bool> isNeededTrack;// Are the tracks redundant?
          int PassedTrk = 0;
          // for(std::vector<pat::PackedCandidate>::const_iterator tk_it=input_tracks.begin();
          //     tk_it != input_tracks.end(); tk_it++){

          for(auto tk_it_it=input_tracks.begin(); // edm::View<pat::PackedCandidate>::const_iterator tk_it
              tk_it_it != input_tracks.end(); tk_it_it++){
            if(PassedTrk >= MAX_TRACK){
              fprintf(stderr,"ERROR: number of tracks exceeds the size of array.\n");
              break;
            }
            std::cout<<"HERE1 " <<std::endl;
            std::cout << "tk_it_it index: " << std::distance(input_tracks.begin(), tk_it_it) << std::endl;  
            auto tk_it = (*tk_it_it); // edm::View<pat::PackedCandidate>
            std::cout << "tk_it " << tk_it << std::endl;

            isNeededTrack.push_back(false);
            if (!tk_it) continue;
            std::cout<<"HERE2 " <<std::endl;
            // if(!tk_it->hasTrackDetails()) continue;
            // TrackCutLevel->Fill(0);//number of all tracks
            bool isMuonTrack = false; //remove muon track
            for(auto it=input_muons.begin() ; // std::vector<pat::Muon>::iterator
                it != input_muons.end(); it++){
              if (!it->track().isNonnull())                   continue;
              if((it->type()|(1<<4))==(1<<4)) continue;//Don't clean track w.r.t. calo muon 
              if (fabs(tk_it->pt() -it->track()->pt() )<0.00001 &&
                  fabs(tk_it->eta()-it->track()->eta())<0.00001 &&
                  fabs(tk_it->phi()-it->track()->phi())<0.00001 ){
                isMuonTrack = true;
                break;
              }
            }
            if (isMuonTrack)                                    continue;
            if (abs(tk_it->charge()) != 1) continue;
            if (tk_it->pt()<tkPtCut_)                           continue;
            if (fabs(tk_it->eta())>tkEtaCut_)                   continue;
            if(doTkPreCut_){
              if( !(tk_it->/*pseudoTrack().*/quality(reco::TrackBase::qualityByName("highPurity")))){ 
                continue;}
            }
            isNeededTrack[tk_it_it-input_tracks.begin()] = true;
            PassedTrk++;
          }//end of track preselection}}}
          if(printInfo_) std::cout<<"PassedTrk: "<<PassedTrk<<std::endl;                    
          //printf("-----*****DEBUG:End of track preselection.\n");

          // BInfo section{{{
          int mu1_index = -1;
          int mu1_hindex = -1;
          bool gogogo = false;
          for(auto mu_it1=input_muons.begin(); // std::vector<pat::Muon>::const_iterator
              mu_it1 != input_muons.end(); mu_it1++){
            //check if muon track is non null
            if(!mu_it1->track().isNonnull()) continue;
            //Check if it in MuonInfo and isNeedeMuon
            mu1_hindex = int(mu_it1 - input_muons.begin());
            gogogo = false;
            for(int i=0; i < MuonInfo.size; i++){
              if (mu1_hindex == MuonInfo.handle_index[i] && MuonInfo.isNeededMuon[i]){
                gogogo = true;
                break;
              }
            }
            if (!gogogo) continue;
            //Get the corrisponding index in MuonInfo
            mu1_index ++;
            if (mu_it1->charge()<0) continue;

            int mu2_index = -1;
            int mu2_hindex = -1; 
            for(auto mu_it2=input_muons.begin();
                mu_it2 != input_muons.end(); mu_it2++){
              //check if muon track is non null
              if(!mu_it2->track().isNonnull()) continue;
              //Check if it in MuonInfo and isNeedeMuon
              mu2_hindex = int(mu_it2 - input_muons.begin()); 
              gogogo = false;
              for(int j=0; j < MuonInfo.size; j++){
                if(mu2_hindex == MuonInfo.handle_index[j] && MuonInfo.isNeededMuon[j]){
                  gogogo = true;
                  break;
                }
              }
              if (!gogogo) continue;
              mu2_index ++;
              if (mu_it2->charge()>0) continue;
              // XbujCutLevel->Fill(0);

              TLorentzVector v4_mu1,v4_mu2;
              v4_mu1.SetPtEtaPhiM(mu_it1->pt(),mu_it1->eta(),mu_it1->phi(),MUON_MASS);
              v4_mu2.SetPtEtaPhiM(mu_it2->pt(),mu_it2->eta(),mu_it2->phi(),MUON_MASS);
              if (fabs((v4_mu1+v4_mu2).Mag()-JPSI_MASS)>0.4) continue;
              if((v4_mu1+v4_mu2).Pt()<jpsiPtCut_)continue;

              //Fit 2 muon
              reco::TransientTrack muonPTT(mu_it1->track(), &(*bField) );
              reco::TransientTrack muonMTT(mu_it2->track(), &(*bField) );
              if(!muonPTT.isValid()) continue;
              if(!muonMTT.isValid()) continue;
              // XbujCutLevel->Fill(1);
    
              // const reco::Muon* rmu1 = dynamic_cast<const reco::Muon * >(mu_it1->originalObject());
              // const reco::Muon* rmu2 = dynamic_cast<const reco::Muon * >(mu_it2->originalObject());
              // std::cout<<rmu1<<", "<<rmu2<<std::endl;
              // if(muon::overlap(*rmu1, *rmu2)) continue;
              // XbujCutLevel->Fill(2);
    
              KinematicParticleFactoryFromTransientTrack pFactory;
              ParticleMass muon_mass = MUON_MASS; //pdg mass
              float muon_sigma = muon_mass*1.e-6;
              float chi = 0.;
              float ndf = 0.;
              std::vector<RefCountedKinematicParticle> muonParticles;
              muonParticles.push_back(pFactory.particle(muonPTT,muon_mass,chi,ndf,muon_sigma));
              muonParticles.push_back(pFactory.particle(muonMTT,muon_mass,chi,ndf,muon_sigma));
    
              KinematicParticleVertexFitter   fitter;   
              RefCountedKinematicTree         ujVFT;
              ujVFT = fitter.fit(muonParticles); 
              if (!ujVFT->isValid()) continue;
              XbujCutLevel->Fill(3); 

              ujVFT->movePointerToTheTop();
    
              RefCountedKinematicParticle ujVFP       = ujVFT->currentParticle();
              RefCountedKinematicVertex   ujVFPvtx    = ujVFT->currentDecayVertex();
              ujVFT->movePointerToTheFirstChild();
              KinematicParameters         ujmu1KP     = ujVFT->currentParticle()->currentState().kinematicParameters();
              ujVFT->movePointerToTheNextChild();
              KinematicParameters         ujmu2KP     = ujVFT->currentParticle()->currentState().kinematicParameters();
              double chi2_prob_uj = TMath::Prob(ujVFPvtx->chiSquared(), ujVFPvtx->degreesOfFreedom());
              if(chi2_prob_uj < uj_VtxChiProbCut_) continue;
              XbujCutLevel->Fill(4);

              if (fabs(ujVFP->currentState().mass()-JPSI_MASS)>0.3) continue;

              TLorentzVector uj_4vec,uj_mu1_4vec,uj_mu2_4vec;
              uj_4vec.SetPxPyPzE(ujVFP->currentState().kinematicParameters().momentum().x(),
                                 ujVFP->currentState().kinematicParameters().momentum().y(),
                                 ujVFP->currentState().kinematicParameters().momentum().z(),
                                 ujVFP->currentState().kinematicParameters().energy());
              uj_mu1_4vec.SetPxPyPzE( ujmu1KP.momentum().x(),
                                      ujmu1KP.momentum().y(),
                                      ujmu1KP.momentum().z(),
                                      ujmu1KP.energy());
              uj_mu2_4vec.SetPxPyPzE( ujmu2KP.momentum().x(),
                                      ujmu2KP.momentum().y(),
                                      ujmu2KP.momentum().z(),
                                      ujmu2KP.energy());
              //uj_4vec.Print();

              BInfo.uj_index         [BInfo.uj_size]= BInfo.uj_size;
              BInfo.uj_mass          [BInfo.uj_size]= uj_4vec.Mag();
              BInfo.uj_pt            [BInfo.uj_size]= uj_4vec.Pt();
              BInfo.uj_eta            [BInfo.uj_size]= uj_4vec.Eta();
              BInfo.uj_phi            [BInfo.uj_size]= uj_4vec.Phi();
              BInfo.uj_px            [BInfo.uj_size]= uj_4vec.Px();
              BInfo.uj_py            [BInfo.uj_size]= uj_4vec.Py();
              BInfo.uj_pz            [BInfo.uj_size]= uj_4vec.Pz();
              BInfo.uj_vtxX          [BInfo.uj_size]= ujVFPvtx->position().x();
              BInfo.uj_vtxY          [BInfo.uj_size]= ujVFPvtx->position().y();
              BInfo.uj_vtxZ          [BInfo.uj_size]= ujVFPvtx->position().z();
              BInfo.uj_vtxXErr       [BInfo.uj_size]= ujVFPvtx->error().cxx();
              BInfo.uj_vtxYErr       [BInfo.uj_size]= ujVFPvtx->error().cyy();
              BInfo.uj_vtxZErr       [BInfo.uj_size]= ujVFPvtx->error().czz();
              BInfo.uj_vtxYXErr      [BInfo.uj_size]= ujVFPvtx->error().cyx();
              BInfo.uj_vtxZXErr      [BInfo.uj_size]= ujVFPvtx->error().czx();
              BInfo.uj_vtxZYErr      [BInfo.uj_size]= ujVFPvtx->error().czy();
              BInfo.uj_vtxdof        [BInfo.uj_size]= ujVFPvtx->degreesOfFreedom();
              BInfo.uj_vtxchi2       [BInfo.uj_size]= ujVFPvtx->chiSquared();
              BInfo.uj_rfmu1_index   [BInfo.uj_size]= mu1_hindex;
              BInfo.uj_rfmu2_index   [BInfo.uj_size]= mu2_hindex;
              BInfo.uj_rfmu1_pt      [BInfo.uj_size]= uj_mu1_4vec.Pt();
              BInfo.uj_rfmu1_eta     [BInfo.uj_size]= uj_mu1_4vec.Eta();
              BInfo.uj_rfmu1_phi     [BInfo.uj_size]= uj_mu1_4vec.Phi();
              BInfo.uj_rfmu2_pt      [BInfo.uj_size]= uj_mu2_4vec.Pt();
              BInfo.uj_rfmu2_eta     [BInfo.uj_size]= uj_mu2_4vec.Eta();
              BInfo.uj_rfmu2_phi     [BInfo.uj_size]= uj_mu2_4vec.Phi();

              BInfo.uj_size++;
              muonParticles.clear();

              //////////////////////////////////////////////////////////////////////////
              // RECONSTRUCTION: J/psi + K
              //////////////////////////////////////////////////////////////////////////
              //float mass_window[2] = {4.3, 6.4};
              //float mass_window[2] = {5., 6.};
              float mass_window[2] = {4.5, 6.5};
              if(Bchannel_[0] == 1){
                BranchOut2MuTk(
                               BInfo,
                               input_tracks,
                               thePrimaryV,
                               isNeededTrack,
                               v4_mu1,
                               v4_mu2,
                               muonPTT,
                               muonMTT,
                               B_counter,
                               mass_window,
                               JPSI_MASS,
                               KAON_MASS,
                               1
                               );
              }                            
              //////////////////////////////////////////////////////////////////////////
              // RECONSTRUCTION: J/psi + Pi
              //////////////////////////////////////////////////////////////////////////
              if(Bchannel_[1] == 1){
                BranchOut2MuTk(
                               BInfo,
                               input_tracks,
                               thePrimaryV,
                               isNeededTrack,
                               v4_mu1,
                               v4_mu2,
                               muonPTT,
                               muonMTT,
                               B_counter,
                               mass_window,
                               JPSI_MASS,
                               PION_MASS,
                               2
                               );
              }                            

              //////////////////////////////////////////////////////////////////////////
              // RECONSTRUCTION: J/psi + Ks
              //////////////////////////////////////////////////////////////////////////

              float TkTk_window = 0;
              TkTk_window = 0.3;
              if(Bchannel_[2] == 1){
                BranchOut2MuX_XtoTkTk(
                                      BInfo,
                                      input_tracks,
                                      thePrimaryV,
                                      isNeededTrack,
                                      v4_mu1,
                                      v4_mu2,
                                      muonPTT,
                                      muonMTT,
                                      B_counter,
                                      mass_window,
                                      JPSI_MASS,
                                      KSHORT_MASS,
                                      TkTk_window,
                                      PION_MASS,        
                                      PION_MASS,
                                      3,
                                      1
                                      );
              }                            
                            
              //////////////////////////////////////////////////////////////////////////
              // RECONSTRUCTION: J/psi + K* (K+, Pi-)
              //////////////////////////////////////////////////////////////////////////

              //TkTk_window = 0.4;
              TkTk_window = 0.25;
              if(Bchannel_[3] == 1){
                BranchOut2MuX_XtoTkTk(
                                      BInfo,
                                      input_tracks,
                                      thePrimaryV,
                                      isNeededTrack,
                                      v4_mu1,
                                      v4_mu2,
                                      muonPTT,
                                      muonMTT,
                                      B_counter,
                                      mass_window,
                                      JPSI_MASS,
                                      KSTAR_MASS,
                                      TkTk_window,
                                      KAON_MASS,        
                                      PION_MASS,
                                      4,
                                      0
                                      );
              }
                            
              //////////////////////////////////////////////////////////////////////////
              // RECONSTRUCTION: J/psi + K* (K-, Pi+)
              //////////////////////////////////////////////////////////////////////////

              //TkTk_window = 0.4;
              TkTk_window = 0.25;
              if(Bchannel_[4] == 1){
                BranchOut2MuX_XtoTkTk(
                                      BInfo,
                                      input_tracks,
                                      thePrimaryV,
                                      isNeededTrack,
                                      v4_mu1,
                                      v4_mu2,
                                      muonPTT,
                                      muonMTT,
                                      B_counter,
                                      mass_window,
                                      JPSI_MASS,
                                      KSTAR_MASS,
                                      TkTk_window,
                                      PION_MASS,        
                                      KAON_MASS,
                                      5,
                                      0
                                      );
              }
                            
              //////////////////////////////////////////////////////////////////////////
              // RECONSTRUCTION: J/psi + phi
              //////////////////////////////////////////////////////////////////////////
                            
              TkTk_window = 0.15;
              if(Bchannel_[5] == 1){
                BranchOut2MuX_XtoTkTk(
                                      BInfo,
                                      input_tracks,
                                      thePrimaryV,
                                      isNeededTrack,
                                      v4_mu1,
                                      v4_mu2,
                                      muonPTT,
                                      muonMTT,
                                      B_counter,
                                      mass_window,
                                      JPSI_MASS,
                                      PHI_MASS,
                                      TkTk_window,
                                      KAON_MASS,        
                                      KAON_MASS,
                                      6,
                                      0
                                      );
              }

              //////////////////////////////////////////////////////////////////////////
              // RECONSTRUCTION: J/psi + pi pi <= psi', X(3872), Bs->J/psi f0
              //////////////////////////////////////////////////////////////////////////
              //mass_window[0] = 3;
              //mass_window[1] = 6.4;
              //TkTk_window = 1.6;
              //mass_window[0] = 3.4;
              //mass_window[1] = 4.2;
              mass_window[0] = 2.5;
              mass_window[1] = 4.5;
              TkTk_window = 0;
              if(Bchannel_[6] == 1){
                BranchOut2MuX_XtoTkTk(
                                      BInfo,
                                      input_tracks,
                                      thePrimaryV,
                                      isNeededTrack,
                                      v4_mu1,
                                      v4_mu2,
                                      muonPTT,
                                      muonMTT,
                                      B_counter,
                                      mass_window,
                                      JPSI_MASS,
                                      -1,
                                      TkTk_window,
                                      PION_MASS,        
                                      PION_MASS,
                                      7,
                                      0
                                      );
              }
                            
            }//Mu2
          }//Mu1
          if(printInfo_){
            printf("B_counter: ");
            for(unsigned int i = 0; i < Bchannel_.size(); i++){
              printf("%d/", B_counter[i]);
            }
            printf("\n");
          }//}}}
          //printf("-----*****DEBUG:End of BInfo.\n");

          // TrackInfo section {{{
          // Setup MVA
          /* Under construction !! */

          // Setup Dedx
          /* Under construction !! */

          for(auto tk_it_it=input_tracks.begin();
              tk_it_it != input_tracks.end() ; tk_it_it++){
            int tk_hindex = int(tk_it_it - input_tracks.begin());
            auto tk_it = (*tk_it_it);
            if(TrackInfo.size >= MAX_TRACK){
              fprintf(stderr,"ERROR: number of tracks exceeds the size of array.\n");
              break;
            }

            if(tk_hindex>=int(isNeededTrack.size())) break;
            if (isNeededTrack[tk_hindex]==false) continue;

            //Create list of relative xb candidates for later filling
            std::vector<int> listOfRelativeXbCands1;//1~nXb
            std::vector<int> listOfRelativeXbCands2;//1~nXb
            for(int iXb=0; iXb < BInfo.size; iXb++){
              if(BInfo.rftk1_index[iXb] == tk_hindex){
                listOfRelativeXbCands1.push_back(iXb+1);
              }if(BInfo.rftk2_index[iXb] == tk_hindex){
                listOfRelativeXbCands2.push_back(iXb+1);
              }
            }
            if(dropUnusedTracks_ && listOfRelativeXbCands1.size() == 0 && listOfRelativeXbCands2.size() == 0) continue;//drop unused tracks
                        
            TrackInfo.index          [TrackInfo.size] = TrackInfo.size;
            TrackInfo.handle_index   [TrackInfo.size] = tk_hindex;
            TrackInfo.charge         [TrackInfo.size] = tk_it->charge();
            TrackInfo.pt             [TrackInfo.size] = tk_it->pt();
            TrackInfo.eta            [TrackInfo.size] = tk_it->eta();
            TrackInfo.phi            [TrackInfo.size] = tk_it->phi();
            TrackInfo.ptErr          [TrackInfo.size] = tk_it->/*pseudoTrack().*/ptError();
            TrackInfo.etaErr         [TrackInfo.size] = tk_it->/*pseudoTrack().*/etaError();
            TrackInfo.phiErr         [TrackInfo.size] = tk_it->/*pseudoTrack().*/phiError();
            //TrackInfo.p              [TrackInfo.size] = tk_it->p();
            TrackInfo.striphit       [TrackInfo.size] = tk_it->/*pseudoTrack().*/hitPattern().numberOfValidStripHits();
            TrackInfo.pixelhit       [TrackInfo.size] = tk_it->/*pseudoTrack().*/hitPattern().numberOfValidPixelHits();
            TrackInfo.nStripLayer    [TrackInfo.size] = tk_it->/*pseudoTrack().*/hitPattern().stripLayersWithMeasurement();
            TrackInfo.nPixelLayer    [TrackInfo.size] = tk_it->/*pseudoTrack().*/hitPattern().pixelLayersWithMeasurement();
            TrackInfo.fpbarrelhit    [TrackInfo.size] = tk_it->/*pseudoTrack().*/hitPattern().hasValidHitInPixelLayer(PixelSubdetector::PixelBarrel,1);
            TrackInfo.fpendcaphit    [TrackInfo.size] = tk_it->/*pseudoTrack().*/hitPattern().hasValidHitInPixelLayer(PixelSubdetector::PixelEndcap,1);
            TrackInfo.chi2           [TrackInfo.size] = tk_it->/*pseudoTrack().*/chi2();
            if (chi2Handle.isValid() && !chi2Handle.failedToGet())
              TrackInfo.chi2           [TrackInfo.size] = (float)((*chi2Handle)[ tks->ptrAt( tk_hindex ) ]) * tk_it->/*pseudoTrack().*/ndof();
            TrackInfo.ndf            [TrackInfo.size] = tk_it->/*pseudoTrack().*/ndof();
            TrackInfo.d0             [TrackInfo.size] = tk_it->/*pseudoTrack().*/d0();
            TrackInfo.d0error        [TrackInfo.size] = tk_it->/*pseudoTrack().*/d0Error();
            TrackInfo.dz             [TrackInfo.size] = tk_it->/*pseudoTrack().*/dz();
            TrackInfo.dzerror        [TrackInfo.size] = tk_it->/*pseudoTrack().*/dzError();
            TrackInfo.dxy            [TrackInfo.size] = tk_it->/*pseudoTrack().*/dxy();
            TrackInfo.dxyerror       [TrackInfo.size] = tk_it->/*pseudoTrack().*/dxyError();
            TrackInfo.dxy1           [TrackInfo.size] = tk_it->/*pseudoTrack().*/dxy(RefVtx);
            TrackInfo.dxyerror1      [TrackInfo.size] = TMath::Sqrt(tk_it->/*pseudoTrack().*/dxyError()*tk_it->/*pseudoTrack().*/dxyError() + thePrimaryV.xError()*thePrimaryV.yError());
            TrackInfo.dz1            [TrackInfo.size] = tk_it->/*pseudoTrack().*/dz(RefVtx);
            TrackInfo.dzerror1       [TrackInfo.size] = TMath::Sqrt(tk_it->/*pseudoTrack().*/dzError()*tk_it->/*pseudoTrack().*/dzError() + thePrimaryV.zError()*thePrimaryV.zError());
            TrackInfo.highPurity     [TrackInfo.size] = tk_it->/*pseudoTrack().*/quality(reco::TrackBase::qualityByName("highPurity"));
            TrackInfo.geninfo_index  [TrackInfo.size] = -1;//initialize for later use
            TrackInfo.geninfo_pdgId  [TrackInfo.size] = -1; // initialize for later use
            // TrackInfo.trkMVAVal      [TrackInfo.size] = (*mvaoutput)[tk_it->track()];
            TrackInfo.trkAlgo        [TrackInfo.size] = tk_it->/*pseudoTrack().*/algo();
            TrackInfo.originalTrkAlgo[TrackInfo.size] = tk_it->/*pseudoTrack().*/originalAlgo();
            // if(readDedx_) {
            //     TrackInfo.dedx           [TrackInfo.size] = dEdxTrack1[tk_it->track()].dEdx();
            // }else
            TrackInfo.dedx           [TrackInfo.size] = -1;
            if (dedxHandle.isValid() && !dedxHandle.failedToGet())
              TrackInfo.dedx           [TrackInfo.size] = ((*dedxHandle)[ tks->ptrAt( tk_hindex ) ]).dEdx();
            
            // Gen-match
            // https://github.com/cms-sw/cmssw/blob/CMSSW_11_2_X/CommonTools/UtilAlgos/interface/MatchByDRDPt.h
            if (!iEvent.isRealData())
              {
                genTrackPtr[TrackInfo.size] = -1;
                float currentdR = 1.e+10;//, currentptRel=0.;
                //bool isGenSignal = false;
                reco::GenParticle _deRef;
                for(auto it_gen=gens->begin();
                    it_gen != gens->end(); it_gen++)
                  {
                    int abspdg = abs(int(it_gen->pdgId()));
                    if (it_gen->status() != 1) continue;
                    if (abspdg != 211 &&
                        abspdg != 321 &&
                        abspdg != 2212) continue;
                    if (tk_it->charge() != it_gen->charge()) continue;
                    float ptRel = fabs(tk_it->pt() - it_gen->pt()) / tk_it->pt();
                    if(ptRel >= 0.2) continue; //
                    float deta = tk_it->eta() - it_gen->eta();
                    float dphi = std::abs(tk_it->phi() - it_gen->phi());
                    if(dphi > float(M_PI))
                      dphi -= float(2 * M_PI);
                    float dR = sqrt(deta*deta + dphi*dphi);
                    if(dR >= 0.02) continue; //
                    if(dR < currentdR) {
                      genTrackPtr[TrackInfo.size] = int(it_gen - gens->begin());
                      currentdR = dR;

                      _deRef = (*it_gen);
                      //reco::Candidate* Myself = dynamic_cast<reco::Candidate*>(&_deRef);
                      //isGenSignal = (Functs.GetAncestor(Myself, 5) | Functs.GetAncestor(Myself, 4));
                    }
                  }
                // genTrackPtr[TrackInfo.size] = tk_it->genParticle();                
                /*if (genTrackPtr[TrackInfo.size] >= 0)
                  
                  std::cout<<"\e[32m"<<TrackInfo.size<<": "<<genTrackPtr[TrackInfo.size]
                           <<" -> "<<isGenSignal
                           <<"\e[0m"<<std::endl;*/
              }
            // <--------------

            //Fill correct track index and track quality to correspond Xb candidate
            for(unsigned int iCands=0; iCands < listOfRelativeXbCands1.size(); iCands++){
              BInfo.rftk1_index[listOfRelativeXbCands1[iCands]-1] = TrackInfo.size;
            }
            for(unsigned int iCands=0; iCands < listOfRelativeXbCands2.size(); iCands++){
              BInfo.rftk2_index[listOfRelativeXbCands2[iCands]-1] = TrackInfo.size;
            }

            TrackInfo.size++;
          }//end of TrackInfo}}}
          //printf("-----*****DEBUG:End of TrackInfo.\n");
        }//has nTracks>1
      }//if no Tracks
    }//if no Muons

    // GenInfo section{{{
    if (!iEvent.isRealData()){


      std::vector<const reco::Candidate *> sel_cands;
      for(auto it_gen=gens->begin(); // std::vector<reco::GenParticle>::const_iterator
          it_gen != gens->end(); it_gen++){
        // if (it_gen->status() > 2 && it_gen->status() != 8 && it_gen->status() != 91 && abs(it_gen->pdgId())!=13) continue; //only status 1, 2, 8(simulated), save all muons
        if (it_gen->status() > 2 && it_gen->status() != 8 && abs(it_gen->pdgId())!=13) continue; //only status 1, 2, 8(simulated), save all muons
        if(GenInfo.size >= MAX_GEN){
          fprintf(stderr,"ERROR: number of gens exceeds the size of array.\n");
          break;;
        }
                
        bool isGenSignal = false;
        //save target intermediate state particle
        if (
            abs(int(it_gen->pdgId()/100) % 100) == 4  ||//c menson
            abs(int(it_gen->pdgId()/100) % 100) == 5  ||//b menson
            //abs(it_gen->pdgId()) == 511 ||//B_0
            //abs(it_gen->pdgId()) == 521 ||//B_+-
            //abs(it_gen->pdgId()) == 531 ||//B_s
            //abs(it_gen->pdgId()) == 311 ||//K0
            //abs(it_gen->pdgId()) == 321 ||//K+
            //abs(it_gen->pdgId()) == 310 ||//KS
            //abs(it_gen->pdgId()) == 313 ||//K*0(892)
            //abs(it_gen->pdgId()) == 323 ||//K*+-(892)
            //abs(it_gen->pdgId()) == 333 ||//phi(1020)
            it_gen->pdgId() == 443     ||//Jpsi
            it_gen->pdgId() == 100443  ||//Psi(2S)
            it_gen->pdgId() == 20443   ||//chi_c1(1P)
            it_gen->pdgId() == 9920443 ||//X3872
            it_gen->pdgId() == 553     ||//Upsilon
            it_gen->pdgId() == 100553    //Upsilon(2S)
            ) isGenSignal = true; //b, c

        if (abs(it_gen->pdgId()) == 13) isGenSignal = true;//all mu

        if (
            abs(int(it_gen->pdgId()/100) % 100) == 3  ||//s menson
            abs(it_gen->pdgId()) == 111 ||//pi0
            abs(it_gen->pdgId()) == 113 ||//rho0
            abs(it_gen->pdgId()) == 130 ||//KL
            abs(it_gen->pdgId()) == 211   //pi+
            ){
          reco::GenParticle _deRef = (*it_gen);
          reco::Candidate* Myself = dynamic_cast<reco::Candidate*>(&_deRef);
          //std::cout<<Myself->pdgId()<<"-----------"<<std::endl;
          isGenSignal = (Functs.GetAncestor(Myself, 5) | Functs.GetAncestor(Myself, 4));
        }//all pi and K from b or c meson
        if (!isGenSignal) continue;

        for(int muonIdx = 0; muonIdx < MuonInfo.size; muonIdx++){
          // match by pat::Muon
          if (genMuonPtr[muonIdx] == 0) continue;
          if (it_gen->p4() == genMuonPtr[muonIdx]->p4()){
            // std::cout<<it_gen<<" =?= "<<genMuonPtr[muonIdx]<<std::endl;
            MuonInfo.geninfo_index[muonIdx] = GenInfo.size;
            break;
          }
        }

        for(int trackIdx = 0; trackIdx < TrackInfo.size; trackIdx++){ //# saving gen-track ref index
          // hope for match by pat::GenericParticle
          if (genTrackPtr[trackIdx] < 0 ) continue;
          if (int(it_gen - gens->begin()) == genTrackPtr[trackIdx]) {
            //std::cout<<"\e[33m"<<trackIdx<<": "<<genTrackPtr[trackIdx]<<"\e[0m"<<std::endl;
            TrackInfo.geninfo_index[trackIdx] = GenInfo.size;
            TrackInfo.geninfo_pdgId[trackIdx] = it_gen->pdgId();
            // break;
          }
        }

        GenInfo.index[GenInfo.size]         = GenInfo.size;
        GenInfo.handle_index[GenInfo.size]  = it_gen-gens->begin();
        GenInfo.pt[GenInfo.size]            = it_gen->pt();
        GenInfo.eta[GenInfo.size]           = it_gen->eta();
        GenInfo.phi[GenInfo.size]           = it_gen->phi();
        GenInfo.mass[GenInfo.size]          = it_gen->mass();
        GenInfo.pdgId[GenInfo.size]         = it_gen->pdgId();
        GenInfo.status[GenInfo.size]        = it_gen->status();
        GenInfo.collisionId[GenInfo.size]   = it_gen->collisionId();
        GenInfo.nMo[GenInfo.size]           = it_gen->numberOfMothers();
        GenInfo.nDa[GenInfo.size]           = it_gen->numberOfDaughters();
        GenInfo.mo1[GenInfo.size]           = -1;//To be matched later.
        GenInfo.mo2[GenInfo.size]           = -1;
        GenInfo.da1[GenInfo.size]           = -1;
        GenInfo.da2[GenInfo.size]           = -1;
        GenInfo.da3[GenInfo.size]           = -1;
        GenInfo.da4[GenInfo.size]           = -1;
        GenInfo.size++;
        sel_cands.push_back(&*it_gen);
      }
      //printf("-----*****DEBUG:End of gens loop.\n");

      int geninfo_idx = 0;
      for(std::vector<const reco::Candidate *>::iterator sCands = sel_cands.begin();
          sCands != sel_cands.end(); sCands++){
        geninfo_idx = int(sCands-sel_cands.begin());
        for(int nGenMo = 0; nGenMo < std::min(2,int((*sCands)->numberOfMothers())); nGenMo++){
          //if((*sCands)->numberOfMothers()==1){
          for(std::vector<const reco::Candidate *>::iterator mCands = sel_cands.begin();
              mCands != sel_cands.end(); mCands++){
            if((*sCands)->mother(nGenMo) == *mCands){
              //if((*sCands)->mother(0) == *mCands){
              if(nGenMo == 0) GenInfo.mo1[geninfo_idx] = int(mCands-sel_cands.begin());
              if(nGenMo == 1) GenInfo.mo2[geninfo_idx] = int(mCands-sel_cands.begin());
            }
          }
        }
        for(int nGenDa = 0; nGenDa < std::min(4,int((*sCands)->numberOfDaughters())); nGenDa++){
          for(std::vector<const reco::Candidate *>::iterator mCands = sel_cands.begin();
              mCands != sel_cands.end(); mCands++){
            if((*sCands)->daughter(nGenDa) == *mCands){
              if(nGenDa == 0) GenInfo.da1[geninfo_idx] = int(mCands-sel_cands.begin());
              if(nGenDa == 1) GenInfo.da2[geninfo_idx] = int(mCands-sel_cands.begin());
              if(nGenDa == 2) GenInfo.da3[geninfo_idx] = int(mCands-sel_cands.begin());
              if(nGenDa == 3) GenInfo.da4[geninfo_idx] = int(mCands-sel_cands.begin());
            }
          }
        }
      }
    }//isRealData}}}
    //printf("-----*****DEBUG:End of GenInfo.\n");
    //std::cout<<"Start to fill!\n";

  }//try
  catch (std::exception & err){
    std::cout  << "Exception during event number: " << iEvent.id()
               << "\n" << err.what() << "\n";
  }//catch 
  root->Fill();
  //std::cout<<"filled!\n";
    
  //Made a Bntuple on the fly
  if(makeBntuple_){
    int ifchannel[8];
    for(int ichannel=0; ichannel<8; ichannel++){ ifchannel[ichannel] = Bchannel_[ichannel]; }
    bool REAL = iEvent.isRealData();
    int Btypesize[8]={0,0,0,0,0,0,0,0};
    Bntuple->makeNtuple(ifchannel, Btypesize, REAL, doBntupleSkim_, &EvtInfo, &VtxInfo, &MuonInfo, &TrackInfo, &BInfo, &GenInfo, nt0, nt1, nt2, nt3, nt5, nt6, nt7);
    if(!REAL) Bntuple->fillGenTree(ntGen, &GenInfo);
  }
}

// ------------ method called once each job just after ending the event loop  ------------{{{
void Bfinder::endJob()
{
}

// ------------ method called when starting to processes a run  ------------
void Bfinder::beginRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a run  ------------
void Bfinder::endRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when starting to processes a luminosity block  ------------
void Bfinder::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a luminosity block  ------------
void Bfinder::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void Bfinder::fillDescriptions(edm::ConfigurationDescriptions& descriptions)
{
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}//}}}

//BranchOut2MuTk{{{
void Bfinder::BranchOut2MuTk(
                             BInfoBranches &BInfo, 
                             // std::vector<pat::PackedCandidate> input_tracks,
                             std::vector<const reco::Track*> input_tracks,
                             reco::Vertex thePrimaryV,
                             std::vector<bool> isNeededTrack,
                             TLorentzVector v4_mu1, 
                             TLorentzVector v4_mu2,
                             reco::TransientTrack muonPTT,
                             reco::TransientTrack muonMTT,
                             std::vector<int> &B_counter,
                             float *mass_window,
                             float MuMu_MASS,
                             float Tk_MASS,
                             int channel_number

                             ){
  if(channel_number > (int)Bchannel_.size()){ printf("Exceeding defined # of channel, exit"); return;}
  float chi = 0.;
  float ndf = 0.;
  int tk1_hindex = -1;
  KinematicParticleFactoryFromTransientTrack pFactory;
  ParticleMass muon_mass = MUON_MASS; //pdg mass
  float muon_sigma = Functs.getParticleSigma(muon_mass);

  for(auto tk_it_it1=input_tracks.begin(); //
      tk_it_it1 != input_tracks.end() ; tk_it_it1++){
    tk1_hindex = int(tk_it_it1 - input_tracks.begin());
    if(tk1_hindex>=int(isNeededTrack.size())) break;
    if (!isNeededTrack[tk1_hindex]) continue;
    auto tk_it1 = (*tk_it_it1);

    if (abs(tk_it1->charge()) != 1) continue;
      
    TLorentzVector v4_tk1;
    v4_tk1.SetPtEtaPhiM(tk_it1->pt(),tk_it1->eta(),tk_it1->phi(),KAON_MASS);
  
    //if ((v4_mu1+v4_mu2+v4_tk1).Mag()<mass_window[0]-0.2 || (v4_mu1+v4_mu2+v4_tk1).Mag()>mass_window[1]+0.2) continue;
    if ((v4_mu1+v4_mu2+v4_tk1).Mag()<mass_window[0] || (v4_mu1+v4_mu2+v4_tk1).Mag()>mass_window[1]) continue;
    XbMassCutLevel[channel_number-1]->Fill(0);
    if((v4_mu1+v4_mu2+v4_tk1).Pt()<bPtCut_[channel_number-1])continue;
    XbMassCutLevel[channel_number-1]->Fill(1);
      
    reco::TransientTrack kaonTT(*tk_it1, &(*bField) );
    if (!kaonTT.isValid()) continue;
    XbMassCutLevel[channel_number-1]->Fill(2);
      
    ParticleMass kaon_mass = Tk_MASS;
    float kaon_sigma = Functs.getParticleSigma(kaon_mass);
      
    std::vector<RefCountedKinematicParticle> Xb_candidate;
    Xb_candidate.push_back(pFactory.particle(muonPTT,muon_mass,chi,ndf,muon_sigma));
    Xb_candidate.push_back(pFactory.particle(muonMTT,muon_mass,chi,ndf,muon_sigma));
    Xb_candidate.push_back(pFactory.particle(kaonTT,kaon_mass,chi,ndf,kaon_sigma));
    RefCountedKinematicTree xbVFT;

    double MaximumDoca = Functs.getMaxDoca(Xb_candidate);
    if (MaximumDoca > MaxDocaCut_[channel_number-1]) continue;
    XbMassCutLevel[channel_number-1]->Fill(3);
      
    ParticleMass uj_mass = MuMu_MASS;
    MultiTrackKinematicConstraint *uj_c = new TwoTrackMassKinematicConstraint(uj_mass);
    KinematicConstrainedVertexFitter kcvFitter;
    xbVFT = kcvFitter.fit(Xb_candidate, uj_c);
    if (!xbVFT->isValid()) continue;
    XbMassCutLevel[channel_number-1]->Fill(4);

    xbVFT->movePointerToTheTop();
    RefCountedKinematicParticle     xbVFP       = xbVFT->currentParticle();
    RefCountedKinematicVertex       xbVFPvtx    = xbVFT->currentDecayVertex();
    if (!xbVFPvtx->vertexIsValid()) continue;
    XbMassCutLevel[channel_number-1]->Fill(5);

    std::vector<RefCountedKinematicParticle> xCands  = xbVFT->finalStateParticles();
      
    double chi2_prob = TMath::Prob(xbVFPvtx->chiSquared(),xbVFPvtx->degreesOfFreedom());
    if (chi2_prob < VtxChiProbCut_[channel_number-1]) continue;
    XbMassCutLevel[channel_number-1]->Fill(6);
      
    //if (xbVFP->currentState().mass()<mass_window[0] || xbVFP->currentState().mass()>mass_window[1]) continue;
      
    TLorentzVector xb_4vec,xb_mu1_4vec,xb_mu2_4vec,xb_tk1_4vec,xb_tk2_4vec;
    xb_4vec.SetPxPyPzE(xbVFP->currentState().kinematicParameters().momentum().x(),
                       xbVFP->currentState().kinematicParameters().momentum().y(),
                       xbVFP->currentState().kinematicParameters().momentum().z(),
                       xbVFP->currentState().kinematicParameters().energy());
    xb_mu1_4vec.SetPxPyPzE(xCands[0]->currentState().kinematicParameters().momentum().x(),
                           xCands[0]->currentState().kinematicParameters().momentum().y(),
                           xCands[0]->currentState().kinematicParameters().momentum().z(),
                           xCands[0]->currentState().kinematicParameters().energy());
    xb_mu2_4vec.SetPxPyPzE(xCands[1]->currentState().kinematicParameters().momentum().x(),
                           xCands[1]->currentState().kinematicParameters().momentum().y(),
                           xCands[1]->currentState().kinematicParameters().momentum().z(),
                           xCands[1]->currentState().kinematicParameters().energy());
    xb_tk1_4vec.SetPxPyPzE(xCands[2]->currentState().kinematicParameters().momentum().x(),
                           xCands[2]->currentState().kinematicParameters().momentum().y(),
                           xCands[2]->currentState().kinematicParameters().momentum().z(),
                           xCands[2]->currentState().kinematicParameters().energy());
      
    BInfo.index[BInfo.size]   = BInfo.size;
    BInfo.unfitted_mass[BInfo.size] = (v4_mu1+v4_mu2+v4_tk1).Mag();
    BInfo.unfitted_pt[BInfo.size] = (v4_mu1+v4_mu2+v4_tk1).Pt();
    BInfo.mass[BInfo.size]    = xb_4vec.Mag();
    BInfo.pt[BInfo.size]    = xb_4vec.Pt();
    BInfo.eta[BInfo.size]    = xb_4vec.Eta();
    BInfo.phi[BInfo.size]    = xb_4vec.Phi();
    BInfo.px[BInfo.size]      = xb_4vec.Px();
    BInfo.py[BInfo.size]      = xb_4vec.Py();
    BInfo.pz[BInfo.size]      = xb_4vec.Pz();
    BInfo.pxE[BInfo.size]     = sqrt(xbVFP->currentState().kinematicParametersError().matrix()(3,3));
    BInfo.pyE[BInfo.size]     = sqrt(xbVFP->currentState().kinematicParametersError().matrix()(4,4));
    BInfo.pzE[BInfo.size]     = sqrt(xbVFP->currentState().kinematicParametersError().matrix()(5,5));
    BInfo.MaxDoca[BInfo.size]         = MaximumDoca;

    VertexDistance3D a3d;
    //https://github.com/cms-sw/cmssw/blob/CMSSW_7_5_0/RecoVertex/VertexTools/src/VertexDistance3D.cc
    BInfo.svpvDistance[BInfo.size] = a3d.distance(thePrimaryV,xbVFPvtx->vertexState()).value();
    BInfo.svpvDisErr[BInfo.size] = a3d.distance(thePrimaryV,xbVFPvtx->vertexState()).error();
    if( (BInfo.svpvDistance[BInfo.size]/BInfo.svpvDisErr[BInfo.size]) < svpvDistanceCut_[channel_number-1]) continue;
    XbMassCutLevel[channel_number-1]->Fill(7);

    reco::Vertex::Point vp1(thePrimaryV.position().x(), thePrimaryV.position().y(), 0.);
    reco::Vertex::Point vp2(xbVFPvtx->vertexState().position().x(), xbVFPvtx->vertexState().position().y(), 0.);
    ROOT::Math::SVector<double, 6> sv1(thePrimaryV.covariance(0,0), thePrimaryV.covariance(0,1), thePrimaryV.covariance(1,1), 0., 0., 0.);
    ROOT::Math::SVector<double, 6> sv2(xbVFPvtx->vertexState().error().cxx(), xbVFPvtx->vertexState().error().cyx(), xbVFPvtx->vertexState().error().cyy(), 0., 0., 0.);
    reco::Vertex::Error ve1(sv1);
    reco::Vertex::Error ve2(sv2);
    reco::Vertex v1(vp1, ve1);
    reco::Vertex v2(vp2, ve2);
    BInfo.svpvDistance_2D[BInfo.size] = a3d.distance(v1, v2).value();
    BInfo.svpvDisErr_2D[BInfo.size] = a3d.distance(v1, v2).error();

    BInfo.vtxX[BInfo.size]    = xbVFPvtx->position().x();
    BInfo.vtxY[BInfo.size]    = xbVFPvtx->position().y();
    BInfo.vtxZ[BInfo.size]    = xbVFPvtx->position().z();
    BInfo.vtxXErr[BInfo.size]   = xbVFPvtx->error().cxx();
    BInfo.vtxYErr[BInfo.size]   = xbVFPvtx->error().cyy();
    BInfo.vtxZErr[BInfo.size]   = xbVFPvtx->error().czz();
    BInfo.vtxYXErr[BInfo.size]  = xbVFPvtx->error().cyx();
    BInfo.vtxZXErr[BInfo.size]  = xbVFPvtx->error().czx();
    BInfo.vtxZYErr[BInfo.size]  = xbVFPvtx->error().czy();
    BInfo.vtxdof[BInfo.size]  = xbVFPvtx->degreesOfFreedom();
    BInfo.vtxchi2[BInfo.size] = xbVFPvtx->chiSquared();

    TVector3 svpvVec;
    svpvVec.SetXYZ(BInfo.vtxX[BInfo.size]-EvtInfo.PVx, BInfo.vtxY[BInfo.size]-EvtInfo.PVy, BInfo.vtxZ[BInfo.size]-EvtInfo.PVz);
    TVector3 dVec;
    dVec.SetXYZ(BInfo.px[BInfo.size], BInfo.py[BInfo.size], BInfo.pz[BInfo.size]);
    BInfo.alpha[BInfo.size] = svpvVec.Angle(dVec);
    if( BInfo.alpha[BInfo.size] > alphaCut_[channel_number-1]) continue;
    XbMassCutLevel[channel_number-1]->Fill(8);
      
    BInfo.rftk1_index[BInfo.size] = -2;
    BInfo.rftk2_index[BInfo.size] = -2;
    BInfo.rfuj_index[BInfo.size]  = BInfo.uj_size-1;
    BInfo.rftk1_index[BInfo.size] = tk1_hindex;
    BInfo.rftk2_index[BInfo.size] = tk1_hindex;
      
    BInfo.rfmu1_pt[BInfo.size] =xb_mu1_4vec.Pt();
    BInfo.rfmu1_eta[BInfo.size]=xb_mu1_4vec.Eta();
    BInfo.rfmu1_phi[BInfo.size]=xb_mu1_4vec.Phi();
    BInfo.rfmu2_pt[BInfo.size] =xb_mu2_4vec.Pt();
    BInfo.rfmu2_eta[BInfo.size]=xb_mu2_4vec.Eta();
    BInfo.rfmu2_phi[BInfo.size]=xb_mu2_4vec.Phi();
    BInfo.rftk1_pt[BInfo.size] =xb_tk1_4vec.Pt();
    BInfo.rftk1_eta[BInfo.size]=xb_tk1_4vec.Eta();
    BInfo.rftk1_phi[BInfo.size]=xb_tk1_4vec.Phi();
    BInfo.rftk2_pt[BInfo.size] =-999.;
    BInfo.rftk2_eta[BInfo.size]=-999.;
    BInfo.rftk2_phi[BInfo.size]=-999.;
      
    BInfo.type[BInfo.size] = channel_number;
    B_counter[channel_number-1]++;
      
    Xb_candidate.clear();
    xCands.clear();
    BInfo.size++;
  }//Tk1
}
//}}}

//BranchOut2MuX{{{
void Bfinder::BranchOut2MuX_XtoTkTk(
                                    BInfoBranches &BInfo, 
                                    // std::vector<pat::PackedCandidate> input_tracks,
                                    std::vector<const reco::Track*> input_tracks,
                                    reco::Vertex thePrimaryV,
                                    std::vector<bool> isNeededTrack,
                                    TLorentzVector v4_mu1, 
                                    TLorentzVector v4_mu2,
                                    reco::TransientTrack muonPTT,
                                    reco::TransientTrack muonMTT,
                                    std::vector<int> &B_counter,
                                    float *mass_window,
                                    float MuMu_MASS,
                                    float TkTk_MASS,
                                    float TkTk_window,
                                    float Tk1_MASS,
                                    float Tk2_MASS,
                                    int channel_number,
                                    int fit_option
                                    ){
  if(channel_number > (int)Bchannel_.size()){ printf("Exceeding defined # of channel, exit"); return;}
  float chi = 0.;
  float ndf = 0.;
  KinematicParticleFactoryFromTransientTrack pFactory;
  ParticleMass muon_mass = MUON_MASS; //pdg mass
  float muon_sigma = Functs.getParticleSigma(muon_mass);
  int tk1_hindex = -1;
  int tk2_hindex = -1;

  for(auto tk_it_it1=input_tracks.begin(); // 
      tk_it_it1 != input_tracks.end() ; tk_it_it1++){
    tk1_hindex = int(tk_it_it1 - input_tracks.begin());
    if(tk1_hindex>=int(isNeededTrack.size())) break;
    if (!isNeededTrack[tk1_hindex]) continue;
    auto tk_it1 = (*tk_it_it1);
  // for(int tk1idx = 0; tk1idx < (int)isNeededTrackIdx.size(); tk1idx++){
  //   tk1_hindex = isNeededTrackIdx[tk1idx];  
  //   auto tk_it1 = input_tracks[tk1_hindex];

    if (tk_it1->charge()<0) continue;
    
    for(auto tk_it_it2=input_tracks.begin();
        tk_it_it2 != input_tracks.end() ; tk_it_it2++){
      if (BInfo.size >= MAX_XB) break;
      tk2_hindex = int(tk_it_it2 - input_tracks.begin());
      if(tk2_hindex>=int(isNeededTrack.size())) break;
      if (!isNeededTrack[tk2_hindex]) continue;
      auto tk_it2 = (*tk_it_it2);

    // for(int tk2idx = 0; tk2idx < (int)isNeededTrackIdx.size(); tk2idx++){
    //   tk2_hindex = isNeededTrackIdx[tk2idx];
    //   auto tk_it2 = input_tracks[tk2_hindex];
      
      if (tk_it2->charge()>0) continue;
            
      TLorentzVector v4_tk1,v4_tk2;
      v4_tk1.SetPtEtaPhiM(tk_it1->pt(),tk_it1->eta(),tk_it1->phi(),Tk1_MASS);
      v4_tk2.SetPtEtaPhiM(tk_it2->pt(),tk_it2->eta(),tk_it2->phi(),Tk2_MASS);
      if(TkTk_MASS > 0) {if (fabs((v4_tk1+v4_tk2).Mag()-TkTk_MASS)>TkTk_window) continue;}
      //else {if (fabs((v4_tk1+v4_tk2).Mag())>TkTk_window) continue;}//if no tktk mass constrain, require it to be at least < some mass value
      XbMassCutLevel[channel_number-1]->Fill(0);
            
      //if ((v4_mu1+v4_mu2+v4_tk1+v4_tk2).Mag()<mass_window[0]-0.2 || (v4_mu1+v4_mu2+v4_tk1+v4_tk2).Mag()>mass_window[1]+0.2) continue;
      if ((v4_mu1+v4_mu2+v4_tk1+v4_tk2).Mag()<mass_window[0] || (v4_mu1+v4_mu2+v4_tk1+v4_tk2).Mag()>mass_window[1]) continue;
      XbMassCutLevel[channel_number-1]->Fill(1);
      if((v4_mu1+v4_mu2+v4_tk1+v4_tk2).Pt()<bPtCut_[channel_number-1])continue;
      XbMassCutLevel[channel_number-1]->Fill(2);
            
      reco::TransientTrack tk1PTT(*tk_it1, &(*bField) );
      reco::TransientTrack tk2MTT(*tk_it2, &(*bField) );
      // reco::TransientTrack tk1PTT(tk_it1->pseudoTrack(), &(*bField) );
      // reco::TransientTrack tk2MTT(tk_it2->pseudoTrack(), &(*bField) );
      if (!tk1PTT.isValid()) continue;
      if (!tk2MTT.isValid()) continue;
      XbMassCutLevel[channel_number-1]->Fill(3);
            
      ParticleMass tk1_mass = Tk1_MASS;
      float tk1_sigma = Functs.getParticleSigma(tk1_mass);
      ParticleMass tk2_mass = Tk2_MASS;
      float tk2_sigma = Functs.getParticleSigma(tk2_mass);

      //doing tktk fit
      std::vector<RefCountedKinematicParticle> tktk_candidate;
      KinematicParticleVertexFitter   tktk_fitter;
      RefCountedKinematicTree         tktk_VFT;
      RefCountedKinematicParticle tktk_VFP;
      RefCountedKinematicVertex   tktk_VFPvtx;

      tktk_candidate.push_back(pFactory.particle(tk1PTT,tk1_mass,chi,ndf,tk1_sigma));
      tktk_candidate.push_back(pFactory.particle(tk2MTT,tk2_mass,chi,ndf,tk2_sigma));
      tktk_VFT = tktk_fitter.fit(tktk_candidate);
      if (tktk_VFT->isValid()){
        XbMassCutLevel[channel_number-1]->Fill(4);
        tktk_VFT->movePointerToTheTop();
        tktk_VFP   = tktk_VFT->currentParticle();
        tktk_VFPvtx = tktk_VFT->currentDecayVertex();
        if (tktk_VFPvtx->vertexIsValid()){
          XbMassCutLevel[channel_number-1]->Fill(5);
          double chi2_prob_tktk = TMath::Prob(tktk_VFPvtx->chiSquared(),tktk_VFPvtx->degreesOfFreedom());
          if (chi2_prob_tktk >= VtxChiProbCut_[channel_number-1]){
            XbMassCutLevel[channel_number-1]->Fill(6);
          }
          else if (TkTk_MASS > 0) continue;
        }
        else if (TkTk_MASS > 0) continue;
      }
      else if (TkTk_MASS > 0) continue;

      std::vector<RefCountedKinematicParticle> Xb_candidate;
      Xb_candidate.push_back(pFactory.particle(muonPTT,muon_mass,chi,ndf,muon_sigma));
      Xb_candidate.push_back(pFactory.particle(muonMTT,muon_mass,chi,ndf,muon_sigma));
      if(fit_option == 0 || TkTk_MASS < 0){
        Xb_candidate.push_back(pFactory.particle(tk1PTT,tk1_mass,chi,ndf,tk1_sigma));
        Xb_candidate.push_back(pFactory.particle(tk2MTT,tk2_mass,chi,ndf,tk2_sigma));
      }
      else if(fit_option == 1){
        VirtualKinematicParticleFactory vFactory;
        float tktkchi = tktk_VFPvtx->chiSquared();
        float tktkndf = tktk_VFPvtx->degreesOfFreedom();
        Xb_candidate.push_back(vFactory.particle(tktk_VFP->currentState(),tktkchi,tktkndf,tktk_VFP));
      }
      else{
        std::cout<<"Unknown fit option, exit"<<std::endl;
        return;
      }

      KinematicConstrainedVertexFitter kcvFitter;
      RefCountedKinematicTree xbVFT;
      RefCountedKinematicParticle     xbVFP; 
      RefCountedKinematicVertex       xbVFPvtx;

      double MaximumDoca = Functs.getMaxDoca(Xb_candidate);
      if (MaximumDoca > MaxDocaCut_[channel_number-1]) continue;
      XbMassCutLevel[channel_number-1]->Fill(7);
            
      ParticleMass uj_mass = MuMu_MASS;
      MultiTrackKinematicConstraint *uj_c = new  TwoTrackMassKinematicConstraint(uj_mass);
      xbVFT = kcvFitter.fit(Xb_candidate, uj_c);
      if (!xbVFT->isValid()) continue;
      XbMassCutLevel[channel_number-1]->Fill(8);

      xbVFT->movePointerToTheTop();
      xbVFP       = xbVFT->currentParticle();
      xbVFPvtx    = xbVFT->currentDecayVertex();
      if (!xbVFPvtx->vertexIsValid()) continue;
      XbMassCutLevel[channel_number-1]->Fill(9);
            
      double chi2_prob = TMath::Prob(xbVFPvtx->chiSquared(),xbVFPvtx->degreesOfFreedom());
      if (chi2_prob < VtxChiProbCut_[channel_number-1]) continue;
      XbMassCutLevel[channel_number-1]->Fill(10);
            
      //Cut out a mass window
      //if (xbVFP->currentState().mass()<mass_window[0]|| xbVFP->currentState().mass()>mass_window[1]) continue;
      //
      std::vector<RefCountedKinematicParticle> xCands  = xbVFT->finalStateParticles();
            
      TLorentzVector xb_4vec,xb_mu1_4vec,xb_mu2_4vec,tktk_4vec,xb_tk1_4vec,xb_tk2_4vec,tktk_tk1_4vec, tktk_tk2_4vec;
      xb_4vec.SetPxPyPzE(xbVFP->currentState().kinematicParameters().momentum().x(),
                         xbVFP->currentState().kinematicParameters().momentum().y(),
                         xbVFP->currentState().kinematicParameters().momentum().z(),
                         xbVFP->currentState().kinematicParameters().energy());
      xb_mu1_4vec.SetPxPyPzE(xCands[0]->currentState().kinematicParameters().momentum().x(),
                             xCands[0]->currentState().kinematicParameters().momentum().y(),
                             xCands[0]->currentState().kinematicParameters().momentum().z(),
                             xCands[0]->currentState().kinematicParameters().energy());
      xb_mu2_4vec.SetPxPyPzE(xCands[1]->currentState().kinematicParameters().momentum().x(),
                             xCands[1]->currentState().kinematicParameters().momentum().y(),
                             xCands[1]->currentState().kinematicParameters().momentum().z(),
                             xCands[1]->currentState().kinematicParameters().energy());

      if(fit_option == 0){
        xb_tk1_4vec.SetPxPyPzE(xCands[2]->currentState().kinematicParameters().momentum().x(),
                               xCands[2]->currentState().kinematicParameters().momentum().y(),
                               xCands[2]->currentState().kinematicParameters().momentum().z(),
                               xCands[2]->currentState().kinematicParameters().energy());
        xb_tk2_4vec.SetPxPyPzE(xCands[3]->currentState().kinematicParameters().momentum().x(),
                               xCands[3]->currentState().kinematicParameters().momentum().y(),
                               xCands[3]->currentState().kinematicParameters().momentum().z(),
                               xCands[3]->currentState().kinematicParameters().energy());
      }
            
      if(fit_option == 1){
        xb_tk1_4vec.SetPxPyPzE(xCands[2]->currentState().kinematicParameters().momentum().x(),
                               xCands[2]->currentState().kinematicParameters().momentum().y(),
                               xCands[2]->currentState().kinematicParameters().momentum().z(),
                               xCands[2]->currentState().kinematicParameters().energy());
      }
            
      BInfo.index[BInfo.size]   = BInfo.size;
      BInfo.mass[BInfo.size]    = xb_4vec.Mag();
      BInfo.unfitted_mass[BInfo.size] = (v4_mu1+v4_mu2+v4_tk1+v4_tk2).Mag();
      BInfo.unfitted_pt[BInfo.size] = (v4_mu1+v4_mu2+v4_tk1+v4_tk2).Pt();
      BInfo.pt[BInfo.size]    = xb_4vec.Pt();
      BInfo.eta[BInfo.size]    = xb_4vec.Eta();
      BInfo.phi[BInfo.size]    = xb_4vec.Phi();
      BInfo.px[BInfo.size]      = xb_4vec.Px();
      BInfo.py[BInfo.size]      = xb_4vec.Py();
      BInfo.pz[BInfo.size]      = xb_4vec.Pz();
      BInfo.pxE[BInfo.size]     = sqrt(xbVFP->currentState().kinematicParametersError().matrix()(3,3));
      BInfo.pyE[BInfo.size]     = sqrt(xbVFP->currentState().kinematicParametersError().matrix()(4,4));
      BInfo.pzE[BInfo.size]     = sqrt(xbVFP->currentState().kinematicParametersError().matrix()(5,5));
      BInfo.MaxDoca[BInfo.size]         = MaximumDoca;

      VertexDistance3D a3d;
      //https://github.com/cms-sw/cmssw/blob/CMSSW_7_5_0/RecoVertex/VertexTools/src/VertexDistance3D.cc
      BInfo.svpvDistance[BInfo.size] = a3d.distance(thePrimaryV,xbVFPvtx->vertexState()).value();
      BInfo.svpvDisErr[BInfo.size] = a3d.distance(thePrimaryV,xbVFPvtx->vertexState()).error();
      if( (BInfo.svpvDistance[BInfo.size]/BInfo.svpvDisErr[BInfo.size]) < svpvDistanceCut_[channel_number-1]) continue;
      XbMassCutLevel[channel_number-1]->Fill(11);
           
      reco::Vertex::Point vp1(thePrimaryV.position().x(), thePrimaryV.position().y(), 0.);
      reco::Vertex::Point vp2(xbVFPvtx->vertexState().position().x(), xbVFPvtx->vertexState().position().y(), 0.);
      ROOT::Math::SVector<double, 6> sv1(thePrimaryV.covariance(0,0), thePrimaryV.covariance(0,1), thePrimaryV.covariance(1,1), 0., 0., 0.);
      ROOT::Math::SVector<double, 6> sv2(xbVFPvtx->vertexState().error().cxx(), xbVFPvtx->vertexState().error().cyx(), xbVFPvtx->vertexState().error().cyy(), 0., 0., 0.);
      reco::Vertex::Error ve1(sv1);
      reco::Vertex::Error ve2(sv2);
      reco::Vertex v1(vp1, ve1);
      reco::Vertex v2(vp2, ve2);
      BInfo.svpvDistance_2D[BInfo.size] = a3d.distance(v1, v2).value();
      BInfo.svpvDisErr_2D[BInfo.size] = a3d.distance(v1, v2).error();

      BInfo.vtxX[BInfo.size]    = xbVFPvtx->position().x();
      BInfo.vtxY[BInfo.size]    = xbVFPvtx->position().y();
      BInfo.vtxZ[BInfo.size]    = xbVFPvtx->position().z();
      BInfo.vtxXErr[BInfo.size]        = xbVFPvtx->error().cxx();
      BInfo.vtxYErr[BInfo.size]        = xbVFPvtx->error().cyy();
      BInfo.vtxZErr[BInfo.size]        = xbVFPvtx->error().czz();
      BInfo.vtxYXErr[BInfo.size]     = xbVFPvtx->error().cyx();
      BInfo.vtxZXErr[BInfo.size]     = xbVFPvtx->error().czx();
      BInfo.vtxZYErr[BInfo.size]     = xbVFPvtx->error().czy();
      BInfo.vtxdof[BInfo.size]  = xbVFPvtx->degreesOfFreedom();
      BInfo.vtxchi2[BInfo.size] = xbVFPvtx->chiSquared();

      TVector3 svpvVec;
      svpvVec.SetXYZ(BInfo.vtxX[BInfo.size]-EvtInfo.PVx, BInfo.vtxY[BInfo.size]-EvtInfo.PVy, BInfo.vtxZ[BInfo.size]-EvtInfo.PVz);
      TVector3 dVec;
      dVec.SetXYZ(BInfo.px[BInfo.size], BInfo.py[BInfo.size], BInfo.pz[BInfo.size]);
      BInfo.alpha[BInfo.size] = svpvVec.Angle(dVec);
      if( BInfo.alpha[BInfo.size] > alphaCut_[channel_number-1]) continue;
      XbMassCutLevel[channel_number-1]->Fill(12);
            
      BInfo.rftk1_index[BInfo.size] = -2;
      BInfo.rftk2_index[BInfo.size] = -2;
      BInfo.rfuj_index[BInfo.size]  = BInfo.uj_size-1;
      BInfo.rftk1_index[BInfo.size] = tk1_hindex;
      BInfo.rftk2_index[BInfo.size] = tk2_hindex;
            
      //tktk fit info
      BInfo.tktk_unfitted_mass[BInfo.size]    = (v4_tk1+v4_tk2).Mag();
      BInfo.tktk_unfitted_pt[BInfo.size]    = (v4_tk1+v4_tk2).Pt();
      if(tktk_VFT->isValid() && tktk_VFPvtx->vertexIsValid()){
        std::vector<RefCountedKinematicParticle> tktkCands  = tktk_VFT->finalStateParticles();
        tktk_4vec.SetPxPyPzE(tktk_VFP->currentState().kinematicParameters().momentum().x(),
                             tktk_VFP->currentState().kinematicParameters().momentum().y(),
                             tktk_VFP->currentState().kinematicParameters().momentum().z(),
                             tktk_VFP->currentState().kinematicParameters().energy());            
        tktk_tk1_4vec.SetPxPyPzE(tktkCands[0]->currentState().kinematicParameters().momentum().x(),
                                 tktkCands[0]->currentState().kinematicParameters().momentum().y(),
                                 tktkCands[0]->currentState().kinematicParameters().momentum().z(),
                                 tktkCands[0]->currentState().kinematicParameters().energy());
        tktk_tk2_4vec.SetPxPyPzE(tktkCands[1]->currentState().kinematicParameters().momentum().x(),
                                 tktkCands[1]->currentState().kinematicParameters().momentum().y(),
                                 tktkCands[1]->currentState().kinematicParameters().momentum().z(),
                                 tktkCands[1]->currentState().kinematicParameters().energy());
        BInfo.tktk_mass[BInfo.size]    = tktk_4vec.Mag();
        BInfo.tktk_pt[BInfo.size]      = tktk_4vec.Pt();
        BInfo.tktk_eta[BInfo.size]     = tktk_4vec.Eta();
        BInfo.tktk_phi[BInfo.size]     = tktk_4vec.Phi();
        BInfo.tktk_px[BInfo.size]      = tktk_4vec.Px();
        BInfo.tktk_py[BInfo.size]      = tktk_4vec.Py();
        BInfo.tktk_pz[BInfo.size]      = tktk_4vec.Pz();
        BInfo.tktk_vtxX[BInfo.size]    = tktk_VFPvtx->position().x();
        BInfo.tktk_vtxY[BInfo.size]    = tktk_VFPvtx->position().y();
        BInfo.tktk_vtxZ[BInfo.size]    = tktk_VFPvtx->position().z();
        BInfo.tktk_vtxXErr[BInfo.size] = tktk_VFPvtx->error().cxx();
        BInfo.tktk_vtxYErr[BInfo.size] = tktk_VFPvtx->error().cyy();
        BInfo.tktk_vtxZErr[BInfo.size] = tktk_VFPvtx->error().czz();
        BInfo.tktk_vtxYXErr[BInfo.size]= tktk_VFPvtx->error().cyx();
        BInfo.tktk_vtxZXErr[BInfo.size]= tktk_VFPvtx->error().czx();
        BInfo.tktk_vtxZYErr[BInfo.size]= tktk_VFPvtx->error().czy();
        BInfo.tktk_vtxdof[BInfo.size]  = tktk_VFPvtx->degreesOfFreedom();
        BInfo.tktk_vtxchi2[BInfo.size] = tktk_VFPvtx->chiSquared();
        BInfo.tktk_rftk1_pt[BInfo.size] =tktk_tk1_4vec.Pt();
        BInfo.tktk_rftk1_eta[BInfo.size]=tktk_tk1_4vec.Eta();
        BInfo.tktk_rftk1_phi[BInfo.size]=tktk_tk1_4vec.Phi();
        BInfo.tktk_rftk2_pt[BInfo.size] =tktk_tk2_4vec.Pt();
        BInfo.tktk_rftk2_eta[BInfo.size]=tktk_tk2_4vec.Eta();
        BInfo.tktk_rftk2_phi[BInfo.size]=tktk_tk2_4vec.Phi();
      }
      else{
        BInfo.tktk_mass[BInfo.size]    = -1;
      }
            
      BInfo.rfmu1_pt[BInfo.size] =xb_mu1_4vec.Pt();
      BInfo.rfmu1_eta[BInfo.size]=xb_mu1_4vec.Eta();
      BInfo.rfmu1_phi[BInfo.size]=xb_mu1_4vec.Phi();
      BInfo.rfmu2_pt[BInfo.size] =xb_mu2_4vec.Pt();
      BInfo.rfmu2_eta[BInfo.size]=xb_mu2_4vec.Eta();
      BInfo.rfmu2_phi[BInfo.size]=xb_mu2_4vec.Phi();
      //If option == 1, this momentum is the tktk virtual particle p.
      BInfo.rftk1_pt[BInfo.size] =xb_tk1_4vec.Pt();
      BInfo.rftk1_eta[BInfo.size]=xb_tk1_4vec.Eta();
      BInfo.rftk1_phi[BInfo.size]=xb_tk1_4vec.Phi();
      if(fit_option == 0){
        BInfo.rftk2_pt[BInfo.size] =xb_tk2_4vec.Pt();
        BInfo.rftk2_eta[BInfo.size]=xb_tk2_4vec.Eta();
        BInfo.rftk2_phi[BInfo.size]=xb_tk2_4vec.Phi();
      }
      else if(fit_option == 1){
        BInfo.rftk2_pt[BInfo.size]=-999;
        BInfo.rftk2_eta[BInfo.size]=-999;
        BInfo.rftk2_phi[BInfo.size]=-999;
      }
            
      BInfo.type[BInfo.size] = channel_number;
      B_counter[channel_number-1]++;
            
      Xb_candidate.clear();
      xCands.clear();
      BInfo.size++;
    }//Tk2
  }//Tk1
}
//}}}

//define this as a plug-in
DEFINE_FWK_MODULE(Bfinder);
