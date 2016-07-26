// Include files
// from bost
#include <boost/assign/list_of.hpp>
#include <boost/array.hpp>
#include <boost/foreach.hpp>
// from Gaudi
#include "GaudiKernel/AlgFactory.h"
#include "Event/Track.h"
#include "Event/StateParameters.h"
#include "Math/CholeskyDecomp.h"
// local
#include "PrLineFitterY.h"
#include "PrHybridSeeding.h"
#include "PrPlaneHybridCounter.h"
#include "Event/FTCluster.h"
#include "Event/FTLiteCluster.h"

#ifdef MVA_used
#include "Event/MCTrackInfo.h"
#include "Event/MCProperty.h"
#include "Event/MCHit.h"
#include "Event/MCParticle.h"
#include "Linker/LinkedTo.h"
#include "Linker/LinkedFrom.h"
#include "Linker/AllLinks.h"
#include "Event/LinksByKey.h"
#include "GaudiAlg/GaudiTupleAlg.h"
#else
#include "GaudiAlg/GaudiAlgorithm.h"
#endif

//-----------------------------------------------------------------------------
// Implementation file for class : PrHybridSeeding
// 
// @author Renato Quagliani (rquaglia@cern.ch)
// @date 2015-03-11 
//-----------------------------------------------------------------------------

// Declaration of the Algorithm Factory
DECLARE_ALGORITHM_FACTORY( PrHybridSeeding )

//=============================================================================
// Standard constructor, initializes variable//=============================================================================

PrHybridSeeding::PrHybridSeeding( const std::string& name,
                                  ISvcLocator* pSvcLocator)
:
#ifdef MVA_used
GaudiTupleAlg(name,pSvcLocator),
#else
  GaudiAlgorithm(name,pSvcLocator),
#endif
  m_hitManager(nullptr)
  ,m_debugTool(nullptr)
  ,m_timerTool(nullptr){
  
  //====================================================== Final Track State parameters ================================================
  declareProperty( "ZOutputs", m_zOutputs ={StateParameters::ZBegT , StateParameters::ZMidT, StateParameters::ZEndT});
  declareProperty( "StateErrorX2",    m_stateErrorX2    =   4. );
  declareProperty( "StateErrorY2",    m_stateErrorY2    = 400. );
  declareProperty( "StateErrorTX2",   m_stateErrorTX2   = 6e-5 );
  declareProperty( "StateErrorTY2",   m_stateErrorTY2   = 1e-4 );
  declareProperty( "MomentumToolName",          m_momentumToolName  = "FastMomentumEstimate" );
  
  declareProperty( "InputName"          ,     m_inputName = ""); // Standalone by Default !
    
  //================================================ Recover Track routine parameters ==================================================
  declareProperty( "Recover", m_recover = true);
  declareProperty( "nUsedThreshold", m_nusedthreshold = {2,1,1 });
  //Kill the tracks with [6,5,4] hits in recover container having >nUsedThreshold hits
  declareProperty( "Recover_tolTy", m_recoTolTy = 0.015);
  declareProperty( "Recover_LineFitHigh" , m_recoLineHigh = 90.0);
  declareProperty( "Recover_LineFitLow" , m_recoLineLow = 16.0);
  declareProperty( "Recover_MaxChiDoF" , m_recoFinalChi2 = 10.0);
  declareProperty( "Recover_MaxNCluster", m_recoNCluster = 3);
  declareProperty( "Recover_maxY0", m_recomaxY0 = 3000.);
  declareProperty( "Recover_minTotHits", m_recoMinTotHits = 9);
  //Recovering routine minUV hits for [6,5,4] hits on x/z projection
  declareProperty( "Recover_minUV", m_recover_minUV = {4,5,5});
  declareProperty( "Recover_OutliersChi" , m_recoChiOutlier = 4.5);
  //============================================= Global parameters ===========================================================
  declareProperty( "OutputName"         ,     m_outputName=LHCb::TrackLocation::Seed);
  declareProperty( "HitManagerName"     ,     m_hitManagerName= "PrFTHitManager");
  
  declareProperty( "XOnly"              ,     m_xOnly= false);
  declareProperty( "MinXPlanes"         ,     m_minXPlanes= 4);         
  //XZ projection finding ( start with 6->5->4 ( outlier removal). 4 hits on track are killed by default at the first fit
  //as well as the outlier removed form 5->4.
  declareProperty( "MaxNHits"           ,     m_maxNHits =12);         
  //Force algorithm to find 12 hits track
  declareProperty( "NCases"             ,     m_nCases = 3);           
  //N Cases to run ( set to 2 to have best timing ( penalising lower P tracks )
  declareProperty( "TimingMeasurement"  ,     m_doTiming= false);       
  //Measure timing of the algorithm splitting in cases up/down
  declareProperty( "PrintSettings"      ,     m_printSettings= false);  
  //Print the settings of the algorithm?
  declareProperty( "RemoveClones"       ,     m_removeClones = true);    
  declareProperty( "FracCommon"         ,     m_fracCommon=0.70);        
  //Global Clone removal % shared hits ( to be implemented )
  declareProperty( "RemoveClonesX"      ,     m_removeClonesX = true);   
  //To be optimised ( BloomFilter for instance & use of fraction , update track quality)
  declareProperty( "FlagHits"           ,     m_FlagHits = true);        //To be impoved
  declareProperty( "RemoveFlagged"      ,     m_removeFlagged = true);   //To be improved
  
  //================================Internal track model parameterization parameters===============================================
  
  declareProperty( "dRatio"             , m_dRatio = -0.000262); 
  // Const dRatio value ( By = B0+B1*z ; B1/B0 = dRatio ) for x/z projection fit
  declareProperty( "UseCorrPosition"    ,     m_useCorrPos = true);  
  // Apply dRatioPar correction to dRatio using dRatioPar
  declareProperty( "dRatioPar"          , m_dRatioPar = { 0.000267957, -8.651e-06,  4.60324e-05 } );
  // Variable dRatio value dRatio = -( dRatioPar[0] + dRatioPar[1]*Radius + dRatio[2]*Radius*Radius )
  // where Radius = sqrt( [0.0005*x(zRef)]^{2} + [0.001*y(zRef)]^{2} )
  /*
    P<5 GeV fit returns
    p0                        =  0.000246549   +/-   9.03218e-06
    p1                        =  3.11309e-05   +/-   1.90404e-05
    p2                        =  3.50241e-05   +/-   9.17754e-06
    P>5GeV fit returns
    p0                        =  0.000267957   +/-   3.80583e-06
    p1                        =   -8.651e-06   +/-   1.40445e-05
    p2                        =  4.60324e-05   +/-   9.62633e-06
  */
  declareProperty( "CConst"             , m_ConstC =  2.458e8);  // Const value to compute the backward projection
  //X_{BackProj} = ax - bx * zRef + cx*CConst
  
  //================================Hit Flagging parameters===============================================
  
  declareProperty( "SizeToFlag"            ,m_SizeFlag = {12,11,10});
  //----    If nHits < 12 Flag only tracks having Chi2DoF<Flag_MaxChi2DoF_11Hits[Case]
  declareProperty( "Flag_MaxChi2DoF_11Hits",m_MaxChi2Flag = {0.5,1.0,1.0});
  //----    If nHits < 12 Flag only tracks having |X0| ( Back. Projection ) < Flag_MaxX0_11Hits[Case]
  declareProperty( "Flag_MaxX0_11Hits"     ,m_MaxX0Flag = {100.,8000.,200.});
  
  //================================Find x/z projections parameters=======================================
  
  //-------   First / Last Layer search windows (change by case to recover hit inefficiencies)
  //--- Case0 : 1st = T1-1X Last = T3-2X
  //--- Case1 : 1st = T1-2X Last = T3-1X
  //--- Case2 : 1st = T1-1X Last = T3-1X
  //--- Hits in T3 are collected accrding to txInf = x1st/Z1st;
  declareProperty( "L0_AlphaCorr"     , m_alphaCorrection = {120.64,510.64,730.64});
  declareProperty( "L0_tolHp"         , m_TolFirstLast = {280.0,540.0,1080.0});
  
  //-------   Find x-Projection 3 - hit Combo T1/T2/T3
  //--- 2 Hit combination defines the value of the backward value at 0 ;
  //--- tx picked = (XT3 -XT1)/(ZT3-ZT1)
  //--- x0 = XT1-txPicked*ZT1
  //--- Hits are collected in T2 based on the following parameters.
  declareProperty( "maxParabolaSeedHits", m_maxParabolaSeedHits = 15); // Max Nb Of hits to process in T2-1X and T2-2X given a 2hit combination
  declareProperty( "x0Corr"             , m_x0Corr = {0.002152 , 0.001534,0.001834});  // Rotation angle for tolerance
  declareProperty( "X0SlopeChange"      , m_x0SlopeChange = {400.,500.,500.});
  declareProperty( "ToleranceX0Up"      , m_TolX0SameSign = {0.75,0.75,0.75});
  //--- x0Cut needs to be different from X0SlopeChangeDown & X0SlopeChange
  declareProperty( "x0Cut"              , m_x0Cut = {1500.,4000.,6000.});
  declareProperty( "TolAtX0Cut"         , m_tolAtX0Cut = {4.5,8.0,14.0});
  declareProperty( "X0SlopeChangeDown"  , m_x0SlopeChange2 = {2000.,2000.,2000.});
  //--- Tolerance inferior for |x0| > m_x0SlopeChange2 when x0 = m_x0Cut
  declareProperty( "TolAtX0CutOpp"      , m_tolAtx0CutOppSig = {0.75,2.0,7.0});
  declareProperty( "ToleranceX0Down"    , m_tolX0Oppsig = {0.75,0.75,0.75});

  //-------   Complete the track with remaining layers and fit X/Z projection
  //--- Tolerance for remaining layers ( after the 3 hit combination ) in mm
  declareProperty( "TolXRemaining"      , m_tolRemaining = {1.0,1.0,1.0});
  declareProperty( "maxChi2HitsX"       , m_maxChi2HitsX = {5.5,5.5,5.5});
  declareProperty( "maxChi2DoFX"        , m_maxChi2DoFX = {4.0,5.0,6.0});
    
  //================================Stereo hit search parameters======================================
  declareProperty("maxNbestCluster", m_maxNClusters = {2,4,4});
  // --------- COLLECT UV
  declareProperty("RemoveHole"            , m_removeHole     = true);  
  // Do not collect hits corresponding to tracks passing trough a hole of radius RadiusHole
  declareProperty("RadiusHole"            , m_radiusHole     = 87.0);
  declareProperty("yMin"                  , m_yMin           = -1.0  *Gaudi::Units::mm); 
  // Y min for hough clusters ( -1, 2700 ) by symmetru ( -2700, +1)
  declareProperty("yMin_TrFix"            , m_yMin_TrFix     = -2.0  *Gaudi::Units::mm); 
  // Y min for hough clusters when triangle fixing On
  // Look up tracks and load low modules (y<0) pick only hits at y>-2 mm.
  declareProperty("yMax"                  , m_yMax           = + 2700. *Gaudi::Units::mm);
  declareProperty("yMax_TrFix"            , m_yMax_TrFix     = +30.0 *Gaudi::Units::mm); 
  // Y max for hough clusters when triangle fixing On
  // Look up tracks and load low modules (y<0) pick only hits at y< 30 mm.
  declareProperty("DoAsymm"               , m_doAsymmUV      = true); 
  // Do asymmetric => take into account the stereo angle of the layers
  declareProperty( "TriangleFix"          , m_useFix         = true); 
  // Use triangle fixing
  declareProperty( "TriangleFix2ndOrder"  , m_useFix2ndOrder = true); 
  // Use triangle fixing accessing the PrHit::yMin and PrHit::yMax information
  //--------- SelectHoughCluster

  declareProperty("minUV6"  ,   m_minUV6   =   {4,4,4});
  declareProperty("minUV5"  ,   m_minUV5   =   {5,5,4});
  declareProperty("minUV4"  ,   m_minUV4   =   {6,6,5});
  
  
  //--- if NXZ + NUV (Planes) <=10  XZChi2DoF + YLineChi2DoF < Chi2LowLine
  //--- 10 or 9 hits on track (= #layers by construction)
  declareProperty("Chi2LowLine"     ,    m_Chi2LowLine ={5.0,6.5,7.5} );
  //--- 11 or 12 hits on track
  declareProperty("Chi2HighLine"    ,    m_Chi2HighLine = {30.0, 50.0, 80.0});
  //--- Minimal total number of hits on track
  declareProperty("minTot"          ,    m_minTot = {9,9,9});
  //--- ty tolerance for stereo hit search
  declareProperty( "TolTyOffset"    ,    m_tolTyOffset = {0.0017, 0.0025, 0.0035});
  declareProperty( "TolTySlope"     ,    m_tolTySlope =  {0.0, 0.025 , 0.035});
  //-------- Simultaneous Fit Y
  //--- MaxChi2 Hit >10 Hits
  //--- N Layers > 10 ( 11 and 12 hits ) outliers removed if MaxChi2Hit > maxChi2Hits_11and12Hit
  declareProperty( "maxChi2Hits_11and12Hit" , m_maxChi2HitFullFitHigh = {5.5,5.5,5.5});
  //--- MaxChi2 Hit <11 Hits + Y(0) +  Y(zRef) cut
  //--- N Layers < 11 (9,10) outliers removed if MaxChi2Hit< maxChi2Hits_less11Hit
  declareProperty( "maxChi2Hits_less11Hit" , m_maxChi2HitFullFitLow = {2.5,2.5,2.5});
  //--- If N fired Layers < 11: kill tracks having y(z=0) >50. mm
  declareProperty( "maxYatZeroLow"         , m_maxY0Low             ={50.,60.,70.});
  //--- If N fired Layers < 11 : kill tracks having y(zRef)>500. mm
  declareProperty( "maxYatzRefLow"         , m_maxYZrefLow          ={400.,550.,700.});
  //--- Global Chi2PerDoF
  declareProperty( "maxChi2PerDoF"         , m_maxChi2PerDoF        = {4.0,6.0,7.0});
}
//=============================================================================
// Destructor
//=============================================================================
PrHybridSeeding::~PrHybridSeeding() {}

//=============================================================================
// Initialization
//=============================================================================
StatusCode PrHybridSeeding::initialize() {
 
#ifdef MVA_used
  StatusCode sc = GaudiTupleAlg::initialize();
#ifdef MVA_classification_inFlagHits


  std::vector<std::string> inputVars;

 
  inputVars.push_back("chi2");
  inputVars.push_back("nx");
  inputVars.push_back("ny");
  inputVars.push_back("X0Back");
  inputVars.push_back("worstHitChi2");



  MLPReader0 =  new ReadMLP0(inputVars);
  #ifndef ntuplesGeneration_inFlagHitsCase1
  MLPReader1 =  new ReadMLP1(inputVars);
  #endif

  //5 variables                                                                                                                                                                                            
  for(int i=0;i<5;i++)
    MLPReader_input.push_back(0.0);

#endif

#else
  StatusCode sc = GaudiAlgorithm::initialize();
#endif

  if( sc.isFailure() ) return sc;  // error printed already by GaudiAlgorithm
  if(msgLevel(MSG::DEBUG) ) info() << "==> Initialize" << endmsg;
  m_hitManager= tool<PrHitManager>( m_hitManagerName );
  m_magFieldSvc = svc<ILHCbMagnetSvc>("MagneticFieldSvc", true);
  m_momentumTool = tool<ITrackMomentumEstimate>(m_momentumToolName);
  m_zReference = StateParameters::ZMidT;
  m_hitManager->buildGeometry();
  m_debugTool = 0;
  if ( "" != m_debugToolName ) {
    m_debugTool = tool<IPrDebugTool>( m_debugToolName );
    info()<<"Debug tool "<<m_debugToolName<<" loaded."<<endmsg;
  }else{
    m_wantedKey = -100;  
  }
  if(UNLIKELY(m_doTiming)){
    m_timerTool = tool<ISequencerTimerTool>( "SequencerTimerTool/Timer", this );
    m_timeTotal   = m_timerTool->addTimer( "PrSeeding total" );
    m_timerTool->increaseIndent(); //+1
    m_timeFromForward = m_timerTool->addTimer( "Time from Forward" );
    m_timerTool->increaseIndent(); //+2    
    //Case 0
    m_timeXProjeUp[0] = m_timerTool->addTimer( "Case 0 : X Projections Up");
    m_timerTool->increaseIndent(); //+3
    m_timeCloneXUp[0] = m_timerTool->addTimer( "Case 0 : Clones X Up");
    m_timerTool->increaseIndent(); //+4
    m_timeStereoUp[0] = m_timerTool->addTimer( "Case 0 : AddStereo Up");
    m_timerTool->increaseIndent(); //+5
    m_timeFlagUp[0] = m_timerTool->addTimer(   "Case 0 : Flag Up");
        
    //--- Case 1
    m_timerTool->decreaseIndent(); //+4
    m_timerTool->decreaseIndent(); //+3
    m_timerTool->decreaseIndent(); //+2
    m_timeXProjeUp[1] = m_timerTool->addTimer( "Case 1 : X Projections Up");
    m_timerTool->increaseIndent(); //+3
    m_timeCloneXUp[1] = m_timerTool->addTimer( "Case 1 : Clones X Up");
    m_timerTool->increaseIndent(); //+4
    m_timeStereoUp[1] = m_timerTool->addTimer( "Case 1 : AddStereo Up");
    m_timerTool->increaseIndent(); //+5
    m_timeFlagUp[1] = m_timerTool->addTimer(   "Case 1 : Flag Up");
    m_timerTool->decreaseIndent(); //+4
    m_timerTool->decreaseIndent(); //+3
    m_timerTool->decreaseIndent(); //+2

    //--- Case 2
    m_timeXProjeUp[2] = m_timerTool->addTimer( "Case 2 : X Projections Up");
    m_timerTool->increaseIndent(); //+3
    m_timeCloneXUp[2] = m_timerTool->addTimer( "Case 2 : Clones X Up");
    m_timerTool->increaseIndent(); //+4
    m_timeStereoUp[2] = m_timerTool->addTimer( "Case 2 : AddStereo Up");
    m_timerTool->increaseIndent(); //+5
    m_timeFlagUp[2] = m_timerTool->addTimer(   "Case 2 : Flag Up");
    m_timerTool->decreaseIndent(); //+4
    m_timerTool->decreaseIndent(); //+3
    m_timerTool->decreaseIndent(); //+2
        
    //--- Case 0 Down
    m_timeXProjeDo[0] = m_timerTool->addTimer( "Case 0 : X Projections Down");
    m_timerTool->increaseIndent(); //+3
    m_timeCloneXDo[0] = m_timerTool->addTimer( "Case 0 : Clones X Down");
    m_timerTool->increaseIndent(); //+4
    m_timeStereoDo[0] = m_timerTool->addTimer( "Case 0 : AddStereo Dowm");
    m_timerTool->increaseIndent(); //+5
    m_timeFlagDo[0] = m_timerTool->addTimer(   "Case 0 : Flag Down");
    m_timerTool->decreaseIndent(); //+4
    m_timerTool->decreaseIndent(); //+3
    m_timerTool->decreaseIndent(); //+2
        
    //--- Case 1 Down
    m_timeXProjeDo[1] = m_timerTool->addTimer( "Case 1 : X Projections Down");
    m_timerTool->increaseIndent();//+3
    m_timeCloneXDo[1] = m_timerTool->addTimer( "Case 1 : Clones X Down");
    m_timerTool->increaseIndent();//+4
    m_timeStereoDo[1] = m_timerTool->addTimer( "Case 1 : AddStereo Down");
    m_timerTool->increaseIndent();//+5
    m_timeFlagDo[1] = m_timerTool->addTimer(   "Case 1 : Flag Down");
    m_timerTool->decreaseIndent();//+4
    m_timerTool->decreaseIndent();//+3
    m_timerTool->decreaseIndent();//+2
        
    //--- Case 2 Down
    m_timeXProjeDo[2] = m_timerTool->addTimer( "Case 2 : X Projections Down");
    m_timerTool->increaseIndent(); //+3
    m_timeCloneXDo[2] = m_timerTool->addTimer( "Case 2 : Clones X Down");
    m_timerTool->increaseIndent(); //+4
    m_timeStereoDo[2] = m_timerTool->addTimer( "Case 2 : AddStereo Down");
    m_timerTool->increaseIndent(); //+5
    m_timeFlagDo[2] = m_timerTool->addTimer(   "Case 2 : Flag Down");
    m_timerTool->decreaseIndent(); //+4
    m_timerTool->decreaseIndent(); //+3
    m_timerTool->decreaseIndent(); //+2
    m_timerTool->decreaseIndent();//+1
    m_timeClone[0] = m_timerTool->addTimer(  "RemoveClones Up");
    m_timeClone[1] = m_timerTool->addTimer(  "RemoveClones Down");
    m_timeRecover = m_timerTool->addTimer(   "RecoverTrack");
    m_timeConvert[0] = m_timerTool->addTimer("Convert Tracks Up");
    m_timeConvert[1] = m_timerTool->addTimer("Convert Tracks Down");
    m_timerTool->decreaseIndent();//0
  }
    
  if(UNLIKELY(m_printSettings)){
    info()<<" ===========================   The algorithm is configured to  ============================= "<< endmsg;
    
    if( m_FlagHits){ 
      info()<<"Flag the Hits" << endmsg;
    }else{
      info()<<"Not Flag the Hits"<<endmsg;
    }
    if( m_removeFlagged){
      info()<<"Will Not re-use Flagged"<<endmsg;
    }else{
      info()<<"Will re-use all flagged hits"<<endmsg;
    }
    
    
    if( m_inputName == ""){ 
      info()<<"Standalone Seeding"<<endmsg;
    }else{
      info()<<"Forward tracking tracks as input"<<endmsg;
    }
    if( m_removeClones){
      info()<<"Will Remove Clones"<<endmsg;
    }else{
      info()<<"Will not Remove Clones"<<endmsg;
    }
    if( m_xOnly){
      info()<<"Will use Only X Layers" <<endmsg;
    }
    info()<<" =========================================================================================== "<<endmsg;
    
    if( m_useFix && !m_xOnly) info()<<"WIll Use the triangle Fix for stereo layers"<<endmsg;
    
    info() << "==================================================="<<endmsg
           << "=================GLOBAL SETTINGS==================="<<endmsg
           << " InputName                                      = "<<  m_inputName               << endmsg
           << " OutputName                                     = "<<  m_outputName              << endmsg
           << " HitManagerName                                 = "<<  m_hitManagerName          << endmsg
           << " XOnly                                          = "<<  m_xOnly                   << endmsg
           << " MinXPlanes                                     = "<<  m_minXPlanes              << endmsg
           << " MaxNHits                                       = "<<  m_maxNHits                << endmsg
           << " NCases                                         = "<<  m_nCases                  << endmsg
           << " TimingMeasurement                              = "<<  m_doTiming                << endmsg
           << " dRatio                                         = "<<  m_dRatio                  << endmsg
           <<"  Use dRatioCorrection with XY position ?          "<<  m_useCorrPos              << endmsg
           << " CConst BackProj                                = "<<  m_ConstC                  << endmsg
           << "=================ADD UV Layer Settings============="                             << endmsg
           << "doAsymmUV                                       = "<<  m_doAsymmUV                << endmsg
           << "Use Triangle Fix                                = "<<  m_useFix                  << endmsg
           << "Use SecondOrder `Triangle Fixing                = "<<  m_useFix2ndOrder          << endmsg
           << "*******Fit Settings"<<endmsg
           << "==================Clone Removal and Flag Hits Settings=========="                << endmsg
           << "Remove Clones after X searh       ?             = "<<  m_removeClonesX<<endmsg
           << "Remove Clones after add stereo UV ?             = "<<  m_removeClones<<endmsg
           << "Fraction of common hits to clone kill           = "<<  m_fracCommon <<endmsg
           << "Flag the hits                                   = "<<  m_FlagHits <<endmsg
           << "Remove All Flagged                              = "<<  m_removeFlagged<<endmsg;
    info()<<"=======  Will Run "<<m_nCases<<"  Cases  ===========" <<endmsg;
    
    info()<<" ============================ Find X Projection Settings ======================="<<endmsg;
    info()<<" 1 ==== First-Last Layer Settings  (2 Hit Combo)  "<<endmsg;
    for(unsigned int kk =0;m_alphaCorrection.size()>kk;kk++){
      info()<<"\t Case "<<kk<<"   Rotation                  ="<<m_alphaCorrection[kk]<<endmsg;
      info()<<"\t Case "<<kk<<"   Tolerance after Rotation  ="<<m_TolFirstLast[kk]<<endmsg;
    }
    info()<<" 2 ==== Add Hit in T2 - XLayers Stations  (3 Hit Combo) "<<endmsg;
    info()<<" Allow N Parabola Seed Hits = " << m_maxParabolaSeedHits;
    for(unsigned int kk =0;m_x0Corr.size() > kk; kk++)
    {
      info()<<"\t Case "<<kk<<"       x0Rotation = "<<m_x0Corr[kk]<<endmsg;
      info()<<"\t Case "<<kk<<"    x0SlopeChange = "<<m_x0SlopeChange[kk]<<"    Tolerance at x0SlopeChange "<<  m_TolX0SameSign[kk]<<endmsg;
      info()<<"\t Case "<<kk<<"   x0SlopeChange2 = "<<m_x0SlopeChange2[kk]<<"    Tolerance at x0SlopeChange2"<< m_tolX0Oppsig[kk]<<endmsg;
      info()<<"\t Case "<<kk<<"             x0Max = "<<m_x0Cut[kk]<<"    Tolerance Up at x0Max"<< m_tolAtX0Cut[kk]<< " Tolerance Down at x0Max "<<m_tolAtx0CutOppSig[kk]<<endmsg;
    }
    info()<<" 3 ==== Add Hit in remaining Layers (4/5/6 Hit Combo) ---- "<<endmsg;
    for( unsigned int kk =0; m_tolRemaining.size()>kk;kk++){
      info()<<"\t Case"<<kk<<"    Tolerance = "<<m_tolRemaining[kk]<<endmsg;
    }
    info()<<" 4 ==== Fit XZ Projection  "<<endmsg;
    info()<<"\t minXPlanes "<< m_minXPlanes;
    info()<<"\t dRatio " << m_dRatio<<endmsg;
    for( unsigned int kk =0; m_maxChi2HitsX.size()>kk;kk++){
      info()<<"\t Case "<<kk<<"      MaxChi2Hit X Fit    = "<<m_maxChi2HitsX[kk]<<endmsg;
      info()<<"\t Case "<<kk<<"      MaxChi2Track X Fit  = "<<m_maxChi2DoFX[kk]<<endmsg;
    }
    info()<<" Remove Clones X " << m_removeClonesX<< endmsg;
    info()<<" ========================= Add The UV part =================================== "<<endmsg;
    info()<<" ===> Hit Selection pre-Hough Clustering "<<endmsg;
    info()<<"\t Y Min                 "<< m_yMin      <<endmsg
          <<"\t Y Max                 "<< m_yMax      <<endmsg
          <<"\t UseTrFix              "<<m_useFix     <<endmsg
          <<"\t Y Min TrFix           "<< m_yMin_TrFix<<endmsg
          <<"\t Y Max TrFix           "<< m_yMax_TrFix<<endmsg
          <<"\t Asymmetric UV search  "<< m_doAsymmUV <<endmsg
          <<"\t RemoveHole            "<<m_removeHole <<endmsg
          <<"\t Radius Hole           "<<m_radiusHole <<endmsg
          <<"\t TrFix 2nd Order       "<<m_useFix2ndOrder<<endmsg
          <<"\t Max N Hits            "<<m_maxNHits   <<endmsg
          <<endmsg;
    
    info()<<" ===> Case selection Hough clusters and track fit"<<endmsg;
    for( unsigned int kk = 0; m_tolTyOffset.size()>kk; kk++){
      info()<<"========  \t Case "<<kk<<" ======"<<endmsg;
      info()<<" maxTy-minTy < " <<m_tolTyOffset[kk]<<" + "<< m_tolTySlope[kk] <<" * minTy"<<endmsg;
      info()<<"--> Special nTot Hits <11 ----- "<<endmsg;
      info()<<"max   Chi2DoFX + Chi2Line (when <11) hit < "<<m_Chi2LowLine[kk]<<endmsg;
      info()<<"max   y(z=0)                             < "<<m_maxY0Low[kk]<<endmsg;
      info()<<"max   y(zRef)                            < "<<m_maxYZrefLow[kk]<<endmsg;
      info()<<"max   Chi2 Hit (full fit)                < "<<m_maxChi2HitFullFitLow[kk]<<endmsg;
      info()<<"--> Special nTot Hits >10 ----- "<<endmsg;
      info()<<"max   Chi2DoFX + Chi2Line (when >10) hit < "<<m_Chi2HighLine[kk]<<endmsg;
      info()<<"max   Chi2 Hit (full fit)                < "<<m_maxChi2HitFullFitHigh[kk]<<endmsg;
      info()<<"--> Final Chi2DoF                        < "<<m_maxChi2PerDoF[kk]<<endmsg;
    }
    
    info()<<" ====================== Flag Hits ============================== " << endmsg;
    for( unsigned int kk = 0; m_SizeFlag.size()>kk;kk++){
      info()<<" ======= \t Case "<<kk<<" ======="<<endmsg;
      info()<<"\t Case"<<kk<<" Size Flag"<< m_SizeFlag[kk]<<endmsg;
      info()<<"\t Case"<<kk<<" MaxChi2DoF when 11 hits to flag : "<<m_MaxChi2Flag[kk]<<endmsg;
      info()<<"\t Case"<<kk<<" Max Back Projection when 11 hits: "<<m_MaxX0Flag[kk]<<endmsg;
    }
    
    info()<<" =====================  RECOVER TRACK ROUTINE =================== "<<endmsg;
    if(!m_recover){
      info()<<" WILL NOT RECOVER TRACKS"<<endmsg;
    }else{
      info()<<" WILL RECOVER TRACKS"<<endmsg;
      info()<<" All hits found are marked as used, you will recover the x/z projection not promoted if"<<endmsg;
      info()<<" NX = 6 : used  < "<<m_nusedthreshold[0]<<endmsg;
      info()<<" NX = 5 : used  < "<<m_nusedthreshold[1]<<endmsg;
      info()<<" NX = 4 : used  < "<<m_nusedthreshold[2]<<endmsg;
      info()<<" Max Nb of clusters to process in stereo hit add                   : "<<m_recoNCluster<<endmsg;
      info()<<" Min Nb of hits on recovered full track                            : "<<m_recoMinTotHits<<endmsg;
      info()<<" Min Nb of UV hits to attach at x/z proj recovered with 6 hits     : "<<m_recover_minUV[0]<<endmsg;
      info()<<" Min Nb of UV hits to attach at x/z proj recovered with 5 hits     : "<<m_recover_minUV[1]<<endmsg;
      info()<<" Min Nb of UV hits to attach at x/z proj recovered with 4 hits     : "<<m_recover_minUV[2]<<endmsg;
      info()<<" Hough cluster spread in recovering routine                        : "<<m_recoTolTy <<endmsg;
      info()<<" Max Y(0) of the track recovered                                   : "<<m_recomaxY0<<endmsg;
      info()<<" Recover Line Fit >10 hits                                         : "<<m_recoLineHigh<<endmsg;
      info()<<" Recover Line Fit <10 hits                                         : "<<m_recoLineLow<<endmsg;
      info()<<" Recover ChiDoF                                                    : "<<m_recoFinalChi2<<endmsg;
      info()<<" Recover Chi Outlier                                               : "<<m_recoChiOutlier<<endmsg;
    }
  }
  if( m_nCases >3){
    error()<<"Algorithm does not support more than 3 Cases"<<endmsg;
    return StatusCode::FAILURE;
  }
  return StatusCode::SUCCESS;
}

//=============================================================================
// Main execution
//=============================================================================
StatusCode PrHybridSeeding::execute() {
#ifdef MVA_used
  using namespace Tuples;
#endif
  if(msgLevel(MSG::DEBUG) ) info() << "==> Execute" << endmsg;
  if(UNLIKELY(m_doTiming) ){
    m_timerTool->start( m_timeTotal );
    m_timerTool->start( m_timeFromForward );
  }
  LHCb::Tracks* result = new LHCb::Tracks();
  put( result, m_outputName );
  m_hitManager->decodeData();
  int multiplicity[24];
  int multiplicityTot = 0;
  // UNFLAG USED flag TO ALL HITS !
  debug()<<"UNFLAGGING ALL HITS"<<endmsg;
  for ( unsigned int zone = 0; m_hitManager->nbZones() > zone; ++zone ) {
    auto itH = m_hitManager->getIterator_Begin( zone);
    auto itEnd = m_hitManager->getIterator_End( zone);
    for(  ; itH!= itEnd; ++itH){
      (*itH).setUnUsed();
      multiplicityTot++;
      multiplicity[zone]++;
    }
  }
  //========================================================
  // Remove Seed segment if we pass the  Forward as Input
  //========================================================
  if ( "" != m_inputName){
    debug()<<"Removing Seed Segment from Forward Tracking"<<endmsg;
    //marking hits in forward tracks as used
    for(int i = 0; i < 24; i++){
      auto itH = m_hitManager->getIterator_Begin(i);
      auto itEnd = m_hitManager->getIterator_End(i);
      for(;itH != itEnd;++itH){
        if( (*itH).isUsedInTrack() ) (*itH).setUsed(); //Tuning done in PrForwardTool //TODO make this smarter?
      }
    }
    //convert forward track- T station segments to seed tracks and push them in the final container directly
    LHCb::Tracks* forward = get<LHCb::Tracks>( m_inputName );
    for ( LHCb::Tracks::iterator itT = forward->begin(); forward->end() != itT; ++itT ) {
      if( (*itT)->info( LHCb::Track::PatQuality, 0. ) > 0.8) continue;  //quality from 0-1 (good-bad)
      std::vector<LHCb::LHCbID> ids;
      ids.reserve(20);
      for ( std::vector<LHCb::LHCbID>::const_iterator itId = (*itT)->lhcbIDs().begin();
            (*itT)->lhcbIDs().end() != itId; ++itId ) {
        if( (*itId).isFT() && (*itId).ftID().layer() < 12 ){
          ids.push_back( *itId );
        }
      }
      LHCb::Track* seed = new LHCb::Track;
      seed->setLhcbIDs( ids );
      seed->setType( LHCb::Track::Ttrack );
      seed->setHistory( LHCb::Track::PrSeeding ); //TODO give special flag that it is already a forward (thus long) track!
      seed->setPatRecStatus( LHCb::Track::PatRecIDs );
      seed->addToStates( (*itT)->closestState( 9000. ) );
      result->insert( seed );
    }
  }
  
  //==========================================================
  // Hits are ready to be processed
  //==========================================================
  
  // up / down ( for m_trackCandidates[0/1] )
  for(unsigned int i = 0; i< m_trackCandidates.size() ; i++){
    m_xCandidates[(int)i].clear();
    m_trackCandidates[(int)i].clear();
  }
  if( UNLIKELY(m_doTiming) ) {
    m_timerTool->stop( m_timeFromForward );
  }

  //========================================================
  //------------------MAIN SEQUENCE IS HERE-----------------
  //========================================================
  //----- Loop through lower and upper half
  for( unsigned int part= 0; 2 > part; ++part ){
    //Swap them ? externally icase & inner loop part? m_xCandidates
    for(unsigned int icase = 0; m_nCases > icase ; ++icase){
      if(UNLIKELY(m_doTiming)){
        if( part ==0){
          m_timerTool->start( m_timeXProjeUp[icase]);
        }else{
          m_timerTool->start( m_timeXProjeDo[icase]);
        }
      }
      m_xCandidates[(int)part].clear(); //x candidates up cleaned every Case!
      findXProjections(part,icase);
      if(UNLIKELY(m_doTiming)){
        if( part ==0){
          m_timerTool->stop( m_timeXProjeUp[icase]);
          m_timerTool->start( m_timeCloneXUp[icase]);
        }else{
          m_timerTool->stop( m_timeXProjeDo[icase]);
          m_timerTool->start(m_timeCloneXDo[icase]);
        }
      }
      if(m_removeClonesX) removeClonesX( part, icase, m_xOnly);
      if( UNLIKELY(m_doTiming)){ 
        if( part==0 ){
          m_timerTool->stop( m_timeCloneXUp[icase]);
          m_timerTool->start(m_timeStereoUp[icase]);
        }else{
          m_timerTool->stop( m_timeCloneXDo[icase]);
          m_timerTool->start(m_timeStereoDo[icase]);
        }
      }
      
      //Add The stereo Part
      if(!m_xOnly){
        addStereo( part, icase );
      }
      if(UNLIKELY(m_doTiming)){
        if(part == 0){
          m_timerTool->stop( m_timeStereoUp[icase]);
          m_timerTool->start(  m_timeFlagUp[icase]);
        }else{
          m_timerTool->stop( m_timeStereoDo[icase]);
          m_timerTool->start( m_timeFlagDo[icase]);
        }
      }
      //Flag found Hits at the end of each single case ( exclude the latest one )
      //icase = 0,1,2
      //m_nCase = 1,2,3
      if(m_FlagHits && (icase +1 < m_nCases)   && !m_xOnly){
        flagHits(icase,part);
      }
      if(UNLIKELY(m_doTiming)){
        if( part ==0 ){
          m_timerTool->stop( m_timeFlagUp[icase]);
        }else{
          m_timerTool->stop(m_timeFlagDo[icase]);
        }
      }
    }
  }
  for( unsigned int part = 0 ; part<2; part ++){
    if(UNLIKELY(m_doTiming)){
      m_timerTool->start( m_timeClone[(int)part]);
    }
    //Clone removal ( up/down )
    removeClones( part );
    if(UNLIKELY(m_doTiming)){
      m_timerTool->stop( m_timeClone[ (int) part]);
    }
  }
  if( m_doTiming){
    m_timerTool->start(m_timeRecover);
  }
  if( m_recover){
    m_xCandidates[0].clear();
    m_xCandidates[1].clear();
    RecoverTrack();
  }
  if(UNLIKELY(m_doTiming)){
    m_timerTool->stop( m_timeRecover);
  }
  for( unsigned int part=0 ; part<2; part++){
    if(UNLIKELY(m_doTiming)){
      m_timerTool->start(m_timeConvert[(int)part]);
    }
    //Convert LHCb tracks
    makeLHCbTracks(result, part);
    if(UNLIKELY(m_doTiming)){
      m_timerTool->stop( m_timeConvert[(int)part]);
    }
  }
  if(msgLevel(MSG::DEBUG)) debug()<<"Making LHCb Tracks Done"<<endmsg;
  if(UNLIKELY(m_doTiming)){
    m_timerTool->stop( m_timeTotal);
  }
  return StatusCode::SUCCESS;
}

//=============================================================================
//  Finalize
//=============================================================================
StatusCode PrHybridSeeding::finalize() {
  if( msgLevel(MSG::DEBUG) ) debug() << "==> Finalize" << endmsg;
  #ifdef MVA_used
#ifdef MVA_classification_inFlagHits
  if(MLPReader0 != nullptr){
    delete MLPReader0;
    MLPReader0 = nullptr;
  }
#ifndef ntuplesGeneration_inFlagHitsCase1
  if(MLPReader1 != nullptr){
    delete MLPReader1;
    MLPReader1 = nullptr;
  }
#endif

#endif
  return GaudiTupleAlg::finalize();
#else
  return GaudiAlgorithm::finalize();
#endif
}

inline PrHitZone& PrHybridSeeding::GetZone( unsigned int kk){
  return m_hitManager->zone(kk);
}


//void PrHybridSeeding::CollectUV(  PrHybridSeedTrack& xProje,const std::vector<PrHitZone*>& uvZones, PrHits& myStereo){
void PrHybridSeeding::CollectUV(  PrHybridSeedTrack& xProje,const std::vector<unsigned int>& uvZones, PrHits& myStereo){
  myStereo.clear();
  unsigned int part = xProje.zone();
  for( const unsigned int& kk : uvZones){
    float yMin = -1.*std::fabs(m_yMax);
    float yMax = +1.*std::fabs(m_yMax);
    float dxDy = GetZone(kk).dxDy();
    float zPlane = GetZone(kk).z();

    float xPred = xProje.x(zPlane);
    float recDxDyZ = 1./(dxDy*zPlane);
    float recDxDy = 1./dxDy;
    if(yMin >yMax){
      std::swap( yMin, yMax);
    }
    if(part==0){
      //Upper tracks ( expect y>0 for them )
      yMin = m_yMin;
      //[-1,2500] mm
      if(kk%2==1 && m_useFix){
        //Upper track but lower module looking at
        //y looked in [-2,30] mm
        yMin = +m_yMin_TrFix;
        yMax = +m_yMax_TrFix;
      }
    }else{
      //Lower tracks ( expected y < 0 ( m_yMin need to be negative, m_yMax need to be positive) )
      yMax = -m_yMin;
      //[-2500,+1] mm
      if(kk%2==0 && m_useFix){
        //Lower track but upper module looking at
        //y looked in [-30,+2] mm
        yMin = -m_yMax_TrFix;
        yMax = -m_yMin_TrFix;
      }
    }
    // Compute the [xMin,xMax] value in the rotated reference system based on y tolerance!
    // Due to dxDy being positive or negative ( yMin here is always < yMax ) you can have them swapped
    // x(y) ~ xPred = x(0) + dxDy*y => x(0) = xPred - dxDy*y ; xMin xMax from yMin and yMax and 
    // swap them if opposite in sign
    float xMin = xPred - yMin*dxDy;
    float xMax = xPred - yMax*dxDy;
    if(xMin > xMax){
      std::swap( xMax, xMin);      //revert to properly take into account dxDy
    }
    auto itH = m_hitManager->getIterator_lowerBound( kk   , xMin);
    auto EndLay = m_hitManager->getIterator_End( kk);
    auto itEnd = m_hitManager->get_upperBound_log( itH, EndLay, xMax);
    for( ; itEnd!= itH; ++itH ){
      if( (*itH).isUsed() && m_removeFlagged ) continue;
      float y = ( xPred - (*itH).x() )*recDxDy;
      if(y >yMax || y <yMin ) continue; //may can skip this if
      if(m_useFix2ndOrder){
        if(part == 0){
          //Up track
          if(kk%2 == 0){
            //Up track Up modules
            if( (*itH).yMin()<0. ){ //leaking area to y negative in up modules
              yMin = m_yMin; yMax = m_yMax;   //[-1.0, 2500]
            }else{//y positive but cut away
              yMin = -2.0 + std::fabs((*itH).yMin());yMax = m_yMax;}
          }else{//kk%2 ==1
            //Up track Down modules
            if( (*itH).yMax()<0. ){ continue;
            }else{
              yMin = -1.0 ;yMax = +2.0 + std::fabs((*itH).yMax());
            }
          }
        }else{
          //Down track
          if(kk%2 == 1){
            //Down track Down modules
            if( (*itH).yMax()>0.){ //leaking area to y positive in down module
              yMin = -std::fabs(m_yMax) ; yMax = +1.0; //[-2500, 1.0]
            }else{//y negative but cut away
              yMin = -std::fabs(m_yMax) ;yMax = +2.0-std::fabs((*itH).yMax()); // [ -2500, 2.0 ]
            }
          }else{//kk%2 ==1
            //Down track Up modules
            if( (*itH).yMin()<0.){
              yMin = -2.0 - std::fabs((*itH).yMin());yMax = +1.0;
            }else{
              continue;
            }
          }
        }
        if(msgLevel(MSG::DEBUG)){
          debug()<<"y Hit computed    "<< y<<endmsg;
        }
      }
      if(y>yMax || y< yMin) continue;
      if(m_removeHole){
        float radius =  y*y + xPred*xPred;
        if(  radius < m_radiusHole*m_radiusHole) continue;
      }
      //Absolute value!
      if( part == 0 && y<0.){
        (*itH).setCoord( -1.*std::fabs((xPred-(*itH).x())*recDxDyZ));
      }else if( part == 1 && y>0.){
        (*itH).setCoord( -1.*std::fabs((xPred-(*itH).x())*recDxDyZ));
      }else{
        //part ==1 & y<0 or part ==0 & y>0
        (*itH).setCoord( std::fabs( (xPred - (*itH).x() )*recDxDyZ ) );
      }
      myStereo.push_back( &*itH );
    }
  }
}

void PrHybridSeeding::addStereo( unsigned int part, unsigned int iCase){
  if(msgLevel(MSG::DEBUG)){ debug()<<"Adding stereo hits for part = "<<part<<" and Case = "<<iCase<<endmsg;
  }
  //reco is used here to switch on the track recovering configurables
  bool reco = false;
  if(iCase == 4){
    reco = true;
    iCase =2;
  }
  std::vector<unsigned int> uvZones;
  if(m_useFix){  uvZones.reserve(12); 
  }else{         uvZones.reserve(6); }
  //---- Load the u/v layers to look for hits
  if(m_useFix){  
    for( unsigned int detector_part :{0,1}){
      for( unsigned int uvZoneID :{ s_T1U, s_T1V,s_T2U,s_T2V,s_T3U,s_T3V}){
        if(msgLevel(MSG::DEBUG)){ debug()<<"Push back zone = "<<uvZoneID+detector_part<<endmsg;
          if( GetZone(uvZoneID|detector_part).isX()) always()<<"ERROR LOAD"<<endmsg;
        }
        uvZones.push_back( uvZoneID | detector_part );
      }
    }
  }else{
    for( unsigned int uvZoneID :{s_T1U, s_T1V,s_T2U,s_T2V,s_T3U,s_T3V}){
      uvZones.push_back( uvZoneID | part );
    }
  }
  PrHits myStereo;
  myStereo.reserve(400);
  PrPlaneHybridCounter plCount;
  int maxUVfound = 0;
  PrLineFitterY BestLine( m_zReference );
  float initCoord = std::numeric_limits<float>::max();
  float ClusterSpread[3][DEPTH_HOUGH];
  HitIter initBeg[3][DEPTH_HOUGH];
  HitIter initEnd[3][DEPTH_HOUGH];
  
  //TODO : HoughUtils from Manuel Shiller could be used here
  
  if(msgLevel(MSG::DEBUG)){ 
    debug() << " Will loop on x-candidates case "<<iCase <<" Part = "<<part<<endmsg;
  }
  for( auto itT =m_xCandidates[part].begin();m_xCandidates[part].end() != itT; ++itT){
    if( !(*itT).valid() ) continue;
    CollectUV( (*itT) , uvZones , myStereo);
    std::sort( myStereo.begin(), myStereo.end(), PrHit::LowerByCoord() );//from Bigger to smaller
    //Define the min UV hits based on the amount of hits of the track
    unsigned int minUV = 4;
    if((*itT).size()==6){
      if(!reco){ minUV = m_minUV6[iCase];
      }else{
        minUV = m_recover_minUV[0];
      }
    }else if((*itT).size()==5){
      if(!reco){  minUV = m_minUV5[iCase];
      }else{
        minUV = m_recover_minUV[1];
      }
    }else if((*itT).size()==4){
      if(!reco){ minUV = m_minUV4[iCase];
      }else{
        minUV = m_recover_minUV[2];
      }
    }
    maxUVfound=-1;
    unsigned int firstSpace = m_trackCandidates[(int)part].size();
    //initialize the clusters ( here 9 )
    if(msgLevel(MSG::DEBUG)) always()<<" initializing clusters "<<endmsg;
    for( int i =0; i<3 ; ++i){
      for( int j = 0; j<DEPTH_HOUGH; ++j){
        ClusterSpread[i][j] = initCoord;
        initBeg[i][j] = myStereo.end();
        initEnd[i][j] = myStereo.end();
      }
    }
    PrHits::iterator itBeg = myStereo.begin();
    PrHits::iterator itEnd =  itBeg+4<=myStereo.end() ? itBeg+4: myStereo.end();
    if(msgLevel(MSG::DEBUG)){ debug()<<"Storing clusters"<<endmsg;}
    for( ; itEnd-itBeg>=4; IncreaseIters( itBeg, itEnd,4, myStereo)){
      // always steps of 4 as minimum!
      float tolTy = m_tolTyOffset[iCase] + m_tolTySlope[iCase]*(*itBeg)->coord();
      if( reco ){  tolTy = m_recoTolTy; } //different tolTy when you are in the recovering routine
      for (int i = 0; i < 3; ++i){
        //3 because you have clusters with 4/5/6 hits (more  u/v layers ? CHANGE IT!)
        auto itLast = itEnd+2-i <= myStereo.end() ? itEnd+2-i : myStereo.end();
        // spread6 computed first
        float spread = (*(itLast-1))->coord() - (*itBeg)->coord();
        if(CheckCluster(itBeg, itLast, tolTy)){
          //cluster is fine check
          int j;
          // find the cluster worse than the current cluster
          for (j = 0; j < DEPTH_HOUGH; ++j){
            if(spread < ClusterSpread[i][j]) break;
          }
          if (j == DEPTH_HOUGH) continue;
          // Shift clusters after the found position
          for(int k = DEPTH_HOUGH-1; k > j; --k){
            initBeg[i][k] = initBeg[i][k-1];
            ClusterSpread[i][k] = ClusterSpread[i][k-1];
            initEnd[i][k] = initEnd[i][k-1];
          }
          // insert current cluster at found position
          initBeg[i][j] = itBeg;
          ClusterSpread[i][j] = spread;
          initEnd[i][j] = itLast;
          break;
        }
      }
    }
    //you have a minimal number of total hits to find on track dependent if standard 3 cases or from Recover routine
    unsigned int minTot = m_minTot[iCase];
    if(reco) minTot = m_recoMinTotHits;
    if( minTot <9 ) minTot= 9; 
    //WARNING : protection for the algorithm , do not allow the user to require less than 9 hits.
    if(msgLevel(MSG::DEBUG)) info()<<"Cluster stored"<<endmsg;
    int nXinit = (*itT).size();
    BestLine.setXProj((*itT));
    int nCandidates = 0;
    int maxNClusters = m_maxNClusters[iCase];
    float minChi2LineLow = m_Chi2LowLine[iCase];
    float minChi2LineHigh = m_Chi2HighLine[iCase];
    float maxChi2DoF = m_maxChi2PerDoF[iCase];
    if( reco){
      minTot = m_recoMinTotHits;
      maxNClusters = m_recoNCluster; //max N of clusters to process changed when in reco case.
      minChi2LineLow = m_recoLineLow;
      minChi2LineHigh = m_recoLineHigh;
      maxChi2DoF = m_recoFinalChi2;
    }
    for( int i = 0; i<3; ++i){
      if( maxUVfound == 6 && !reco) break;
      for( int j = 0; j<DEPTH_HOUGH; ++j){
        if( maxUVfound == 6 && !reco) break;
        if(msgLevel(MSG::DEBUG)){ debug()<<"Processing cluster i ="<<i<<"  j = "<<j<<endmsg; }
        //you are at the cluster Number (1+j) +(3*i)
        //By construction the first clusters you process are the ones with smaller spread
        PrHits::iterator itBeg = myStereo.begin();
        PrHits::iterator itEnd = myStereo.end();
        int nLay = 0;
        itBeg = initBeg[i][j];
        itEnd = initEnd[i][j];
        if( itBeg == myStereo.end() && itEnd == myStereo.end()){ break;}
        ExtendCluster( itBeg, itEnd , iCase, myStereo , nLay ,reco);
        //maxUVfound are the minimal UV layers to find since track quality is based on nHits
        if( nXinit + nLay < (int)minTot ) continue;
        //if you already have less than minTot layers, look to next "extended" cluster"
        if( !reco){  //recovering track routine accept also same size hits!
          if(  (int)nLay <= (int)maxUVfound || (int)nLay < (int)minUV )  continue; //<= if at the end you have only >= maxUVfound.
        }else{ 
          if(  (int)nLay < (int)maxUVfound || (int)nLay < (int)minUV ) continue; //<= in reco you allow to find same size tracks.
        }
        //You increased the clusters and check to have at least minUV hits and number of fired layers in current clusters can kill current candidates produced
        nCandidates++;
        if( nCandidates > maxNClusters ) break;
        //Set the plane counter
        plCount.set( itBeg , itEnd );
        if(UNLIKELY( msgLevel(MSG::DEBUG))){ 
          debug()<<" Plane counter diff UV = "<< plCount.nbDifferentUV() <<" n Lay from extending cluter   "<< nLay<<endmsg;
        }
        nLay = plCount.nbDifferentUV();
        if( nLay == 4 && !plCount.isOKUV()){
          continue;
        }
        BestLine.setXProj( (*itT) );
        PrHybridSeedTrack temp(*itT);
        bool fit = false;
        if( plCount.nbSingleUV() == plCount.nbDifferentUV()){
          if(msgLevel(MSG::DEBUG)){ 
            info()<<"Fitting line when same amount of UV in Planes"<<endmsg;
          }
          fit = BestLine.fit( itBeg,itEnd );
          if( fit && LineOK( minChi2LineLow, minChi2LineHigh, BestLine, temp) ){
            for( auto itH = itBeg; itEnd != itH; ++itH){
              temp.addHit( (*itH) );
            }
            temp.setYParam( BestLine.ay0(), BestLine.by());
          }
        }else{ 
          // More than one hit per layer in the cluster
          // Remove outliers from multiple hit in layers is forced here
          fit = BestLine.fit2( itBeg, itEnd, plCount );
          for( auto itH = itBeg; itEnd!= itH; ++itH){
            temp.addHit( (*itH) );
          }
          while( BestLine.nHitsLine() > plCount.nbDifferentUV() ){
            BestLine.fit2( temp.hits().begin(), temp.hits().end() , plCount);
            if(msgLevel(MSG::DEBUG)){ debug()<<" Removing outliers"<<endmsg; }
            auto worstMultiple = temp.hits().end();
            float worstchi = std::numeric_limits<float>::min();
            float chi_hit = std::numeric_limits<float>::min();
            for( auto line_hit = temp.hits().begin(); line_hit!=temp.hits().end(); ++line_hit){
              if( (*line_hit) ->isX() ) continue;
              if( plCount.nbInPlane( (*line_hit)->planeCode() ) < 2) continue;
              chi_hit = BestLine.chi2hit( (*line_hit));
              if( chi_hit > worstchi ){
                worstchi = chi_hit ;
                worstMultiple = line_hit;
              }
            }
            plCount.removeHitInPlane( (*worstMultiple)->planeCode() );
            //it will just change the m_planeCode[int] not the whole internal hit counting ! )
            temp.hits().erase( worstMultiple );
            fit = BestLine.fit2( temp.hits().begin(), temp.hits().end(), plCount);
          }
          //fit status is the Success of the fit ( no requirements on Chi2 )
          if( fit ){
            fit = LineOK( minChi2LineLow, minChi2LineHigh, BestLine, (*itT));
            temp.setYParam( BestLine.ay0(), BestLine.by());
          }
        }
        if( temp.size() < minTot || !fit ) continue;
        bool ok = false;
        temp.setnXnY( (*itT).size(), plCount.nbDifferentUV());
        ok = fitSimultaneouslyXY(temp,iCase);
        while(!ok && temp.size() > minTot ){
          ok = removeWorstAndRefit( temp, iCase);
          if( temp.size()>m_maxNHits) ok = false;
        }
        if( ok  && temp.size() >= minTot){
          setChi2(temp);
          if( temp.chi2PerDoF() < maxChi2DoF ){
            temp.setCase(iCase);
            int nUVfound = temp.size()-nXinit;
            if( nUVfound > maxUVfound && !reco){
              if( (
                   ( temp.size() < 11 &&
                     (std::fabs(temp.y(0.))< m_maxY0Low[iCase] &&
                      std::fabs(temp.y(m_zReference)) < m_maxYZrefLow[iCase] )
                     ) ||
                   ( temp.size() > 10 ) ) ){
                if(m_removeClones){
                  std::sort(temp.hits().begin(), temp.hits().end(), compLHCbID());
                }
                maxUVfound = (int)temp.size()-nXinit;
                m_trackCandidates[(int)part].push_back( temp );
              }
            }else if( m_recover && reco && nUVfound >= maxUVfound &&  std::fabs( temp.y(0.)) < m_recomaxY0  && temp.size() >= minTot){
              maxUVfound = (int)temp.size()-nXinit;
              m_trackCandidates[(int)part].push_back( temp );
              //ordered by spread take one small spread passing selections.
            }
          }
        }
      }
    }
    if( m_recover && !reco && m_trackCandidates[(int)part].size() == firstSpace ){
      m_trackToRecover[(int)part].push_back( (*itT) );
    }
    //if in reco case keep them all.
    if(  m_trackCandidates[(int)part].size() > firstSpace+1  ){
      for( unsigned int kk = firstSpace; m_trackCandidates[(int)part].size() > kk -1 ; ++kk ){
        if(!m_trackCandidates[(int)part][kk].valid()) continue;
        for( unsigned int ll = kk + 1; m_trackCandidates[(int)part].size() > ll; ++ll ){
          if(!m_trackCandidates[(int)part][ll].valid())continue;
          if(m_trackCandidates[(int)part][ll].size() < m_trackCandidates[(int)part][kk].size()){
            m_trackCandidates[(int)part][ll].setValid( false );
          }
          else if ( m_trackCandidates[(int)part][ll].size() > m_trackCandidates[(int)part][kk].size()){
            m_trackCandidates[(int)part][kk].setValid( false );
          }
          else if( m_trackCandidates[(int)part][kk].chi2() < m_trackCandidates[(int)part][ll].chi2()) {
            m_trackCandidates[(int)part][ll].setValid( false );
          }
          else{
            m_trackCandidates[(int)part][kk].setValid( false );
          }
        }//end loop track Candidates for removal bad tracks
      }//loop candidate removal
    }
  }
}


void PrHybridSeeding::RecoverTrack(){
  for( int part = 0; part< 2; part++){
    for( auto itT1 = m_trackCandidates[(int)part].begin() ; itT1!= m_trackCandidates[(int)part].end() ; ++itT1){
      if(!(*itT1).valid()) continue;
      for( auto itH = (*itT1).hits().begin(); itH!= (*itT1).hits().end() ; ++itH){
        (*itH)->setUsed();
      }
    }
  }
  PrPlaneHybridCounter plCount;
  for( int part = 0; part<2; part++){
    for( auto itT1 = m_trackToRecover[(int)part].begin(); itT1!= m_trackToRecover[(int)part].end() ; ++itT1){
      int nused_threshold = 1;
      bool recover_it = false;
      if( (*itT1).size() == 6){
        nused_threshold = m_nusedthreshold[0];
      }else if( (*itT1).size() == 5){
        nused_threshold = m_nusedthreshold[1];
      }else if( (*itT1).size() == 4){
        nused_threshold = m_nusedthreshold[2];
      }
      int nUsed = plCount.nUsed( (*itT1).hits().begin(), (*itT1).hits().end());
      if( nUsed < nused_threshold ){
        if( nUsed ==0){
          recover_it = true;
        }else if( nUsed > 0){
          while( nUsed !=0){
            auto it = remove_if( (*itT1).hits().begin(), (*itT1).hits().end(), [](const PrHit* hit){ return hit->isUsed();});
            (*itT1).hits().erase(it, (*itT1).hits().end());
            nUsed--;
          }
          recover_it = plCount.isOKXTrack( (*itT1).hits().begin(), (*itT1).hits().end());
          if( recover_it) {
            recover_it = fitXProjection( (*itT1) , 2 );
            setChi2X( (*itT1));
            recover_it = true;
          }
        }
        if(recover_it){
          (*itT1).setRecovered(true);
        }
      }
    }
    std::sort(m_trackToRecover[(int)part].begin(), m_trackToRecover[(int)part].end(), PrHybridSeedTrack::GreaterBySize());
    unsigned int maxCommon =0;
    for( auto itT1 = m_trackToRecover[(int)part].begin(); itT1!= m_trackToRecover[(int)part].end(); ++itT1){
      //maxCommon = (*itT1).size();
      if( !(*itT1).recovered() ) continue;
      maxCommon = std::floor( 0.7*(*itT1).size() ); //70 % of lower number of hits track
      for( auto itT2 = itT1+1; itT2!= m_trackToRecover[(int)part].end(); ++itT2){
        if( !(*itT2).recovered()) continue;
        if( !CloseEnough( (*itT1), (*itT2), 5.0)) continue;
        if( areClones((*itT1), (*itT2), maxCommon) ){
          (*itT1).setRecovered(false);
          break;
        }
      }
      if( (*itT1).recovered() ){
        m_xCandidates[int(part)].push_back( (*itT1));
      }
    }
    addStereo( part, 4);
    m_trackToRecover[(int)part].clear();
  }
}



//=====================================================
//Add Stereo method and support methods
//=====================================================
bool PrHybridSeeding::CheckCluster( PrHits::iterator& itBeg, PrHits::iterator& itEnd, float tol){
  return ( (*(itEnd-1))->coord() - (*itBeg)->coord() < tol);
}

void PrHybridSeeding::IncreaseIters( PrHits::iterator& itBeg, PrHits::iterator& itEnd, unsigned int minUV , PrHits& myStereo){
  ++itBeg;
  itEnd = (itBeg + minUV <=myStereo.end()? itBeg+minUV : myStereo.end());
}
bool PrHybridSeeding::LineOK( float minChi2Low, float minChi2High, PrLineFitterY& line, PrHybridSeedTrack& xProje ){
  const int nHits = line.nHitsLine() + xProje.size();
  //Do some kind of weighting here par[0] par[1] 
  const float Chi2DoFLineXProj = xProje.chi2PerDoF() + line.Chi2DoF();
  //Should check here ( maybe bottleneck for efficiency at large slope in Y )
  if( nHits > 10 && Chi2DoFLineXProj < minChi2High) return true;
  if( nHits < 11 && Chi2DoFLineXProj < minChi2Low ) return true;
  return false;
}

void PrHybridSeeding::removeClonesX(unsigned int part,unsigned int icase, bool xOnly){
  std::sort(m_xCandidates[(int)part].begin(), m_xCandidates[(int)part].end(),PrHybridSeedTrack::GreaterBySize());
  unsigned int maxCommon = 0;
  if( ! m_xOnly ){
    for( auto itT1 = m_xCandidates[part].begin(); m_xCandidates[part].end() !=itT1; ++itT1){
      if(!(*itT1).valid() || xOnly) continue;
      for( auto itT2 = itT1 + 1; m_xCandidates[part].end() !=itT2; ++itT2 ){
        if(!(*itT2).valid())                continue;
        if(((*itT1).bloomfilter() & (*itT2).bloomfilter()).empty()){
          continue;
        }
        if(! CloseEnough( (*itT1), (*itT2), 2.0) ) continue;
        int Compare = (*itT1).size()+(*itT2).size(); // 36 , 30 ,30  ( 6 vs 6 6 vs 5 )
        if( Compare <10 ){
          maxCommon = 1;
        }else if( Compare ==10 ){
          maxCommon = 2;
        }else if( Compare >10){
          maxCommon = 3;
        }
        if( areClones( (*itT1), (*itT2), maxCommon)){
          (*itT1).setValid(false);
          break;
        }
      }
    }
  }else{
    unsigned int maxCommon = 0;
    for( auto itT1 = m_xCandidates[part].begin(); m_xCandidates[part].end() !=itT1; ++itT1){
      if( !(*itT1).valid()) continue;
      for( auto itT2 = itT1 +1; m_xCandidates[part].end()!=itT2; ++itT2){
        if(!(*itT2).valid() ) continue;
        if( CloseEnough( (*itT1), (*itT2), 2.0)) continue;
        int Compare = (*itT1).size()+(*itT2).size();
        if( Compare <10 ){      maxCommon = 1;}
        else if( Compare ==10 ){maxCommon = 2;}
        else if( Compare >10){  maxCommon = 3;}
        if( areClones( (*itT1), (*itT2), maxCommon)){
          (*itT1).setValid(false);
          break;
        }
      }
    }
    if( icase == m_nCases-1){ // xOnly is automatic!
      for( auto itT1 = m_xCandidates[part].begin() ; m_xCandidates[part].end()!= itT1; ++itT1){
        if((*itT1).valid()) m_trackCandidates[(int)part].push_back( (*itT1));
      }
    }
  }
}

void PrHybridSeeding::removeClones( unsigned int part){
  std::sort(m_trackCandidates[(int)part].begin(), m_trackCandidates[(int)part].end(),PrHybridSeedTrack::GreaterBySize());
  //small size in front with larger chi2
  for( auto itT1 = m_trackCandidates[(int)part].begin(); m_trackCandidates[(int)part].end()!=itT1; ++itT1){
    if( !(*itT1).valid() ) continue;
    for ( auto itT2 = itT1 + 1; m_trackCandidates[(int)part].end() !=itT2; ++itT2 ) {
      if( !(*itT2).valid()) continue;
      if(((*itT1).bloomfilter() & (*itT2).bloomfilter()).empty()) continue;
      if( !CloseEnough( (*itT1),(*itT2), 2.0) ) continue;
      unsigned int maxCommon = std::floor( m_fracCommon * std::min( (*itT1).size(), (*itT2).size() ));
      if( areClones((*itT1), (*itT2), maxCommon)){
        (*itT1).setValid(false);
        break;
      }
    }
  }
}
void PrHybridSeeding::flagHits(unsigned int icase, unsigned int part){
  

#ifdef ntuplesGeneration_inFlagHitsCase0
  // TMVA::Reader *reader = new TMVA::Reader();                                                                                                                                                           
  for(auto track_ptr = m_trackCandidates[(int)part].begin(); m_trackCandidates[(int)part].end()!=track_ptr ; ++track_ptr){
    PrHybridSeedTrack track= *track_ptr;
    LHCb::MCParticle* partic = nullptr;
    double efficiency = -1;
    int nAssocHits = -10;
    bool assoc = AssocTrack(track,efficiency,partic,nAssocHits);

    if(track.Case()!=icase) continue;
    if(!(track).valid()) continue;

    //you can also use the track.setWorstHitChi2(maxChi2); (but it really doesn't matter in the ntuples generation phase to reduce timing)                                                                
    auto worst= track.hits().begin();
    float maxChi2 = std::numeric_limits<float>::min();
    float chi_hit = std::numeric_limits<float>::min();
    //  if( track.nx() <4 || track.ny()<4) return false;                                                                                                                                                  
    for( auto itH = track.hits().begin(); track.hits().end() != itH; ++itH ){
      chi_hit = track.chi2( (*itH));
      if( chi_hit  > maxChi2){
	maxChi2 = chi_hit;
	worst = itH;
      }
    }
    float worstHitChi2 = maxChi2;
    Tuple nTupleXY0 = nTuple("XYtracks0","XYtracks0");


    nTupleXY0->column("assoc",assoc);
    nTupleXY0->column("Case",(double)icase);
    nTupleXY0->column("chi2",(double)track.chi2());
    nTupleXY0->column("nx",(double)track.nx());
    nTupleXY0->column("ny",(double)track.ny());
    nTupleXY0->column("X0Back",(double)abs(track.X0()));
    nTupleXY0->column("worstHitChi2",(double)worstHitChi2);
    nTupleXY0->column("nXHits",(double)track.hits().size());

    nTupleXY0->column("xT1",(double)track.xT1());
    nTupleXY0->column("xT2",(double)track.xT2());
    nTupleXY0->column("xT3",(double)track.xT3());

    nTupleXY0->column("ax",(double)track.ax());
    nTupleXY0->column("bx",(double)track.bx());
    nTupleXY0->column("cx",(double)track.cx());
    nTupleXY0->column("ay",(double)track.ay());
    nTupleXY0->column("by",(double)track.by());
    nTupleXY0->column("dRatio",(double)track.dRatio());

    nTupleXY0->write();

    if(assoc)
      {for(auto it = (track).hits().begin();(track).hits().end()!=it; ++it){
	  if( (int)part != (int) (*it)->zone()%2) continue;
	  (*it)->setUsed();}}
  }
#endif

#ifdef MVA_classification_inFlagHits

  if(icase==0){
    float MLP0_scut=0.7; //or 0.48 //based on s/sqrt(s+b)                                                                                                                         
    float MLP0_bcut=0.2;

  for(auto track_ptr = m_trackCandidates[(int)part].begin(); m_trackCandidates[(int)part].end()!=track_ptr ; ++track_ptr){
   

    if((*track_ptr).Case()!=icase) continue;
    if(!(*track_ptr).valid()) continue;


    //MLPReader_input[0]=(double)icase;
    MLPReader_input[0]=(double)(*track_ptr).chi2();
    MLPReader_input[1]=(double)(*track_ptr).nx();
    MLPReader_input[2]=(double)(*track_ptr).ny();
    MLPReader_input[3]=(double)abs((*track_ptr).X0());
    MLPReader_input[4]=(double)(*track_ptr).worstHitChi2();



    float quality = MLPReader0->GetMvaValue(MLPReader_input);

    if(quality>MLP0_scut)
      {
        for(auto it = (*track_ptr).hits().begin();(*track_ptr).hits().end()!=it; ++it){
          if( (int)part != (int) (*it)->zone()%2) continue;
          (*it)->setUsed();}
      }
    else
      if(quality<MLP0_bcut)
	(*track_ptr).setValid(false);
  }
  }
  else //if(icase==1)
    {
#ifdef ntuplesGeneration_inFlagHitsCase1
      //here we generate ntuples for case 1 after using the mva in case 0 
      for(auto track_ptr = m_trackCandidates[(int)part].begin(); m_trackCandidates[(int)part].end()!=track_ptr ; ++track_ptr){
	PrHybridSeedTrack track= *track_ptr;
	LHCb::MCParticle* partic = nullptr;
	double efficiency = -1;
	int nAssocHits = -10;
	bool assoc = AssocTrack(track,efficiency,partic,nAssocHits);

	if(track.Case()!=icase) continue;
	if(!(track).valid()) continue;

	//you can also use the track.setWorstHitChi2(maxChi2); (but it really doesn't matter in the ntuples generation phase to reduce timing)                                                             

	auto worst= track.hits().begin();
	float maxChi2 = std::numeric_limits<float>::min();
	float chi_hit = std::numeric_limits<float>::min();
	//  if( track.nx() <4 || track.ny()<4) return false;                                                                                                                                               

                                                                                                                                                                                                          
	for( auto itH = track.hits().begin(); track.hits().end() != itH; ++itH ){
	  chi_hit = track.chi2( (*itH));
	  if( chi_hit  > maxChi2){
	    maxChi2 = chi_hit;
	    worst = itH;
	  }
	}
	float worstHitChi2 = maxChi2;
	Tuple nTupleXY1 = nTuple("XYtracks1","XYtracks1");


	nTupleXY1->column("assoc",assoc);
	nTupleXY1->column("Case",(double)icase);
	nTupleXY1->column("chi2",(double)track.chi2());
	nTupleXY1->column("nx",(double)track.nx());
	nTupleXY1->column("ny",(double)track.ny());
	nTupleXY1->column("X0Back",(double)abs(track.X0()));
	nTupleXY1->column("worstHitChi2",(double)worstHitChi2);
	nTupleXY1->column("nXHits",(double)track.hits().size());

	nTupleXY1->column("xT1",(double)track.xT1());
	nTupleXY1->column("xT2",(double)track.xT2());
	nTupleXY1->column("xT3",(double)track.xT3());

	nTupleXY1->column("ax",(double)track.ax());
	nTupleXY1->column("bx",(double)track.bx());
	nTupleXY1->column("cx",(double)track.cx());
	nTupleXY1->column("ay",(double)track.ay());
	nTupleXY1->column("by",(double)track.by());
	nTupleXY1->column("dRatio",(double)track.dRatio());

	nTupleXY1->write();

	if(assoc)
	  {for(auto it = (track).hits().begin();(track).hits().end()!=it; ++it){
	      if( (int)part != (int) (*it)->zone()%2) continue;
	      (*it)->setUsed();}}
      }


#else //if the ntuples in case 1 are already generated and the weights added, the algo is ready to use mva classification with 2 different cuts ( case 0 and case 1) 

      //WE'RE IN CASE == 1
      float MLP1_scut=0.5; //based on s/sqrt(s+b)                                                                                                                                                        
      float MLP1_bcut=0.2;                                                                                                                                                                                 

      for(auto track_ptr = m_trackCandidates[(int)part].begin(); m_trackCandidates[(int)part].end()!=track_ptr ; ++track_ptr){


	if((*track_ptr).Case()!=icase) continue;
	if(!(*track_ptr).valid()) continue;


	//MLPReader_input[0]=(double)icase;                                                                                                                                                              
	MLPReader_input[0]=(double)(*track_ptr).chi2();
	MLPReader_input[1]=(double)(*track_ptr).nx();
	MLPReader_input[2]=(double)(*track_ptr).ny();
	MLPReader_input[3]=(double)abs((*track_ptr).X0());
	MLPReader_input[4]=(double)(*track_ptr).worstHitChi2();



	float quality = MLPReader1->GetMvaValue(MLPReader_input);

	if(quality>MLP1_scut)
	  {
	    for(auto it = (*track_ptr).hits().begin();(*track_ptr).hits().end()!=it; ++it){
	      if( (int)part != (int) (*it)->zone()%2) continue;
	      (*it)->setUsed();}
	  }
	else
	  if(quality<MLP1_bcut)
	    (*track_ptr).setValid(false);
      }
#endif
      //ntuples generation for case1
    }


#else //else MVA_classification_inFlagHits                                                                                                                                                                

  std::sort( m_trackCandidates[(int)part].begin() , m_trackCandidates[(int)part].end() , PrHybridSeedTrack::LowerBySize()); //bigger size is in front                                                     
  for(auto track = m_trackCandidates[(int)part].begin(); m_trackCandidates[(int)part].end()!=track ; ++track){
    // important the sorting for this break                                                                                                                                                               
    if((*track).size() < m_SizeFlag[icase]) break;
    if(!(*track).valid()) continue;
    if((*track).Case()!=icase) continue;
    // Important the sorting of before for this break                                                                                                                                                     
    if(!( ((*track).size()==11
           && (*track).chi2PerDoF()< m_MaxChi2Flag[icase]
           && std::fabs((*track).X0())<m_MaxX0Flag[icase])
          || ((*track).size()==12 )  )
       ) continue;
    for(auto it = (*track).hits().begin();(*track).hits().end()!=it; ++it){
      if( (int)part != (int) (*it)->zone()%2) continue;
      (*it)->setUsed();
    }
  }
#endif //endif MVA_classification_inFlagHits                                                                                                                                                              

}


//=========================================================================
//  Convert to LHCb tracks
//=========================================================================
void PrHybridSeeding::makeLHCbTracks( LHCb::Tracks* result, unsigned int part){
  for( auto track = m_trackCandidates[(int)part].begin();m_trackCandidates[(int)part].end() != track; ++track ){
    if ( !(*track).valid() ) continue;
    if ( msgLevel(MSG::DEBUG) ) debug()<<"Creating LHCb Track"<<endmsg;
    //***** EXACT SAME IMPLEMENTATION OF PatSeedingTool *****//
    LHCb::Track* out = new LHCb::Track;
    out->setType( LHCb::Track::Ttrack );
    out->setHistory( LHCb::Track::PrSeeding );
    out->setPatRecStatus( LHCb::Track::PatRecIDs);
    for( auto itH = (*track).hits().begin(); (*track).hits().end()!=itH; ++itH){
      out->addToLhcbIDs( (*itH)->id());
    }
    LHCb::State temp(
                     Gaudi::TrackVector( (*track).x( m_zReference), (*track).y( m_zReference),
                                         (*track).xSlope( m_zReference), (*track).ySlope(), 0.),
                     m_zReference, LHCb::State::AtT);
    double qOverP , sigmaQOverP;
    float m_momentumScale = 35.31328;
    //float m_momentumScale = 40.3751; (PatSeeding is also using this value)
    const float scaleFactor = m_magFieldSvc->signedRelativeCurrent();
    if( m_momentumTool->calculate( &temp, qOverP, sigmaQOverP, true).isFailure()){
      if( std::abs( scaleFactor) <1.e-4){
        qOverP = (((*track).cx() < 0) ? -1. : 1.) *
          ((scaleFactor < 0.) ? -1. : 1.) / Gaudi::Units::GeV;
        sigmaQOverP = 1. / Gaudi::Units::MeV;
      }
      else{
        qOverP = (*track).cx() * m_momentumScale / (-1*scaleFactor);
        sigmaQOverP = 0.5 * qOverP;
      }
    }
    temp.setQOverP( qOverP);
    Gaudi::TrackSymMatrix& cov = temp.covariance();
    cov(0,0) = m_stateErrorX2;
    cov(1,1) = m_stateErrorY2;
    cov(2,2) = m_stateErrorTX2;
    cov(3,3) = m_stateErrorTY2;
    cov(4,4) = sigmaQOverP * sigmaQOverP;
    for( const float z: m_zOutputs ){
      temp.setX((*track).x(z));
      temp.setY((*track).y(z));
      temp.setZ(z);
      temp.setTx((*track).xSlope(z));
      temp.setTy((*track).ySlope());
      out->addToStates( temp );
    }
    out->setChi2AndDoF( (*track).chi2(), (unsigned)(*track).nDoF());
    result->insert( out);
  }
}


void PrHybridSeeding::solveParabola(const PrHit* hit1,const PrHit* hit2,const PrHit* hit3,float& a1, float& b1,float& c1){
  const float z1 = hit1->z() - m_zReference;//
  const float z2 = hit2->z() - m_zReference;//
  const float z3 = hit3->z() - m_zReference;//
  const float x1 = hit1->x();
  const float x2 = hit2->x();
  const float x3 = hit3->x();
  //const float e = m_dRatio;
  const float corrZ1 = 1.+m_dRatio*z1;
  const float corrZ2 = 1.+m_dRatio*z2;
  const float corrZ3 = 1.+m_dRatio*z3;
  const float det = (z1*z1)*corrZ1*z2 + z1*(z3*z3)*corrZ3 + (z2*z2)*corrZ2*z3 - z2*(z3*z3)*corrZ3 - z1*(z2*z2)*corrZ2 - z3*(z1*z1)*corrZ1;
  const float recdet = 1./det;
  if( std::fabs(det) < 1e-8 ){
    b1 = (x3-x1)/(z3-z1);
    c1 = x1 - b1*z1;
    a1 = 0.f;
    return;
  }
  const float det1 = (x1)*z2 + z1*(x3) + (x2)*z3 - z2*(x3) - z1*(x2) - z3*(x1);
  const float det2 = (z1*z1)*corrZ1*x2 + x1*(z3*z3)*corrZ3 + (z2*z2)*corrZ2*x3 - x2*(z3*z3)*corrZ3 - x1*(z2*z2)*corrZ2 - x3*(z1*z1)*corrZ1;
  const float det3 = (z1*z1)*corrZ1*z2*x3 + z1*(z3*z3)*corrZ3*x2 +(z2*z2)*corrZ2*z3*x1 - z2*(z3*z3)*corrZ3*x1 - z1*(z2*z2)*corrZ2*x3 - z3*(z1*z1)*corrZ1*x2;
  a1 = recdet*det1;
  b1 = recdet*det2;
  c1 = recdet*det3;
}

//=========================================================================
//  Fit the track, return OK if fit sucecssfull
//=========================================================================
bool PrHybridSeeding::fitSimultaneouslyXY( PrHybridSeedTrack& track , unsigned int iCase){
  if( track.nx() <4 || track.ny() <4 ) return false;
  float mat[15];
  float rhs[5];
  // unsigned int nHitsX = 0;
  // unsigned int nHitsStereo = 0;
  const float zRef = m_zReference;
  for (unsigned int loop = 0; 3 > loop ; ++loop ){
    if(loop ==1 && m_useCorrPos){
      float RadiusPosition = std::sqrt((0.0005*0.0005*track.ax()*track.ax() + 0.001*0.001*track.y(zRef)*track.y(zRef)));
      float dRatioPos = -1.*(m_dRatioPar[0] + m_dRatioPar[1]*RadiusPosition + m_dRatioPar[2]*RadiusPosition*RadiusPosition);
      track.setdRatio(dRatioPos);
    }
    std::fill(mat,mat+15,0.);
    std::fill(rhs,rhs+5,0.);
    for( auto itH = track.hits().begin(); track.hits().end() != itH; ++itH ){
      const float w = (*itH)->w();
      const float dxdy = (*itH)->dxDy();
      const float yOnTrack = track.yOnTrack( (*itH) ) ;
      const float   dz = 0.001*((*itH)->z( yOnTrack ) - zRef );
      const float dRatio = track.dRatio();
      const float deta = dz*dz*(1. + dz*dRatio);
      const float wdz = w * dz;
      const float weta = w * deta;
      const float wdxdy = w * dxdy;
      const float wdxdydz = wdxdy * dz;
      const float dist = track.distance(*itH);
      mat[0] += w;
      mat[1] += wdz;  mat[2] += wdz * dz;
      mat[3] += weta; mat[4] += weta * dz; mat[5] += weta * deta;
      mat[6] -= wdxdy;mat[7] -= wdxdydz;   mat[8] -= wdxdy * deta;  mat[9] += wdxdy * dxdy;
      mat[10] -= wdxdydz; mat[11] -= wdxdydz * dz;mat[12] -= wdxdydz * deta;mat[13] += wdxdydz * dxdy; mat[14] += wdxdydz * dz * dxdy;
      // fill right hand side
      rhs[0] += w * dist;
      rhs[1] += wdz * dist;
      rhs[2] += weta * dist;
      rhs[3] -= wdxdy * dist;
      rhs[4] -= wdxdydz * dist;
    }//Loop over Hits to fill the matrix
    // decompose matrix, protect against numerical trouble
    //track.setnXnY( nHitsX, nHitsStereo );
    //if(nHitsX < 4 || nHitsStereo < 4) return false;
    ROOT::Math::CholeskyDecomp<float, 5> decomp(mat);
    if (!decomp) return false;
    decomp.Solve(rhs);
    rhs[1]*=1.e-3;
    rhs[2]*=1.e-6;
    rhs[4]*=1.e-3;
    rhs[3]-=rhs[4]*zRef;
    //crazy values!
    if( (std::fabs(rhs[0]) > 1e4 || std::fabs(rhs[1]) > 5. ||
         std::fabs(rhs[2]) > 1e-3 || std::fabs(rhs[3]) > 1e4 || std::fabs(rhs[4]) > 1.)) return false;
    track.updateParameters(rhs[0],rhs[1],rhs[2],rhs[3],rhs[4]);
  }
  float maxChi2 = std::numeric_limits<float>::min();
  float chi2_onHit = std::numeric_limits<float>::max();
  for( auto itH = track.hits().begin(); track.hits().end() != itH; ++itH ){
    chi2_onHit = track.chi2( *itH);
    if ( chi2_onHit > maxChi2 ){
      maxChi2 = chi2_onHit;
    }
  }//Set Max Chi2DoF
  float X0 = track.ax() - track.bx()*m_zReference+track.cx()*m_ConstC;
  track.setX0(X0);
  if( ( track.recovered()) ){
    if( maxChi2 < m_recoChiOutlier) return true;
    return false;
  }else{
    if( (track.size()>10) && maxChi2<m_maxChi2HitFullFitHigh[iCase]) return true;
    if( maxChi2<m_maxChi2HitFullFitLow[iCase] && track.size()<11) return true;
  }
  return false;
}

//=======================================
//Fit Only X Projection
//=======================================
bool PrHybridSeeding::fitXProjection(PrHybridSeedTrack& track, unsigned int iCase ){
  if(msgLevel(MSG::DEBUG)) debug()<<"Fitting"<<endmsg;
  if(track.size()<m_minXPlanes) return false;
  float mat[6];
  float rhs[3];
  const float zRef = m_zReference;
  for(unsigned int loop = 0;3>loop;++loop){
    std::fill(mat,mat+6,0.);
    std::fill(rhs,rhs+3,0.);
    for( auto itH = track.hits().begin(); track.hits().end() != itH; ++itH ){
      const float dRatio = track.dRatio();
      float w = (*itH)->w();//squared
      const float dz= 0.001*((*itH)->z() - zRef);
      const float deta = dz*dz*(1. + dRatio*dz);
      const float dist = track.distance( *itH );
      // if(msgLevel(MSG::DEBUG)) debug()<<"Loop \t"<<loop<<"\n Distance From Hit \t"<<dist<<endmsg;
      mat[0]+= w;
      mat[1]+= w * dz;   mat[2]+= w * dz * dz;
      mat[3]+= w * deta; mat[4]+= w * dz * deta;  mat[5]+= w * deta * deta;
      rhs[0]+= w * dist;
      rhs[1]+= w * dist * dz;
      rhs[2]+= w * dist * deta;
    }
    ROOT::Math::CholeskyDecomp<float,3> decomp(mat);
    if(!decomp){
      return false;
    }
    //Solve linear system
    decomp.Solve(rhs);
    rhs[1]*=1.e-3;
    rhs[2]*=1.e-6;
    if(msgLevel(MSG::DEBUG)){ debug()<<"Loop \t"<<loop<<"\n a = \t"<<rhs[0]<<"\n b = \t"<<rhs[1]<<"\n c = \t"<<rhs[2]<<endmsg;}
    // protect against unreasonable track parameter corrections
    // (check that out)
    if(std::fabs(rhs[0]) > 1.e4 || std::fabs(rhs[1]) > 5. ||
       std::fabs(rhs[2]) > 1.e-3 ) return false;
    //Small corrections
    track.updateParameters(rhs[0],rhs[1],rhs[2],0.,0.);
    if(loop==0){
      track.setdRatio(m_dRatio);
    }
    //Put back later faster maybe
    if(loop >0 && std::abs(rhs[0]) < 5e-5 && std::abs(rhs[1]) < 5e-8 &&
       std::abs(rhs[2]) < 5e-11){
      break;
    }
  }
  //Compute some values on the track
  float maxChi2 = 0.;
  for ( auto itH = track.hits().begin(); track.hits().end() != itH; ++itH )  //Loop over all hits in PrHybridSeedTrack
  {
    float chi2_onHit = track.chi2( *itH );
    if ( chi2_onHit > maxChi2 ){
      maxChi2 = chi2_onHit;
    }
  }
  return (maxChi2 < m_maxChi2HitsX[iCase]) ;
  return false;
}



bool PrHybridSeeding::removeWorstAndRefit(PrHybridSeedTrack& track, unsigned int iCase){
  //maybe useless?
  if(msgLevel(MSG::DEBUG)){ debug()<<"Removing Worst UV layer"<<endmsg;}
  auto worst= track.hits().begin();
  float maxChi2 = std::numeric_limits<float>::min();
  float chi_hit = std::numeric_limits<float>::min();
  if( track.nx() <4 || track.ny()<4) return false;
  for( auto itH = track.hits().begin(); track.hits().end() != itH; ++itH ){
    chi_hit = track.chi2( (*itH));
    if( chi_hit  > maxChi2){
      maxChi2 = chi_hit;
      worst = itH;
    }
  }
  track.hits().erase(worst);
  if( (*worst)->isX()){
    track.setnXnY( track.nx()-1, track.ny());
  }else{
    track.setnXnY( track.nx(), track.ny()-1);
  }
  if(msgLevel(MSG::DEBUG)) debug()<<"returning"<<endmsg;
  return fitSimultaneouslyXY(track, iCase);
    
  return false;
}

//=========================================================================
//  Remove the worst hit and refit.
//=========================================================================
bool PrHybridSeeding::removeWorstAndRefitX ( PrHybridSeedTrack& track , unsigned int iCase)
{
  if(msgLevel(MSG::DEBUG)){ debug()<<"NHits X projection passed = "<<track.size()<<endmsg;}
  if(track.size()<=m_minXPlanes) return false;
  if(msgLevel(MSG::DEBUG)) debug()<<"Removing Worst and Refitting X"<<endmsg;
  float maxChi2 = std::numeric_limits<float>::min();
  auto worst = track.hits().begin();
  for ( auto itH = track.hits().begin(); track.hits().end() != itH; ++itH ){
    float chi2 = track.chi2( *itH );
    if( chi2 > maxChi2 ){
      maxChi2 = chi2;
      worst = itH;
    }
  }
  track.hits().erase(worst);
  return fitXProjection(track, iCase);
}
void PrHybridSeeding::setChi2( PrHybridSeedTrack& track ){
  float chi2 = 0.;
  int   nDoF = -5;
  for (PrHit* hit : track.hits()){
    chi2 +=  track.chi2( hit);
    nDoF += 1;
  }
  //setChi2 is called als with U/V Hits
  track.setChi2( chi2, nDoF );
}

//=========================================================================
//  Set the chi2 of the track
//=========================================================================
void PrHybridSeeding::setChi2X ( PrHybridSeedTrack& track ) {
  float chi2 = 0.;
  int   nDoF = -3;  // Fitted a parabola
  for (PrHit* hit : track.hits()) {
    if(msgLevel(MSG::DEBUG)){  
      if(hit->dxDy() !=0){ debug()<<"You were picking up Stereo Layers!!!"<<endmsg;  } 
    }
    const float Chi2Hit = track.chi2( hit);
    chi2+=  Chi2Hit;
    nDoF += 1;
  }
  if(msgLevel(MSG::DEBUG)){ debug()<<"Chi2 Set for track = \t"<<chi2<<endmsg;}
  track.setChi2( chi2, nDoF );
}




void PrHybridSeeding::findXProjections(unsigned int part, unsigned int iCase){
  PrHits parabolaSeedHits; 
  // List of hits in both x-layers in T2 from a 2-hit combo
  // pushed back in this container
  parabolaSeedHits.reserve(12);
  
  //xHitsLists is the xHits you have
  std::vector<PrHits> xHitsLists;
  xHitsLists.reserve(100);
  PrHits xHits;
  xHits.reserve(10);
  //----------  Identifying layers [start]
  //----------  Just do the 1st one here 1st layer and last one
  unsigned int firstZoneId = -1;
  unsigned int lastZoneId = -1;
  std::array<unsigned int , PrHitManager::Numbers::NFTXLayers> zones;
  std::array<unsigned int , PrHitManager::Numbers::NFTXLayers> planeCode;
  std::array<float        , PrHitManager::Numbers::NFTXLayers> zLays;
  if(0 == iCase){
    firstZoneId = s_T1X1 | part;
    lastZoneId  = s_T3X2 | part;
  }else if ( 1 == iCase ){
    firstZoneId = s_T1X2 | part;
    lastZoneId  = s_T3X1 | part;
  }else if ( 2 == iCase ){
    firstZoneId = s_T1X1 | part;
    lastZoneId  = s_T3X1 | part;
  }else if ( 3 == iCase ){
    firstZoneId = s_T1X2 | part;
    lastZoneId  = s_T3X2 | part;
  }
  if(msgLevel(MSG::DEBUG)){ debug()<<"\t Loading Case Hit in first and last Zone"<<endmsg;}
  
  //Array[0] = first layer in T1 for 2-hit combo
  auto firstZone = GetZone( firstZoneId );
  auto lastZone = GetZone( lastZoneId);
  auto T2_X1 = GetZone( s_T2X1 | part);
  auto T2_X2 = GetZone( s_T2X2 | part);
  
  zones[0]    = firstZoneId;
  zLays[0]    = firstZone.z();
  planeCode[0]= firstZone.planeCode();
  //Array[1] = last  layer in T3 for 2-hit combo
  zones[1]= lastZoneId ;
  zLays[1]= lastZone.z();
  planeCode[1]= lastZone.planeCode();
  if(msgLevel(MSG::DEBUG)){ debug()<<"Hits in last and first Zone Loaded"<<endmsg;}
  //Array[2] = T2-1st x-layers
  zones[2]=     s_T2X1 | part;
  zLays[2] =    T2_X1.z();
  planeCode[2]= T2_X1.planeCode();
  //Array[3] = T2-2nd x-layers
  zones[3] =           s_T2X2|part;
  zLays[3] =   T2_X2.z();
  planeCode[3]= T2_X2.planeCode();
  if(msgLevel(MSG::DEBUG)){ debug()<<"Hits in last and first Zone Loaded"<<endmsg;}
  unsigned int i=4;
  //Add extra layer HERE, if needed!!!!
  for(unsigned int xZoneId : {s_T1X1, s_T1X2, s_T3X1, s_T3X2}){
    xZoneId =xZoneId | part;
    if(xZoneId != firstZoneId && xZoneId != lastZoneId){
      auto Zone = GetZone(xZoneId);
      zones[i] = xZoneId; zLays[i] = Zone.z(); planeCode[i] = Zone.planeCode();
      i++;
    }
  }
  float Zlast  = zLays[1];//zones[1]->z();
  float Zfirst = zLays[0];//zones[0]->z();
  float invZfirst = 1./Zfirst;// zones[0]->z();
  float invZlZf = 1./( Zlast-Zfirst );
  float slope = (m_tolAtX0Cut[iCase]-m_TolX0SameSign[iCase])/(m_x0Cut[iCase]-m_x0SlopeChange[iCase]);
  float slopeopp =  (m_tolAtx0CutOppSig[iCase] -m_tolX0Oppsig[iCase])/(m_x0Cut[iCase]-m_x0SlopeChange2[iCase]);
  Boundaries Bounds;
  //============Initialization
  std::array< HitIterat, PrHitManager::Numbers::NFTXLayers> endsZone;
  std::array< HitIterat, PrHitManager::Numbers::NFTXLayers> begZone;
  for(unsigned int i=0;i<PrHitManager::Numbers::NFTXLayers;i++){
    //Always: first=0 last=1 middles=2,3 remains=4,5 (order of loops)
    InitializeBB(Bounds[i],zones[i]);
    begZone[i]  = m_hitManager->getIterator_Begin( zones[i]);
    endsZone[i] = m_hitManager->getIterator_End( zones[i]);
  }
  Bounds[1].end = Bounds[1].begin;//initializing both of iterators
  PrHit *fHit;
  std::array<bool,PrHitManager::Numbers::NFTXLayers> first;  
  for( unsigned int j = 0; j<PrHitManager::Numbers::NFTXLayers; j++){ first[j] = true; }
  //Used to cache the position to do a "look-around search"
  std::array<float,PrHitManager::Numbers::NFTXLayers> xminPrev, xmaxPrev; 
  for(; Bounds[0].begin != Bounds[0].end; ++Bounds[0].begin) //for a hit in first layer
  {
    fHit = &(*(Bounds[0].begin));
    if( fHit->isUsed() && m_removeFlagged) continue;
    float tx_inf = fHit->x() * invZfirst;
    float xProjeInf = tx_inf *Zlast ;
    float tolHp = m_TolFirstLast[iCase];
    float maxXl = xProjeInf + tx_inf*m_alphaCorrection[iCase] +tolHp;
    float minXl = xProjeInf + tx_inf*m_alphaCorrection[iCase] -tolHp;
    Bounds[1].end    = m_hitManager->get_upperBound_log( Bounds[1].end, endsZone[1] , maxXl);
    Bounds[1].begin  = m_hitManager->get_lowerBound_lin( Bounds[1].begin, Bounds[1].end, minXl);
    xminPrev[1] = minXl;
    xmaxPrev[1] = maxXl;
    
    PrHit* lHit;
    float tx_pickedcombination;
    float x0;
    float CorrX0;
    float x0new;
    float xMax;
    float xMin;
    float max = 0.;
    float min = 0.;
    float xProjected;
    float xProjectedCorrected;
    unsigned int nfired = 1;
    for(auto Lhit=Bounds[1].begin; Lhit != Bounds[1].end; ++Lhit){
      lHit = &(*Lhit);
      if(lHit->isUsed() && m_removeFlagged) continue;
      //Define the boundaries in T2 to look for parabolaSeedHits 
      //( push them back now or later? ) (dependent on x(z=0) )
      //from the 2-hit combo
      tx_pickedcombination = (lHit->x() - fHit->x())*invZlZf;
      x0 = fHit->x() - tx_pickedcombination*Zfirst;
      CorrX0 = m_x0Corr[iCase]*x0;
      x0new = x0*(1.+m_x0Corr[iCase]);
      if(x0>0.){
        min = x0 > m_x0SlopeChange[iCase]?  -slope*( x0 - m_x0SlopeChange[iCase]) - m_TolX0SameSign[iCase] : -m_TolX0SameSign[iCase];
        max = x0 > m_x0SlopeChange2[iCase]? slopeopp*( x0 - m_x0SlopeChange2[iCase]) + m_tolX0Oppsig[iCase] : +m_tolX0Oppsig[iCase];
      }else{
        max = x0 < - m_x0SlopeChange[iCase]? -slope*( x0 + m_x0SlopeChange[iCase]) + m_TolX0SameSign[iCase]: m_TolX0SameSign[iCase];
        min = x0 < - m_x0SlopeChange2[iCase]? slopeopp*( x0 + m_x0SlopeChange2[iCase]) - m_tolX0Oppsig[iCase]: -m_tolX0Oppsig[iCase] ;
      }
      parabolaSeedHits.clear();
      nfired = 2;
      for( int i = 2; i<4; ++i){
        xProjected = x0 + zLays[i]*tx_pickedcombination ;
        xProjectedCorrected = xProjected + CorrX0 ;
        xMin =   xProjectedCorrected + min ;//may add a 1.0 mm here?
        xMax =   xProjectedCorrected + max ;//may add a 1.0 mm here?
        if( xMin > xMax){
          if(msgLevel(MSG::DEBUG)){  debug()<<" xMin is bigger than xMax in T2 layers : "<<i<<endmsg;}
          std::swap(xMin, xMax);
        }
        if( first[i] ){
          Bounds[i].begin = m_hitManager->get_lowerBound_log( Bounds[i].begin, Bounds[i].end, xMin );
          Bounds[i].end = m_hitManager->get_upperBound_lin( Bounds[i].begin, Bounds[i].end, xMax );
          first[i] = false;
        }else{
          LookAroundMax( Bounds[i], xmaxPrev[i], xMax, endsZone[i], begZone[i] );
          LookAroundMin( Bounds[i], xminPrev[i], xMin, Bounds[i].end, begZone[i] );
        }
        xminPrev[i] = xMin;
        xmaxPrev[i] = xMax;
        bool fired = false;
        if( Bounds[i].end - Bounds[i].begin > 0){
          for( auto itM = Bounds[i].begin ; itM!= Bounds[i].end; ++itM){
            if( (*itM).isUsed() && m_removeFlagged) continue;
            if( msgLevel(MSG::DEBUG)){
              debug()<<"push back parabola seed hit"<<endmsg;
            }
            parabolaSeedHits.push_back( &(*itM));
            fired = true;
          }
        }
        nfired += (fired);
      }
      if( nfired < 3 ){
        if(msgLevel(MSG::DEBUG)){
          debug()<<"Not enough fired layers : N fired = "<<nfired<<endmsg;
        }
        continue;
      }
      //=======================================================
      // We have 1 Hit in 1st 1 Hit in last and a
      // vector of Hits for in-between
      //=======================================================
      xHitsLists.clear();
      //=======================================================
      //Sort the ParabolaSeedHits for in-between layers by increasing distance from the Projected Corrected position only when we have more than 1 ParabolaSeedHit
      //=======================================================
      if(parabolaSeedHits.size()>1){
        //Principle of the Lambda funtion, Hits sorted wrt distance from linear Projection 1st-3rd layer
        std::sort( parabolaSeedHits.begin(),parabolaSeedHits.end(),
                   [x0new,tx_pickedcombination](const PrHit* lhs, const PrHit* rhs)
                   ->bool{
                     return std::fabs( lhs->x() - ( x0new + lhs->z()*tx_pickedcombination)) < std::fabs( rhs->x() - (x0new + rhs->z()*tx_pickedcombination));} );
      }
      if(msgLevel(MSG::DEBUG)){ info()<<"The Lambda Function Sorting end"<<endmsg;}
      unsigned int maxParabolaSeedHits = m_maxParabolaSeedHits;
      if(parabolaSeedHits.size()<m_maxParabolaSeedHits){
        maxParabolaSeedHits = parabolaSeedHits.size();
      }else{
        maxParabolaSeedHits = m_maxParabolaSeedHits;
      }
      for(unsigned int i = 0; i<maxParabolaSeedHits;++i){ //for a give hit in 1st layer, and a given hit in last layer
        float a = 0;
        float b = 0;
        float c = 0;
        solveParabola(fHit,parabolaSeedHits[i],lHit,a,b,c); //Extrapolation with dRatio
        //===================================================
        // Look in all the other layers except the
        // 1st/last/zone except the parabolaSeedHit
        //===================================================
        //Loop on all the xZones
        float dz ;
        float xAtZ;
        float xMaxAtZ;
        float xMinAtZ;
        PrHit *bestProj;
        float bestDist;
        HitIterat itH;
        HitIterat itE;
        xHits.clear();
        for(int j =2 ;j< PrHitManager::Numbers::NFTXLayers;++j){ //look in zones except the 1st and last zones
          if(msgLevel(MSG::DEBUG)){debug()<<"Selecting ParSeedHits"<<endmsg;}
          if((int) planeCode[j] == (int)parabolaSeedHits[i]->planeCode()) continue;
          dz = zLays[j] - m_zReference;
          xAtZ= a * dz * dz * (1. + m_dRatio* dz) + b * dz + c; //Prediction of the position with the cubic correction!
          xMaxAtZ = xAtZ + m_tolRemaining[iCase]; // * (1. + dz*invZref ); Make them z-dependent and/or x0 dependent to cure for track-model error?
          xMinAtZ = xAtZ - m_tolRemaining[iCase]; // * (1. + dz*invZref ); 
          // Parabola here is larger than tolerances in x0! ( 0.75 mm up / down when x0 ~0 ) 
          // may we want ot make the m_tolRemaining[iCase] ( x0 dependent too ?? )
          // here you are re-using the tolerance defined before from the linear extrapolation of 2-hit combo
          itH = m_hitManager->get_lowerBound_log(Bounds[j].begin, Bounds[j].end, xMinAtZ);
          itE = m_hitManager->get_upperBound_lin( itH, Bounds[j].end , xMaxAtZ);
          //Linear search of them. 
          //auto itE = LookAroundMax_it( Bounds[j].end  , xmaxPrev[j], xMaxAtZ, endsZone[j], begZone[j]);
          //auto itE = LookAroundMax_it( Bounds[j].end, xmaxPrev[j], xMaxAtZ, endsZone[j], begZone[j]);
          bestProj = nullptr;
          bestDist = 10.0; // 1 cm around predicted position
          PrHit* hit;
          for(; itH != itE ; ++itH){
            hit = &*itH;
            if(hit->isUsed() && m_removeFlagged) continue;
            if(std::fabs(hit->x() - xAtZ)  <  bestDist){
              bestDist = std::fabs(hit->x() - xAtZ);
              bestProj = hit;
            }
          }
          if(bestProj != nullptr){
            xHits.push_back( bestProj);
          }
        }//end loop xZones excluding parabolaseed hit
        //in xHits are not present the first layer and last + parabola seed hits
        if(msgLevel(MSG::DEBUG)){ debug()<<"End Loop in between zones to pick up Projection of parabola"<<endmsg;}
        //Up to here xHits contains only the remaining layers hits!
        if( xHits.size() < 2 ) continue; //next parabolaSeedHits ; you must find 2 hits in remaining layers for this given parabola seed hits!
        xHits.push_back(fHit);
        xHits.push_back(parabolaSeedHits[i]);
        xHits.push_back(lHit);
        //sort the hits in xHits by LHCbID ( important for clone removal purpose )!
        std::sort(xHits.begin(), xHits.end(), compLHCbID());
        xHitsLists.push_back( xHits);
      }//End loop parabolaSeedHits
      if(msgLevel(MSG::DEBUG)){ 
        debug()<<"End Loop For pick up Parabola Hits and build the xHitsLists"<<endmsg;
        //-------- Remove Duplicates from search in parabolaSeedHits
        debug()<<"xHitsLists size before removing duplicates: "<<xHitsLists.size()<<endmsg;
      }
      if(xHitsLists.size() == 0){
        continue;
      }
      if(xHitsLists.size() > 1){
        //---Remove Duplicates in the HitsList
        if(msgLevel(MSG::DEBUG)){debug()<<"Remove duplicate xHitsLists"<<endmsg;}
        std::stable_sort( xHitsLists.begin(), xHitsLists.end() );
        xHitsLists.erase( std::unique(xHitsLists.begin(), xHitsLists.end()), xHitsLists.end());
      }
      if(msgLevel(MSG::DEBUG)){debug()<<"xHitsLists size after removing duplicates: "<<xHitsLists.size()<<endmsg;}
      //----- You have list of xHits : fit them!
      for(PrHits& xHits : xHitsLists){
        if(msgLevel(MSG::DEBUG)){ debug()<<"Fit Track"<<endmsg;}
        //------    Create the track
        PrHybridSeedTrack temp_track( part , m_zReference , xHits);
        //------    Set the dRatio value for the track here
        temp_track.setdRatio(m_dRatio);
        

#ifdef ntuplesGeneration_inX

        if(UNLIKELY(m_debug) ){   debug()<<"Attempting to Fit the following Track"<<endmsg; printTrack(temp_track);  }
        //-----     OK is the status of the fit                                                                                                                                                           

        bool OK = false;
        if(temp_track.size()>m_minXPlanes){
          //-----   no N = m_minXPlanes hits at first fit                                                                                                                                                  

          OK = fitXProjection(temp_track,iCase);
        }
        while(!OK && temp_track.size()>m_minXPlanes){
          if(temp_track.size() <=m_minXPlanes){
	    OK = false;
            break;
          }
          nIter++;
          if( nIter==1 && temp_track.size() == 5){
            //fit for entering tracks with 5 hits was failing , do not remove hit and refit it!                                                                                                           

            OK = false;
            break;
          }
          if( temp_track.size() > m_minXPlanes){
            OK = removeWorstAndRefitX(temp_track,iCase);
          }
	}

	LHCb::MCParticle* partic = nullptr;
        double efficiency = -1;
        int nAssocHits = -10;
        if(OK){
	  bool assoc = AssocTrack(temp_track,efficiency,partic,nAssocHits);


	  Tuple nTupleXZ = nTuple("XZtracks","XZtracks");
	  PrPlaneHybridCounter counter;
	  counter.set(temp_track.hits().begin(),temp_track.hits().end());


	  nTupleXZ->column("assoc",assoc);
	  nTupleXZ->column("Case",(double)iCase);
	  nTupleXZ->column("chi2",(double)temp_track.chi2());
	  nTupleXZ->column("nx",(double)temp_track.nx());
	  nTupleXZ->column("ny",(double)temp_track.ny());
	  nTupleXZ->column("X0Back",(double)abs(temp_track.X0()));
	  nTupleXZ->column("worstHitChi2",(double)temp_track.worstHitChi2());
	  nTupleXZ->column("nXHits",(double)temp_track.hits().size());

	  nTupleXZ->column("ax",(double)temp_track.ax());
	  nTupleXZ->column("bx",(double)temp_track.bx());
	  nTupleXZ->column("cx",(double)temp_track.cx());
	  nTupleXZ->column("ay",(double)temp_track.ay());
	  nTupleXZ->column("by",(double)temp_track.by());
	  nTupleXZ->column("dRatio",(double)temp_track.dRatio());



	  nTupleXZ->column("nbSingleX",(double)counter.nbSingleX());
	  nTupleXZ->column("nXHits",(int)temp_track.hits().size());
	  //nTupleXZ->column("ax",temp_track.ax());                                                                                                                                                       
	  //nTupleXZ->column("bx",temp_track.bx());                                                                                                                                                       
	  //nTupleXZ->column("cx",temp_track.cx());                                                                                                                                                       


	  nTupleXZ->write();
	  //      setChi2X(temp_track);                                                                                                                                                                  
	  temp_track.setCase( iCase );    // useful later in clone killing step                                                                                                                           

	  temp_track.setXT2(temp_track.x( (m_zones[12]->z())));
	  temp_track.setXT3(temp_track.x( (m_zones[23]->z())));
	  temp_track.updateIDs();  //for bloom filter                                                                                                                                                     

	  nTupleXZ->column("xT1",(double)temp_track.xT1());
	  nTupleXZ->column("xT2",(double)temp_track.xT2());
	  nTupleXZ->column("xT3",(double)temp_track.xT3());

	  if(assoc)
	    m_xCandidates[(int)part].push_back(temp_track); //The X Candidate is created                                                                                                                
        }
#else






	//------    Algorithm allows to go down to 4 hits only from 6 hits track refit by construction. 
        int nIter = 0;
        if(msgLevel(MSG::DEBUG) ){   debug()<<"Attempting to Fit the following Track"<<endmsg; printTrack(temp_track);  }
        //-----     OK is the status of the fit
        bool OK = false;
        if(temp_track.size()>m_minXPlanes){ 
          //-----   no N = m_minXPlanes hits at first fit
          OK = fitXProjection(temp_track,iCase);
        }
        while(!OK && temp_track.size()>m_minXPlanes){
          if(temp_track.size() <=m_minXPlanes){
            OK = false;
            break;
          }
          nIter++;
          if( nIter==1 && temp_track.size() == 5){ 
            //fit for entering tracks with 5 hits was failing , do not remove hit and refit it!
            OK = false;
            break;
          }
          if( temp_track.size() > m_minXPlanes){
            OK = removeWorstAndRefitX(temp_track,iCase);
          }
        }
        //---- Fit has to be fine to make an x-candidate
        if( OK ){
          setChi2X(temp_track);
        }
        if( OK &&
            temp_track.size() >= m_minXPlanes &&
            (( temp_track.chi2PerDoF() < m_maxChi2DoFX[iCase]))){
          temp_track.setCase( iCase );
          temp_track.setXT1( temp_track.x( GetZone(0).z()  ));
          temp_track.setXT2( temp_track.x( GetZone(12).z() ));
          temp_track.setXT3( temp_track.x( GetZone(23).z() ));
          temp_track.updateIDs();//for bloom filter
          m_xCandidates[(int)part].push_back(temp_track); //The X Candidate is created
        }
#endif
      }//end Loop xHist:xHitsLists
    }//end loop Last Zone given a firsZone selected
  }//end loop first zone
}

//Update the bounds (Hiterator Pairs) and old Min => new Min searching linearly around current boundary begin
inline void PrHybridSeeding::LookAroundMin( IterPairs& bound, float& oMin,float& nMin,  HitIterat& eZone ,HitIterat& bZone){
  if( nMin < oMin){ bound.begin  = m_hitManager->get_lowerBound_lin_reverse( bZone, bound.begin, nMin);
  }else{ bound.begin  = m_hitManager->get_lowerBound_lin( bound.begin, eZone, nMin); }
  oMin = nMin;
}

//Just grab the iterator searching around current bound begin
inline HitIterat PrHybridSeeding::LookAroundMin_it( HitIterat& currBeg, float& oMin,float& nMin,  HitIterat& eZone ,HitIterat& bZone){
  if( nMin < oMin){   return m_hitManager->get_lowerBound_lin_reverse( bZone, currBeg, nMin);
  }else{              return m_hitManager->get_lowerBound_lin( currBeg, eZone, nMin);  }
  // oMin = nMin
}

//Update the bounds (HitIterator Pairs) and old Max => new Max searching linearly around current boundary end
inline void PrHybridSeeding::LookAroundMax( IterPairs& bound, float& oMax,float& nMax,  HitIterat& eZone ,HitIterat& bZone){
  if( nMax > oMax){   bound.end = m_hitManager->get_upperBound_lin( bound.end, eZone, nMax);
  }else{              bound.end = m_hitManager->get_upperBound_lin_reverse(  bZone , bound.end, nMax);  }
  oMax = nMax;
}

//Just grab the iterator searching around current bound end
inline HitIterat PrHybridSeeding::LookAroundMax_it( HitIterat& currEnd, float& oMax,float& nMax,  HitIterat& eZone ,HitIterat& bZone){
  if( nMax > oMax){  return m_hitManager->get_upperBound_lin( currEnd, eZone, nMax); 
  }else{             return m_hitManager->get_upperBound_lin_reverse(  bZone , currEnd, nMax); }
  // oMax = nMax;
}



inline bool PrHybridSeeding::areClones( PrHybridSeedTrack& tr1, PrHybridSeedTrack& tr2, unsigned int& maxCommon){
  unsigned int nCommon = maxCommon;
  auto itH1 = tr1.hits().begin();
  auto itH2 = tr2.hits().begin();
  auto itEnd1 = tr1.hits().end();
  auto itEnd2 = tr2.hits().end();
  while( nCommon!=0 && itH1!=itEnd1 && itH2!=itEnd2){
    if((*itH1)->id() == (*itH2)->id()){
      --nCommon;
      ++itH1;
      ++itH2;
    }
    else if( (*itH1)->id() < (*itH2)->id() ){
      ++itH1;
    }
    else{
      ++itH2;
    }
  }
  return(nCommon==0);
}

inline bool PrHybridSeeding::CloseEnough( PrHybridSeedTrack& tr1, PrHybridSeedTrack& tr2, float distance){
  if(  std::fabs( tr1.xT1() - tr2.xT1() ) < abs(distance) ||
       std::fabs( tr1.xT2() - tr2.xT2() ) < abs(distance) ||
       std::fabs( tr1.xT3() - tr2.xT3() ) < abs(distance) )  return true;
  return false;
}

void PrHybridSeeding::ExtendCluster(PrHits::iterator& itBeg, PrHits::iterator& itEnd , unsigned int iCase,PrHits&myStereo , int &nLay, bool reco){
  if(msgLevel(MSG::DEBUG)){info()<<"Pop counter"<<endmsg;}
  int planes = 0;
  for( auto Hit = itBeg; itEnd!=Hit; ++Hit){
    planes|= 1<<(*Hit)->planeCode();
  }
  nLay = __builtin_popcount( planes);
  float tolTy = m_tolTyOffset[iCase] + (*itBeg)->coord()*m_tolTySlope[iCase] ;
  if(reco) tolTy = m_recoTolTy;
  while( itEnd != myStereo.end() && ((*itEnd)->coord() - (*itBeg)->coord() ) < tolTy && nLay<6){
    //as soon as you find a sixth element exit the loop (in other case extend up to tolerance
    if(UNLIKELY( msgLevel(MSG::DEBUG))){ info()<<" Counting layers "<<endmsg;}
    planes|= 1<<(*itEnd)->planeCode();
    ++itEnd;
    nLay = __builtin_popcount( planes);
  }
}

//=========================================================================
//  Print the whole track
//=========================================================================
void PrHybridSeeding::printTrack ( PrHybridSeedTrack& track ) {
  for ( auto itH = track.hits().begin(); track.hits().end() != itH; ++itH ) {
    info() << format( "dist %7.3f dy %7.2f chi2 %7.2f ", track.distance( *itH ), track.deltaY( *itH ), track.chi2( *itH ) );
    printHit( *itH );
  }
}
void PrHybridSeeding::InitializeBB(IterPairs &BB,unsigned int ZoneNumber){
  BB.begin = m_hitManager->getIterator_Begin(ZoneNumber);
  BB.end   = m_hitManager->getIterator_End(ZoneNumber);
}

//=========================================================================
//  Print the information of the selected hit
//=========================================================================
void PrHybridSeeding::printHit( const PrHit* hit, std::string title ) {
  info() << "  " << title << " "
         << format( " Plane%3d zone%2d z0 %8.2f x0 %8.2f coord %8.3f used%2d ",
                    hit->planeCode(), hit->zone(), hit->z(), hit->x(),
                    hit->coord(), hit->isUsed() );
  info() << endmsg;
}
void PrHybridSeeding::printHit( const PrHit& hit, std::string title){
  printHit(&hit,title);
}



#ifdef MVA_used
//Match the hit with the key of the MCParticle                                                                                                                                                             

bool PrHybridSeeding::matchKey( PrHit* hit, int key){
  //Given the key of the MCParticle, check if the Hit has some relation to it.                                                                                                                             
  //Basically: GetMCParticles( PrHit -> FTCluster associated -> ChannelID of the cluster ) -> key () == key*(of the particle searched for)                                              

LinkedTo<LHCb::MCParticle> fLink( evtSvc() , msgSvc(), LHCb::FTClusterLocation::Default );
//the LinkedTo is linking LHCb::MCParticle to FTClusters which are identified by FTChannelID, and each hit has the FTChannelID                                                                            
  LHCb::FTChannelID id0 = hit->id().ftID();    // id0 is the ID of the current Hit                                                                                                                         
  LHCb::MCParticle* partic = fLink.first(id0); // fLink.first( id0) returns the first MCParticle which has produced the cluster                                                                            
  while(0!=partic){                            // we enter a loop iterating over the clusters, until the key is equal to the key of the MCParticle you extract from the Hit                                
    if(key == partic->key() ) return true;     // int key here is the unique identifier of the MCParticle                                                                                                  
    partic=fLink.next();
  }
  return false;
}


bool PrHybridSeeding::AssocTrack(PrHybridSeedTrack track, double& efficiency,LHCb::MCParticle*& particle, int& nHits){
  std::map<int,int> Counter; //map key and number of associated	
  std::map<LHCb::MCParticle*,int> CounterParticle;
  Counter.clear();
  for(PrHits::iterator hit = track.hits().begin();hit!=track.hits().end();++hit){
    LHCb::FTChannelID id = (*hit)->id().ftID();
    LinkedTo<LHCb::MCParticle> ftLink1(evtSvc(),msgSvc(),LHCb::FTClusterLocation::Default);
    LHCb::MCParticle* part = ftLink1.first(id);
    while(0!=part){
      int key = part->key();
      Counter[key]++;
      CounterParticle[part]++;
      part=ftLink1.next();
    }
  }
  efficiency = -1;
  nHits = -1;
  if(!CounterParticle.empty()){
    std::vector<std::pair<LHCb::MCParticle*,int> > vpar(std::begin(CounterParticle),std::end(CounterParticle));
    std::sort(std::begin(vpar),std::end(vpar),[](const std::pair<LHCb::MCParticle*,int>& a, const std::pair<LHCb::MCParticle*,int>& b){return a.second<b.second;});
    std::pair<LHCb::MCParticle*,int> parmax = vpar.back();
    particle = parmax.first;
    efficiency = (double)parmax.second/(double)track.hits().size();
    nHits = parmax.second;
  }
  Counter.clear();
  for(PrHits::iterator hit = track.hits().begin();hit!=track.hits().end();++hit){
    LHCb::FTChannelID id = (*hit)->id().ftID();
    LinkedTo<LHCb::MCParticle> ftLink1(evtSvc(),msgSvc(),LHCb::FTClusterLocation::Default);
    LHCb::MCParticle* part = ftLink1.first(id);
    while(0!=part){
      int key = part->key();
      Counter[key]++;
      CounterParticle[part]++;
      part=ftLink1.next();
    }
  }
  efficiency = -1;
  nHits = -1;
  if(!CounterParticle.empty()){
    std::vector<std::pair<LHCb::MCParticle*,int> > vpar(std::begin(CounterParticle),std::end(CounterParticle));
    std::sort(std::begin(vpar),std::end(vpar),[](const std::pair<LHCb::MCParticle*,int>& a, const std::pair<LHCb::MCParticle*,int>& b){return a.second<b.second;});
    std::pair<LHCb::MCParticle*,int> parmax = vpar.back();
    particle = parmax.first;
    efficiency = (double)parmax.second/(double)track.hits().size();
    nHits = parmax.second;
  }
  Counter.clear();
  for(PrHits::iterator hit = track.hits().begin();hit!=track.hits().end();++hit){
    LHCb::FTChannelID id = (*hit)->id().ftID();
    LinkedTo<LHCb::MCParticle> ftLink1(evtSvc(),msgSvc(),LHCb::FTClusterLocation::Default);
    LHCb::MCParticle* part = ftLink1.first(id);
    while(0!=part){
      int key = part->key();
      Counter[key]++;
      CounterParticle[part]++;
      part=ftLink1.next();
    }
  }
  efficiency = -1;
  nHits = -1;
  if(!CounterParticle.empty()){
    std::vector<std::pair<LHCb::MCParticle*,int> > vpar(std::begin(CounterParticle),std::end(CounterParticle));
    std::sort(std::begin(vpar),std::end(vpar),[](const std::pair<LHCb::MCParticle*,int>& a, const std::pair<LHCb::MCParticle*,int>& b){return a.second<b.second;});
    std::pair<LHCb::MCParticle*,int> parmax = vpar.back();
    particle = parmax.first;
    efficiency = (double)parmax.second/(double)track.hits().size();
    nHits = parmax.second;
  }
  Counter.clear();
  // Require the 90% of the hits on the track to be truth matched to a MCParticle for the Assoc flag. \
  
  // Anyway you can look directly to this value since the method returns an updated efficiency value \
  
  if( efficiency>0.9 )return true;
  if( efficiency<0.9 )return false;
  return false;
}





#endif
