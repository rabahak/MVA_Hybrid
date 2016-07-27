#include "TFile.h"
#include "TCut.h"
#include "TTree.h"
#include <iostream>

// What do we do : this is a little script to reduce the data for MVA input
// and count the luminosity just because we can.

void hybrid_trainandtest_treeGeneration(){cout<<"===Setting trees started==="<<endl;
  cout<<"Methods to call is: "<<endl;
 
  cout<<"SB_Reduction()"<<endl;
}

void SB_Reduction(){


  //Here you read the root file                                                                                                                                                   

  // TFile * f = new TFile ("/Users/rabah/Desktop/Git_python/lambdab/mc-15154001.root");                                                                                          
  TFile * f = new TFile ("/Users/rabah/Desktop/Hybrid_MVA/inFlagHits/Allev_CasesSeparation_recoveryOFF/case0/Allev_recoveryOFF0.root");
  assert(f);
  // Here you load the tree.                                                                                                                                                      

  //  TTree * tree = (TTree*)f->Get("Tuple_Bu2LLK_eeLine2/DecayTree;1");                                                                                                          
  TTree * tree = (TTree*)f->Get("PrHybridSeeding/XYtracks");
  assert(tree);
  // Here you ask how many evenys you have in your tree                                                                                                                          
  // cout << "Total number of entries : "  << tree->GetEntries () << endl;                                                                                                        



  // Signal/Background categories                                                                                                                                                 
  TCut b1="assoc==0";
  TCut b2="Case==0";
  TCut bcut = b1 && b2; //background

  TCut s1="assoc==1";
  TCut s2="Case==0";
  TCut scut = s1  && s2 ; //signal                                                                                                                               

  //  cout <<"For now we use the following cuts : " << ApplyThis << endl;                                                                                                         
  //TFile * signal_file = new TFile("S_Lambdab_muons.root", "recreate");                                                                                                          
  TFile* signal_file = TFile::Open("SignalXY.root","recreate");
  TTree * signal_copy = tree->CopyTree( scut);
  signal_copy->Write(0,TObject::kOverwrite);
  signal_file->Write(0,TObject::kOverwrite);
  signal_file->Close();
  /// Same thing for the background                                                                                                                                               
  //TFile * backg_file = new TFile("B_Lambdab_muons.root", "recreate");                                                                                                           
  TFile* backg_file = TFile::Open("BackgroundXY.root","recreate");
  TTree * backg_copy = tree->CopyTree( bcut);
  backg_copy->Write(0,TObject::kOverwrite);
  backg_file->Write(0,TObject::kOverwrite);
  backg_file->Close();

  cout << "========We are done with the first reduction.======== " << endl;
  return;
}
/*
//NOT FINISHED
void var_Reduction(){
  ///Users/rabah/Desktop/Hybrid_MVA/test.root
  TFile * signal_file = new TFile ("/Users/rabah/Desktop/Hybrid_MVA/SignalXY.root");
  TFile * backg_file = new TFile ("/Users/rabah/Desktop/Hybrid_MVA/BackgroundXY.root");
  assert(signal_file);
  assert(backg_file);
  TTree * signal = (TTree*)signal_file->Get("XYtracks");
  assert(signal);
  TTree * background = (TTree*)backg_file->Get("XYtracks");
  
  assert(background);

  Bool_t assoc_b;
  Double_t Case_i,nbSingleX_i;
  Double_t dRatio_f,Chi2PerDoF_f,ax_f,bx_f,cx_f,ay_f,by_f,X0Back_f;

 
  
  signal->SetBranchAddress("Case_i", &Case_i);
  signal->SetBranchAddress("nbSingleX_i", &nbSingleX_i);
 
  signal->SetBranchAddress("dRatio_f", &dRatio_f);
  signal->SetBranchAddress("Chi2PerDoF_f", &Chi2PerDoF_f);
  signal->SetBranchAddress("ax_f", &ax_f);
  signal->SetBranchAddress("bx_f", &bx_f);
  signal->SetBranchAddress("cx_f", &cx_f);
  signal->SetBranchAddress("ay_f", &ay_f);
  signal->SetBranchAddress("by_f", &by_f);
  signal->SetBranchAddress("X0Back_f", &X0Back_f);

  

 

  background->SetBranchAddress("Case_i", &Case_i);
  background->SetBranchAddress("nbSingleX_i", &nbSingleX_i);

  background->SetBranchAddress("dRatio_f", &dRatio_f);
  background->SetBranchAddress("Chi2PerDoF_f", &Chi2PerDoF_f);
  background->SetBranchAddress("ax_f", &ax_f);
  background->SetBranchAddress("bx_f", &bx_f);
  background->SetBranchAddress("cx_f", &cx_f);
  background->SetBranchAddress("ay_f", &ay_f);
  background->SetBranchAddress("by_f", &by_f);
  background->SetBranchAddress("X0Back_f", &X0Back_f);
 
  Int_t strain_events;
  Int_t stest_events;
  Int_t btrain_events;
  Int_t btest_events;
  Bool_t testtrain_choice;

 

  cout<<"There's "<<signal->GetEntries()<<" signal events and "<<background->GetEntries()<<" background events."<<endl;
  cout<<"if you want to split up in half press 1, if no press 0: ";cin>>testtrain_choice;
  if(testtrain_choice)
    {
      if(signal->GetEntries()%2==0)
	{
	  strain_events=signal->GetEntries()/2;
	  stest_events=signal->GetEntries()/2;
	}
      else
	{
	  strain_events=(signal->GetEntries()-1)/2;
          stest_events=(signal->GetEntries()-1)/2;
	}

      if(background->GetEntries()%2==0)
        {
          btrain_events=background->GetEntries()/2;
          btest_events=background->GetEntries()/2;
        }
      else
	{
          btrain_events=(background->GetEntries()-1)/2;
          btest_events=(background->GetEntries()-1)/2;
        }

    }
  else
    {
      cout<<"number of Signal events requested for train: ";cin>>strain_events;
      cout<<"number of Signal events requested for test: ";cin>>stest_events;
      cout<<"number of Backg  events requested for train: ";cin>>btrain_events;
      cout<<"number of Backg  events requested for test: ";cin>>btest_events;
    }
  if(strain_events+stest_events>signal->GetEntries()){cout<<"Impossible signal request !"<<endl;return;}
  if(btrain_events+btest_events>background->GetEntries()){cout<<"Impossible backg request !"<<endl;return;}


  TTree * new_signaltrain=new TTree("signaltrain","signal train tree");
  TTree * new_signaltest=new TTree("signaltest","signal test tree");
  TTree * new_backgroundtrain=new TTree("backgroundtrain","background train tree");
  TTree * new_backgroundtest=new TTree("backgroundtest","background test tree");
  new_signaltrain->SetDirectory(0);
  new_signaltest->SetDirectory(0);
  new_backgroundtrain->SetDirectory(0);  
  new_backgroundtest->SetDirectory(0);
  //==SIGNAL
  new_signaltrain->SetBranchAddress("Case_i", &Case_i);
  new_signaltrain->SetBranchAddress("nbSingleX_i", &nbSingleX_i);

  new_signaltrain->SetBranchAddress("dRatio_f", &dRatio_f);
  new_signaltrain->SetBranchAddress("Chi2PerDoF_f", &Chi2PerDoF_f);
  new_signaltrain->SetBranchAddress("ax_f", &ax_f);
  new_signaltrain->SetBranchAddress("bx_f", &bx_f);
  new_signaltrain->SetBranchAddress("cx_f", &cx_f);
  new_signaltrain->SetBranchAddress("ay_f", &ay_f);
  new_signaltrain->SetBranchAddress("by_f", &by_f);
  new_signaltrain->SetBranchAddress("X0Back_f", &X0Back_f);

  new_signaltest->SetBranchAddress("Case_i", &Case_i);
  new_signaltest->SetBranchAddress("nbSingleX_i", &nbSingleX_i);

  new_signaltest->SetBranchAddress("dRatio_f", &dRatio_f);
  new_signaltest->SetBranchAddress("Chi2PerDoF_f", &Chi2PerDoF_f);
  new_signaltest->SetBranchAddress("ax_f", &ax_f);
  new_signaltest->SetBranchAddress("bx_f", &bx_f);
  new_signaltest->SetBranchAddress("cx_f", &cx_f);
  new_signaltest->SetBranchAddress("ay_f", &ay_f);
  new_signaltest->SetBranchAddress("by_f", &by_f);
  new_signaltest->SetBranchAddress("X0Back_f", &X0Back_f);




  new_backgroundtrain->SetBranchAddress("Case_i", &Case_i);
  new_backgroundtrain->SetBranchAddress("nbSingleX_i", &nbSingleX_i);

  new_backgroundtrain->SetBranchAddress("dRatio_f", &dRatio_f);
  new_backgroundtrain->SetBranchAddress("Chi2PerDoF_f", &Chi2PerDoF_f);
  new_backgroundtrain->SetBranchAddress("ax_f", &ax_f);
  new_backgroundtrain->SetBranchAddress("bx_f", &bx_f);
  new_backgroundtrain->SetBranchAddress("cx_f", &cx_f);
  new_backgroundtrain->SetBranchAddress("ay_f", &ay_f);
  new_backgroundtrain->SetBranchAddress("by_f", &by_f);
  new_backgroundtrain->SetBranchAddress("X0Back_f", &X0Back_f);
  
  new_backgroundtest->SetBranchAddress("Case_i", &Case_i);
  new_backgroundtest->SetBranchAddress("nbSingleX_i", &nbSingleX_i);

  new_backgroundtest->SetBranchAddress("dRatio_f", &dRatio_f);
  new_backgroundtest->SetBranchAddress("Chi2PerDoF_f", &Chi2PerDoF_f);
  new_backgroundtest->SetBranchAddress("ax_f", &ax_f);
  new_backgroundtest->SetBranchAddress("bx_f", &bx_f);
  new_backgroundtest->SetBranchAddress("cx_f", &cx_f);
  new_backgroundtest->SetBranchAddress("ay_f", &ay_f);
  new_backgroundtest->SetBranchAddress("by_f", &by_f);
  new_backgroundtest->SetBranchAddress("X0Back_f", &X0Back_f);

//==BACKGROUND
 

  
  for(Int_t i = 0; i<signal->GetEntries(); i++)   
    {signal->GetEntry(i);
      if(i<strain_events)
	new_signaltrain->Fill();
      else
	if(i<strain_events+stest_events)
	  new_signaltest->Fill();
    }
  //  for(Int_t i = 0; i<s_events; i++)  new_signal->Fill();
  for(Int_t i = 0; i<background->GetEntries(); i++)  
    {background->GetEntry(i);
      if(i<btrain_events)
	new_backgroundtrain->Fill();
      else
	if(i<btrain_events+btest_events)
	  new_backgroundtest->Fill();
    }
  //for(Int_t i = 0; i<b_events; i++)  new_background->Fill();
  TFile * strain=TFile::Open("signaltrain.root","recreate");
  new_signaltrain->Write(0,TObject::kOverwrite);
  strain->Close();

  TFile * stest=TFile::Open("signaltest.root","recreate");
  new_signaltest->Write(0,TObject::kOverwrite);
  stest->Close();
 
  TFile * btrain=TFile::Open("backgroundtrain.root","recreate");
  new_backgroundtrain->Write(0,TObject::kOverwrite);
  btrain->Close();
  
  TFile * btest=TFile::Open("backgroundtest.root","recreate");
  new_backgroundtest->Write(0,TObject::kOverwrite);
  btest->Close();

  return;
}
*/
