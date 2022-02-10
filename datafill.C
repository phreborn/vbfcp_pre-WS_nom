#include "utils.h"

using namespace TMVA::Experimental;

void datafill()
{

  char *cf_cats = (char*)"cats.cfg";
  //char *cf_bins = (char*)"binnings.cfg";
  map<TString, string> catCuts;
  //map<TString, vector<float>> catBins;
  getCatCuts(cf_cats, catCuts); for(auto c : catCuts) cout<<c.first<<c.second<<endl;
  //getCatBinning(cf_bins, catBins); for(auto c : catBins) { cout<<"dafad"<<c.first<<endl; for(auto b : c.second)  cout<<b<<", "; cout<<endl<<endl; }

  string config = "config";
  string baseCuts = "";
  string blindCut = "";
  readConfigFile(config.data(), "cuts", baseCuts);
  readConfigFile(config.data(), "blindSel", blindCut);

  baseCuts = TString(baseCuts.data()).ReplaceAll("??", "").Data();
  blindCut = TString(blindCut.data()).ReplaceAll("??", "").Data();

  bool doBlind = true;
  string blind_tmp = "";
  readConfigFile(config.data(), "doBlind", blind_tmp);
  if(blind_tmp.find("TRUE") != std::string::npos) doBlind = true;
  else doBlind = false;

  map<TString, pair<float, float>> bins;
  bins["b1"] = make_pair(-999999999, -2);
  bins["b2"] = make_pair(-2, -1);
  bins["b3"] = make_pair(-1, 0);
  bins["b4"] = make_pair(0, 1);
  bins["b5"] = make_pair(1, 2);
  bins["b6"] = make_pair(2, 99999999);

  std::map<TString, std::vector<float>> cats;
  cats["TT"] = {0.14, 1., 0.23, 1.};
  cats["TL"] = {0.14, 1., -1., 0.23};
  cats["LT"] = {-1, 0.14, 0.05, 1.};
  cats["LL"] = {-1, 0.14, -1., 0.05};

  TString dirpath = "/scratchfs/atlas/huirun/atlaswork/VBF_CP/h026/";
  std::string path_str = dirpath.Data();
  std::vector<std::string> sub_dirs = getDirBinsSortedPath(path_str);

  std::vector<std::string> files(0);

  for(auto d : sub_dirs){
    if(d==".") continue;
    if(d.find("data") == std::string::npos) continue;
    if(d.find(".root") != std::string::npos) continue;
    cout<<"d: "<<path_str+d<<endl;
    std::vector<std::string> fs = getDirBinsSortedPath(path_str+d+"/");
    for(auto f : fs){
      if(f==".") continue;
      if(f.find("slim") != std::string::npos) continue;
      cout<<"f: "<<path_str+"/"+d+"/"+f<<endl;
      files.push_back(path_str+d+"/"+f);
      
    }
  }

  cout<<"\n"<<endl;

//  TString sample = "dataSB";
  TString sample = "fulldata";

  cout<<"processing "<<sample<<".."<<endl;

  Int_t N_j,N_photon,cutflow,Category;
  Float_t m_yy,pT_y1,pT_y2,m_jj_30,DeltaEta_jj,Zepp,oo1,oo2,WeightDtilde1,WeightDtilde2,weight,xsec_kF_eff,total_weights;
  Bool_t isDalitz,isPassed,isPassedIsolation,isPassedPID,isPassedTriggerMatch;

  for(auto cat : catCuts){
    TString catName = cat.first; cout<<"=== "<<catName<<" ==="<<endl;
    string catCut = cat.second;

    catCut = TString(catCut.data()).ReplaceAll("??", "").Data();

    TH1F *h_myy = new TH1F("m_yy_"+sample, "", 110, 105, 160);
    TH1F *h_oo = new TH1F("oo_"+sample, "", 40, -20, 20);

    TChain ch("output", "output");
  
    for(auto f : files){
      TString filepath = f.data();
  
      cout<<filepath<<endl;

      ch.Add(filepath);
 

    }

    ROOT::RDataFrame df(ch, {"m_yy"});

    string allCuts = "";
    string allCuts_noBlind = Form("%s && %s", baseCuts.data(), catCut.data()); cout<<"all cuts: "<<allCuts<<endl;
    if(doBlind) allCuts = Form("%s && %s && %s", baseCuts.data(), catCut.data(), blindCut.data());
    else allCuts = allCuts_noBlind;

    auto df_cut = df.Filter(allCuts);
    df_cut.Foreach([&h_myy, &h_oo](float myy, float oo){ h_myy->Fill(myy/1000); h_oo->Fill(oo); }, {"m_yy", "oo1"}); cout<<df_cut.Sum("weight").GetValue()<<endl;
    df_cut.Snapshot("CollectionTree", Form("tree_%s_OO_%s.root", sample.Data(), catName.Data()), {"weight", "m_yy", "oo1"});
    //df_cut.Snapshot("CollectionTree", Form("hists_%s_%s.root", sample.Data(), catName.Data()), {"weight", "m_yy", "oo1", "BDTout_ggH", "BDTout_yy"});
  
    TFile *f_out = new TFile(Form("hists_%s_%s.root", sample.Data(), catName.Data()),"update");

    f_out->cd();
    h_myy->Write();
    h_oo->Write();
    //f_out->Close();
  
    delete h_myy;
    delete h_oo;

  }

}
