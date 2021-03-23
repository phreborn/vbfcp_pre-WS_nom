#include <iterator>
#include <iostream>
#include <vector>
#include <stdio.h>
#include <dirent.h>

#include <sstream>
#include <fstream>

#include <string.h>

using namespace std;

static bool readConfigFile(const char * cfgfilepath, const string & key, string & value)
{
    fstream cfgFile;
    cfgFile.open(cfgfilepath);
    if( ! cfgFile.is_open())
    {
        cout<<"can not open cfg file!"<<endl;
        return false;
    }
    char tmp[1000];
    while(!cfgFile.eof())
    {
        cfgFile.getline(tmp,1000);
        string line(tmp);
        size_t pos = line.find(':');
        if(pos==string::npos) continue;
        string tmpKey = line.substr(0,pos);
        if(key==tmpKey)
        {
            value = line.substr(pos+1);
        }
    }
    return false;
}


int fileNameFilter(const struct dirent *cur) {
    std::string str(cur->d_name);
    if (str.find("..") != std::string::npos) {
        return 0;
    }
    return 1;
}

std::vector<std::string> getDirBinsSortedPath(std::string dirPath) {
    struct dirent **namelist;
    std::vector<std::string> ret;
    int n = scandir(dirPath.c_str(), &namelist, fileNameFilter, alphasort);
    if (n < 0) {
        return ret;
    }
    for (int i = 0; i < n; ++i) {
        std::string filePath(namelist[i]->d_name);
        ret.push_back(filePath);
        free(namelist[i]);
    };
    free(namelist);
    return ret;
}

vector<string> readInLines(const char * cfgfilepath){
    vector<string> readin;

    fstream cfgFile;
    cfgFile.open(cfgfilepath);
    if( ! cfgFile.is_open())
    {
        cout<<"can not open cfg file!"<<endl;
    }
    char tmp[1000];
    while(!cfgFile.eof())
    {
        cfgFile.getline(tmp,1000);
        string line(tmp);
        if(line.find("#") != string::npos) continue;
        if(line=="") continue;
        readin.push_back(line);
    }

    return readin;
}

map<TString, string> sepKeyValue(string cfg){
  map<TString, string> cats;

  vector<string> lines = readInLines(cfg.data());
  for(auto l : lines){
    int pos = l.find(":");
    string sKey = l.substr(0,pos);
    string sValues = l.substr(pos+1);
    TString tsKey = sKey.data();
    tsKey.ReplaceAll(" ", "");

    cats[tsKey] = sValues;
  }

  return cats;
}

void getCatCuts(string cfg, map<TString, string> &catCuts){
  catCuts = sepKeyValue(cfg);
}

void getCatBinning(string cfg, map<TString, vector<float>> &catBins){
  map<TString, string> cats = sepKeyValue(cfg);

  for(auto c : cats){
    string sBinEdges = c.second;
    TString tsCat = c.first;
 
    TString ts_tmp;
    vector<float> f_vec;
    char *p;
    char *buff = (char*) sBinEdges.data();
    char *sep = (char*)",";
    p = strsep(&buff, sep);
    while(p!=NULL){
      ts_tmp = p;
      ts_tmp.ReplaceAll(" ", "");
      if(ts_tmp!="") f_vec.push_back((float)atof(ts_tmp.Data())); //cout<<ts_tmp<<endl;
      p = strsep(&buff, sep);
    }

    catBins[tsCat] = f_vec;
  }
}

TString getMCSampleName(int mcID){
  string name;
  readConfigFile("/scratchfs/atlas/chenhr/atlaswork/VBF_CP/syst/MCSamples.config", Form("SampleName.%d", mcID), name);
  while(name.find(" ")!=std::string::npos) { name.replace(name.find(" "), 1, ""); }
  return name.data();
}

TH1F *getCutFlowHist(int mcID, TFile* file){
  TString suffix = "_noDalitz_weighted";
  TString cutFlowName = Form("CutFlow_%s%s", getMCSampleName(mcID).Data(), suffix.Data());
  TH1F *cutFlow = (TH1F*) file->Get(cutFlowName);
  return cutFlow;
}

double getSumOfWeights(int mcID, TFile* file){
  double NxAOD = getCutFlowHist(mcID, file)->GetBinContent(1);
  double NDxAOD = getCutFlowHist(mcID, file)->GetBinContent(2);
  double WDxAOD = getCutFlowHist(mcID, file)->GetBinContent(3);

  double weightSum = WDxAOD*NxAOD/NDxAOD;
  cout<<"xAOD, DxAOD, allEvt: "<<NxAOD<<", "<<NDxAOD<<", "<<WDxAOD<<endl;
  return weightSum;
}
