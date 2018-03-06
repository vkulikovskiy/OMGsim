#include <TFormula.h>

#include <KM3ParametrizedVolume.hh>
#include <sys/types.h>
#include <dirent.h>
#include <iostream>
#include <libconfig.h++>
#include <regex>

#include <G4VisAttributes.hh>


using namespace std;
using namespace libconfig;


namespace km3net
{

  std::map <std::string, long>         KM3ParametrizedVolume::fIntegerParameters;
  std::map <std::string, double>       KM3ParametrizedVolume::fFloatingParameters;
  std::map <std::string, std::string>  KM3ParametrizedVolume::fStringParameters;
  std::map <std::string, std::string>  KM3ParametrizedVolume::fFormulaParameters;


  G4VisAttributes* KM3ParametrizedVolume::fgold        = new G4VisAttributes(G4Colour (255./255,185./255,15./255,1));
  G4VisAttributes* KM3ParametrizedVolume::fsilver      = new G4VisAttributes(G4Colour (192./255,192./255,192./255,1));
  G4VisAttributes* KM3ParametrizedVolume::fbrown      = new G4VisAttributes(G4Colour (102./255,51./255,0./255,1));
  G4VisAttributes* KM3ParametrizedVolume::forange      = new G4VisAttributes(G4Colour (255./255,100./255,0./255,1));
  G4VisAttributes* KM3ParametrizedVolume::fcyan        = new G4VisAttributes(G4Colour(0.,1.,1.));
  G4VisAttributes* KM3ParametrizedVolume::ftransparent = new G4VisAttributes(G4Colour(1.,1.,1.,1.));

  //string KM3ParametrizedVolume::ConfigFilePrefix = "          ";

  KM3ParametrizedVolume::KM3ParametrizedVolume ()
  {
    fFormula = new TFormula;
  }

  KM3ParametrizedVolume::~KM3ParametrizedVolume ()
  {
    delete fFormula;
  }


  void KM3ParametrizedVolume::SetParametersListFromMac (string SParList) {
    size_t doubledotPos=SParList.find (';');
    string subSParList=SParList;
    while (doubledotPos != string::npos) {
      subSParList=SParList.substr (0,doubledotPos);
      cout << "match to search " << subSParList << endl;
      ReadAStringParameterFromMac (subSParList);
      SParList.erase (0,doubledotPos+1);
      doubledotPos=SParList.find(';');
    }
    if (SParList.size ()) {
      ReadAStringParameterFromMac (SParList);
      subSParList=SParList.substr (0,doubledotPos);
    }
  }

  void KM3ParametrizedVolume::ReadAStringParameterFromMac (std::string SPar) {
    std::cmatch m;
    std::regex_match ( SPar.data (), m, std::regex(" *([a-zA-Z_0-9]+) *= *([0-9]+) *") );
    string Stmp;
    for (unsigned i=1; i<m.size(); ++i)
      {
        if (i%2 == 1){
          std::cout << "match: adding to parameter list: " << m[i];
          Stmp = m[i];
        }else {
          cout << " with value " << stof(m.str(i)) << endl;
          fFloatingParameters[Stmp] = stof(m.str(i));
        }
      }
  }

  void KM3ParametrizedVolume::SetDataFilesSources (string Path) {
    size_t doubledotPos=Path.find (':');
    string subPath=Path;
    while (doubledotPos != string::npos) {
      subPath=Path.substr (0,doubledotPos);
      LoadDirectory (subPath);
      Path.erase (0,doubledotPos+1);
      doubledotPos=Path.find(':');
    }
    if (Path.size ())
      LoadDirectory (subPath);
  }

  void KM3ParametrizedVolume::LoadDirectory (string DIRNAME)
  {
    struct dirent *de=NULL;
    DIR *d=NULL;

    d=opendir(DIRNAME.data ());
    if(d == NULL)
      {
        cerr << "Couldn't open directory" << DIRNAME << endl;
        return;
      }

    // Loop while not NULL
    while( (de = readdir(d)) )
      {
        string filename=de->d_name;
        if (de->d_type == 8 && filename.substr (0,ConfigFilePrefix.size ()) == ConfigFilePrefix &&
            filename.rfind (".dat") == filename.size () - 4)
          {
            clog << "KM3ParametrizedVolume: Loading file " << DIRNAME  << '/' << de->d_name<< endl;
            ReadParameterList (DIRNAME+"/"+filename);
            ReadConfiguration (DIRNAME+"/"+filename);
          }
      }
    closedir(d);
  }

  bool KM3ParametrizedVolume::ReadBaseParameters (Config* cfg, VolumeParameters* thePTP)
  {
    return ReadBaseParameters (&cfg->getRoot (), thePTP);
  }

  bool KM3ParametrizedVolume::ReadBaseParameters (Setting* stg, VolumeParameters* thePTP)
  {
    stg->lookupValue("CheckOverlaps"      , thePTP->CheckOverlaps);
    return
      stg->lookupValue("Prefix"       , thePTP->Prefix);
  }

  vector <double> KM3ParametrizedVolume::ReadThisArray (const libconfig::Setting& STG) const
  {
    vector <double> tmpVect;
    if  (STG.getType () != Setting::TypeArray)
      {
        cerr << "ERROR: KM3ParametrizedVolume::ReadThisArray: I'm asked to read a setting which is not an array, the SettingType is " << STG.getType () << ", I'll do nothing."<< endl;
        return tmpVect;
      }
    for (int index=0; index<STG.getLength (); index++)
      tmpVect.push_back (STG[index]);

    return tmpVect;
  }

  void KM3ParametrizedVolume::ReadParameterList (std::string FILENAME)
  {
    string ParSuffix="__";

    Config cfg;
    cfg.setAutoConvert (true);
    cfg.readFile (FILENAME.data ());

    Setting& stgRoot= cfg.getRoot ();

    int SettingSize=stgRoot.getLength ();

    for (int it=0; it<SettingSize; it++)//Setting::iterator itr=stgRoot.begin(); itr!=stgRoot.end(); itr++)
        {
          Setting& CurrentSetting=stgRoot[it];
          string SettingName=CurrentSetting.getName ();
          int strSize=SettingName.size ();
          if (strSize > 2 && SettingName.substr(strSize-2,strSize-1) == ParSuffix)
            {
              if (CurrentSetting.getType () == Setting::TypeInt ||
                  CurrentSetting.getType () == Setting::TypeInt64)
                fIntegerParameters[SettingName] = long(CurrentSetting);
              else if (CurrentSetting.getType () == Setting::TypeFloat)
                fFloatingParameters[SettingName] = double(CurrentSetting);
              else if (CurrentSetting.getType () == Setting::TypeString)
                {
                  string SubFormula = isAFormula (CurrentSetting);
                  if (SubFormula.size ())
                    fFormulaParameters[SettingName]=SubFormula;//ReadThisFormula (CurrentSetting);
                  else
                    fStringParameters[SettingName] = (const char*)CurrentSetting;
                }
              else if (CurrentSetting.getType () == Setting::TypeArray)
                fScalarArrayParameters[SettingName]=ReadThisArray (CurrentSetting);
              //FormulaLooping ();
            }
        }
  }

  void KM3ParametrizedVolume::FormulaLooping ()
  {
    for (auto& aStrPair: fStringParameters)
      {
        string SubFormula = isAFormula (aStrPair.second);
        if (SubFormula.size ())
          {
            fFormulaParameters[aStrPair.first]=SubFormula;//ReadThisFormula (CurrentSetting);
            fStringParameters.erase (aStrPair.first);
          }
      }
  }

  string KM3ParametrizedVolume::isAFormula (const Setting& STG)
  {
    string aStr = (const char*)STG;
    return isAFormula (aStr);
  }

  string KM3ParametrizedVolume::isAFormula (string& STR)
  {
    string sformula = ReadThisFormula (STR);
    int isOk = fFormula->Compile (sformula.data ());
    if (isOk == 0)
      return sformula;
    return "";

  }
  string KM3ParametrizedVolume::ReadThisFormula (const Setting& STG)
  {
    string strToTest = (const char*)STG;
    return ReadThisFormula (strToTest);
  }

  string KM3ParametrizedVolume::ReadThisFormula (string& strToTest)
  {
  beginWhileTest:
    if (strToTest.find ("__") != string::npos)
      {
        for (const auto& parPair: fFloatingParameters)
          {
            size_t posToReplace=0;
            posToReplace = strToTest.find (parPair.first);
            if (posToReplace == string::npos)
              continue;
            string substitute = to_string (parPair.second);
            cout << strToTest << " " << substitute << " ";
            strToTest.replace (posToReplace, parPair.first.size (), substitute);
            goto beginWhileTest;
          }
        for (const auto& parPair: fIntegerParameters)
          {
            size_t posToReplace=0;
            posToReplace = strToTest.find (parPair.first);
            if (posToReplace == string::npos)
              continue;
            string substitute = to_string (parPair.second);
            strToTest.replace (posToReplace, parPair.first.size (), substitute);
            goto beginWhileTest;
          }
        for (const auto& parPair: fFormulaParameters)
          {
            size_t posToReplace=0;
            posToReplace = strToTest.find (parPair.first);
            if (posToReplace == string::npos)
              continue;
            string substitute = parPair.second;
            strToTest.replace (posToReplace, parPair.first.size (), substitute);
            goto beginWhileTest;
          }
      }
    return strToTest;
  }

  vector<double> KM3ParametrizedVolume::GetScalarArrayParameter (const Setting& STG) const
  {
    if (STG.getType () == Setting::TypeArray)
      return ReadThisArray(STG);

    if (STG.getType () != Setting::TypeString)
      {
        //cerr << "WARNING:" << STG.getName () << " parameter is not and int neither a string, return empty vector."<<endl;
        return vector<double>();
      }

    string Name = (const char*) STG;

    if (fScalarArrayParameters.count(Name))
      {
        return fScalarArrayParameters.at(Name);
      }

    cerr << "WARNING:" << STG.getName () << " does not exist in table, might not be inizializated, or the parameter name does not end with \"__\"" << endl;
    return vector<double>();
  }

  long KM3ParametrizedVolume::GetIntegerParameter (const Setting& STG, float x, float y, float z)
  {
    if (STG.getType () == Setting::TypeInt ||
        STG.getType () == Setting::TypeInt64)
      return STG;
    if (STG.getType () != Setting::TypeString)
      {
        cerr << "WARNING:" << STG.getName () << " parameter is not and int neither a string, return max of double."<<endl;
        return 0xffffffffffffffff;
      }

    string Name = (const char*) STG;

    if (fIntegerParameters.count(Name))
      {
        return fIntegerParameters[Name];
      }

    if (fFormulaParameters.count(Name))
      {
        fFormulaParameters[STG.getName ()]=ReadThisFormula (STG);
        fFormula->Compile (fFormulaParameters[STG.getName ()].data ());
        return fFormula->Eval (x,y,z);
      }

    cerr << "WARNING:" << STG.getName () << " does not exist in table, might not be inizializated, or the parameter name does not end with \"__\"" << endl;
    return 0xffffffffffffffff;
  }

  double KM3ParametrizedVolume::GetFloatingParameter (const Setting& STG, float x, float y, float z)
  {
    if (STG.getType () == Setting::TypeFloat ||
        STG.getType () == Setting::TypeInt ||
        STG.getType () == Setting::TypeInt64)
      return STG;
    if (STG.getType () != Setting::TypeString)
      {
        cerr << "WARNING:" << STG.getName () << " parameter is not and double/float neither a string, return max of double."<<endl;
        return 0xffffffffffffffff;
      }

    string Name = (const char*) STG;

    if (fFloatingParameters.count(Name))
      {
        return fFloatingParameters[Name];
      }

    if (fIntegerParameters.count(Name))
      {
        return fIntegerParameters[Name];
      }

    if (fFormulaParameters.count(Name))
      {
        fFormula->Compile (fFormulaParameters[Name].data ());
        return fFormula->Eval (x,y,z);
      }

    if (STG.getType () == Setting::TypeString)
      {
        string SubFormula = isAFormula (STG);
        if (SubFormula.size ())
          {
            fFormulaParameters[STG.getName ()]=ReadThisFormula (STG);
            fFormula->Compile (fFormulaParameters[STG.getName ()].data ());
            return fFormula->Eval (x,y,z);
          }
      }


    cerr << "WARNING:" << Name << " does not exist in table, might not be inizializated, or the parameter name does not end with \"__\"" << endl;
    return 0xffffffffffffffff;
  }

  string KM3ParametrizedVolume::GetStringParameter (const Setting& STG)
  {
    if (STG.getType () != Setting::TypeString)
      {
        cerr << "WARNING:" << STG.getName () << " parameter is not a string."<<endl;
        return "";
      }

    string SValue = (const char*)STG;

    if (!fStringParameters.count(SValue))
      {
        return SValue;
      }
    return fStringParameters[SValue];
  }

  long KM3ParametrizedVolume::GetIntegerParameter (string NAME)
  {
    if (fIntegerParameters.count(NAME))
      return fIntegerParameters[NAME];
    return 0xffffffffffffffff;
  }

  double KM3ParametrizedVolume::GetFloatingParameter (string NAME)
  {
    if (fFloatingParameters.count(NAME))
      return fFloatingParameters[NAME];
    return 0xffffffffffffffff;
  }

  string KM3ParametrizedVolume::GetFormulaParameter (string NAME)
  {
    if (fFormulaParameters.count(NAME))
      return fFormulaParameters[NAME];
    return "";
  }

  string KM3ParametrizedVolume::GetStringParameter (string NAME)
  {
    if (fStringParameters.count(NAME))
      return fStringParameters[NAME];
    return "";
  }

  const KM3ParametrizedVolume::KM3Volume* KM3ParametrizedVolume::GetVolume (string NAME)
  {
    if (!fParametersDict.count (NAME))
      {
        cerr << "KM3ParametrizedVolume::KM3Volume warning, cannot create " << NAME <<", it does not exist in the dictionary. return NULL."<< endl;
        return NULL;
      }
    if (!fVolumeDict.count (NAME))
      {
        fVolumeDict[NAME] = ConstructTheVolume (fParametersDict[NAME]);
      }
    return fVolumeDict[NAME];
  }

  G4LogicalVolume* KM3ParametrizedVolume::GetLogicalVolume (string NAME)
  {
    if (!fVolumeDict.count(NAME))
      {
        if (this->GetVolume (NAME) == NULL)
          return NULL;
      }
    G4LogicalVolume* toreturn = fVolumeDict[NAME]->Construct ();
    if (!fLogicalVolDict.count(toreturn))
      fLogicalVolDict[toreturn]=NAME;
    return toreturn;
  }

  KM3ParametrizedVolume::KM3Volume::~KM3Volume ()
  {
  }

}
