#include <KM3RunManager.hh>
#include <KM3RunMessenger.hh>

#include <iostream>

using namespace std;

namespace km3net
{

  KM3RunManager::KM3RunManager ()
  {
    fMessenger = new KM3RunMessenger (this);
  }

  KM3RunManager::~KM3RunManager ()
  {
    delete (fMessenger);
  }


}
