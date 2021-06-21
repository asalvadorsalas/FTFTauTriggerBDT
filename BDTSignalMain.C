#include "TCanvas.h"
#include "TROOT.h"
#include "TGraphErrors.h"
#include "TF1.h"
#include "TLegend.h"
#include "TArrow.h"
#include "TLatex.h"
#include <vector>
#include "BDTSignalClass.cxx"

#ifndef __CINT__       
int main(int argc, char* argv[]){
  const char* inputfile = (argc > 1)? argv[1] : "test.root";

  std::cout << "Inputfile    : "<< inputfile << std::endl;
  BDTSignalClass Ana(inputfile,"");
  Ana.Loop();

  return 0;     
}
#endif