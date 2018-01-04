#include "CalcAngularRes.hh"


int main(int argc, char * argv[]) {

  AngularResolution Res;

  Res.CalcAngularResolution(argv[1], argv[2]);

  return 0;
}
