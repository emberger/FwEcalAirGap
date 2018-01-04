#include "EResShowerProf.hh"




int main(int argc, char* argv[]){

        EResShowerProf E;

        std::string p=argv[2];
        Double_t cut=std::stod(argv[1]);

        E.EResAndProfile(cut, p);

        return 0;
}
