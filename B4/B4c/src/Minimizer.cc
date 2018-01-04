#include "Minimizer.hh"


double Fcn::operator()(const std::vector<double> &par) const
{

        double chisq = 0.;
        //std::cout<<COGs[currentEvent].size()<<std::endl;
        if(mode=="photon") {
                if(seg==1 || seg==2) {
                        for(int i = 0; i<COGs.size(); i++) {

                                double termx = 0;
                                double termy = 0;

                                termx += ((std::get<0>(COGs[i]) - (par[0] * std::get<2>(COGs[i])+par[1])) / std::get<3>(COGs[i]))
                                         *std::get<5>(COGs[i]);

                                termy += ((std::get<1>(COGs[i]) - (par[2] * std::get<2>(COGs[i])+par[3])) / std::get<4>(COGs[i]))
                                         *std::get<5>(COGs[i]);
                                chisq +=termx*termx+termy*termy;
                        }
                }
                if(seg==3||seg==4) {
                        for(int i = 0; i<COGs.size(); i++) {
                                double termx = 0;
                                double termy = 0;

                                termx += ((std::get<0>(COGs[i]) - (par[0] * std::get<1>(COGs[i])+par[1])) / std::get<4>(COGs[i]))
                                         *std::get<5>(COGs[i]);

                                termy += ((std::get<2>(COGs[i]) - (par[2] * std::get<1>(COGs[i])+par[3])) / std::get<3>(COGs[i]))
                                         *std::get<5>(COGs[i]);
                                chisq +=termx*termx+termy*termy;
                        }

                }
                if(seg==5||seg==6) {
                        for(int i = 0; i<COGs.size(); i++) {
                                double termx = 0;
                                double termy = 0;

                                termx += ((std::get<1>(COGs[i]) - (par[0] * std::get<0>(COGs[i])+par[1])) / std::get<3>(COGs[i]))
                                         *std::get<5>(COGs[i]);

                                termy += ((std::get<2>(COGs[i]) - (par[2] * std::get<0>(COGs[i])+par[3])) / std::get<4>(COGs[i]))
                                         *std::get<5>(COGs[i]);
                                chisq +=termx*termx+termy*termy;
                        }
                }
        }

        if(mode=="photon3D") {

                for(int i = 0; i<COGs.size(); i++) {
                        Double_t dist=0;
                        //  Double_t k= ((std::get<0>(COGs[i])-par[0])/(par[3]))+((std::get<1>(COGs[i])-par[1])/(par[4]))+((std::get<2>(COGs[i])-par[2])/(par[5]));
                        Double_t k=(-1*par[0]*par[3] + std::get<0>(COGs[i])*par[3]
                                    -1*par[1]*par[4] + std::get<1>(COGs[i])*par[4]
                                    -1*par[2]*par[5] + std::get<2>(COGs[i])*par[5])/(par[3]*par[3]+par[4]*par[4]+par[5]*par[5]);
                        TVector3 A(par[0]+k*par[3], par[1]+k*par[4],par[2]+k*par[5]);   //calculate point on line, corresponding to perpendicular distance

                        TVector3 X(std::get<0>(COGs[i]),std::get<1>(COGs[i]),std::get<2>(COGs[i]));
                        TVector3 d=A-X;           //calculate distance vector
                        dist=d.Mag();
                        chisq+=dist*std::get<5>(COGs[i]);       //add distance to chi2
                }
                return chisq;
        }

        else if(mode=="pion") {

                // Double_t X_d=((std::get<0>(ParamsPions[0])*par[0] + std::get<1>(ParamsPions[0]))-
                //               (std::get<0>(ParamsPions[1])*par[0] + std::get<1>(ParamsPions[1])));
                //
                // Double_t Y_d=((std::get<2>(ParamsPions[0])*par[0] + std::get<3>(ParamsPions[0]))-
                //               (std::get<2>(ParamsPions[1])*par[0] + std::get<3>(ParamsPions[1])));
                //
                // chisq=TMath::Sqrt( X_d*X_d + Y_d*Y_d);

                Double_t dx=trackph1.first.X()-trackph2.first.X()+par[0]*trackph1.second.X()-par[1]*trackph2.second.X();
                Double_t dy=trackph1.first.Y()-trackph2.first.Y()+par[0]*trackph1.second.Y()-par[1]*trackph2.second.Y();

                Double_t dz=trackph1.first.Z()-trackph2.first.Z()+par[0]*trackph1.second.Z()-par[1]*trackph2.second.Z();

                chisq=dx*dx + dy*dy + dz*dz;

        }

        else if(mode=="approach") {
                Double_t sigPh1=(1/TMath::Sqrt(l_ph1[3]))*delta1;
                Double_t sigPh2=(1/TMath::Sqrt(l_ph2[3]))*delta2;

                //std::cout<<"sig1: "<<sigPh1<<" sig2: "<<sigPh2<<std::endl;


                Double_t delR1=TMath::Sqrt(((l_ph1[0]-par[0])*(l_ph1[0]-par[0])) + ((l_ph1[1]-par[1])*(l_ph1[1]-par[1])) + ((l_ph1[2]-par[2])*(l_ph1[2]-par[2])));
                Double_t delR2=TMath::Sqrt(((l_ph2[0]-par[0])*(l_ph2[0]-par[0])) + ((l_ph2[1]-par[1])*(l_ph2[1]-par[1])) + ((l_ph2[2]-par[2])*(l_ph2[2]-par[2])));

                chisq=(delR1/sigPh1)*(delR1/sigPh1)+(delR2/sigPh2)*(delR2/sigPh2);

        }

        else if(mode=="InvMass") {
                Double_t sigPh1=(1/TMath::Sqrt(l_ph1[3]))*delta1;
                Double_t sigPh2=(1/TMath::Sqrt(l_ph2[3]))*delta2;

                //std::cout<<"sig1: "<<sigPh1<<" sig2: "<<sigPh2<<std::endl;


                Double_t delR1=TMath::Sqrt(((l_ph1[0]-par[0])*(l_ph1[0]-par[0])) + ((l_ph1[1]-par[1])*(l_ph1[1]-par[1])) + ((l_ph1[2]-par[2])*(l_ph1[2]-par[2])));
                Double_t delR2=TMath::Sqrt(((l_ph2[0]-par[0])*(l_ph2[0]-par[0])) + ((l_ph2[1]-par[1])*(l_ph2[1]-par[1])) + ((l_ph2[2]-par[2])*(l_ph2[2]-par[2])));


                TVector3 dir_ph1(std::get<0>(cogph1)-par[0],std::get<1>(cogph1)-par[1],std::get<2>(cogph1)-par[2]);

                dir_ph1=dir_ph1.Unit();

                TLorentzVector lv_ph1(dir_ph1.X()*l_ph1[3],dir_ph1.Y()*l_ph1[3],dir_ph1.Z()*l_ph1[3],l_ph1[3]);



                TVector3 dir_ph2(std::get<0>(cogph2)-par[0],std::get<1>(cogph2)-par[1],std::get<2>(cogph2)-par[2]);

                dir_ph2=dir_ph2.Unit();

                TLorentzVector lv_ph2(dir_ph2.X()*l_ph2[3],dir_ph2.Y()*l_ph2[3],dir_ph2.Z()*l_ph2[3],l_ph2[3]);

                TLorentzVector lv_pi;
                lv_pi=lv_ph1+lv_ph2;

                Double_t curr_invmass=lv_pi.M();

                //std::cout<<par[0]<<":"<<par[1]<<":"<<par[2]<<std::endl;

                chisq=((delR1/sigPh1)*(delR1/sigPh1))+((delR2/sigPh2)*(delR2/sigPh2))+(((curr_invmass-134.9766)/0.001)*((curr_invmass-134.9766)/0.001));

        }

        else if(mode=="InvMass2") {

                TVector3 dir_ph1(par[0]-par[3], par[1]-par[4], par[2]-par[5]);
                dir_ph1=dir_ph1.Unit();

                TLorentzVector lv_ph1(dir_ph1.X()*l_ph1[3],dir_ph1.Y()*l_ph1[3],dir_ph1.Z()*l_ph1[3],l_ph1[3]);

                TVector3 dir_ph2(par[0]-par[6], par[1]-par[7], par[2]-par[8]);
                dir_ph2=dir_ph2.Unit();

                TLorentzVector lv_ph2(dir_ph2.X()*l_ph2[3],dir_ph2.Y()*l_ph2[3],dir_ph2.Z()*l_ph2[3],l_ph2[3]);

                TLorentzVector lv_pi=lv_ph1+lv_ph1;

                Double_t curr_invmass=lv_pi.M();

                //std::cout<<par[0]<<":"<<par[1]<<":"<<par[2]<<std::endl;

                chisq=(((curr_invmass-134.9766)/0.1)*((curr_invmass-134.9766)/0.1));


        }

        return chisq;

}

double Fcn::Up() const
{
        return 1;
}
