// Implementation file for Pi0Event and components
// create dictionary with rootcint RootClasses_dict.cpp -c include/B4ROOTEvent.hh include/LinkDef.h
//

#include "B4ROOTEvent.hh"
//#include "Stiostream.h"




TClonesArray* B4ROOTEvent::aHits = 0;

// TPi0Event class

B4ROOTEvent::B4ROOTEvent() {
								if (!aHits) aHits = new TClonesArray("B4ROOTHit", 1000);
								m_Hits=aHits;

								m_MomentumPh1.SetXYZ(0.,0.,0.);
								m_GunPos.SetXYZ(0.,0.,0.);


								m_EventNo= 0;

								m_GapEnergy= 0;

								m_NHits= 0;

								m_InnerAbsoThickness= 0;
								m_InnerGapThickness= 0;
								m_InnerLayerno= 0;
								m_InnertilesizeX= 0;
								m_InnertilesizeY= 0;
								Innerabsofirst= true;

								m_OuterAbsoThickness= 0;
								m_OuterGapThickness= 0;
								m_OuterLayerno= 0;
								m_OutertilesizeX= 0;
								m_OutertilesizeY= 0;
								Outerabsofirst= true;

								//m_airgapThickness=0;
								m_calsizeXY= 0;
}

B4ROOTEvent::B4ROOTEvent(const B4ROOTEvent& origEvent) : TObject(origEvent) {

								m_InnerGapThickness=origEvent.m_InnerGapThickness;
								m_InnerAbsoThickness=origEvent.m_InnerAbsoThickness;
								m_InnerLayerno=origEvent.m_InnerLayerno;
								m_InnertilesizeX=origEvent.m_InnertilesizeX;
								m_InnertilesizeY=origEvent.m_InnertilesizeY;
								Innerabsofirst=origEvent.Innerabsofirst;

								m_OuterGapThickness=origEvent.m_OuterGapThickness;
								m_OuterAbsoThickness=origEvent.m_OuterAbsoThickness;
								m_OuterLayerno=origEvent.m_OuterLayerno;
								m_OutertilesizeX=origEvent.m_OutertilesizeX;
								m_OutertilesizeY=origEvent.m_OutertilesizeY;
								Outerabsofirst=origEvent.Outerabsofirst;

								m_calsizeXY=origEvent.m_calsizeXY;
								//m_airgapThickness=origEvent.m_airgapThickness;

								m_MomentumPh1=origEvent.m_MomentumPh1;
								m_GunPos=origEvent.m_GunPos;
								m_Hits=origEvent.m_Hits;
								m_EventNo=origEvent.m_EventNo;
								m_GapEnergy=origEvent.m_GapEnergy;
								m_NHits=origEvent.m_NHits;
}

B4ROOTEvent::~B4ROOTEvent() {
								Clear();
}

void B4ROOTEvent::Clear(const Option_t*) {

								m_MomentumPh1.SetXYZ(0.,0.,0.);
								m_GunPos.SetXYZ(0.,0.,0.);
								m_Hits->Clear();

								m_EventNo= 0;

								m_GapEnergy= 0;

								m_NHits= 0;

								m_InnerAbsoThickness= 0;
								m_InnerGapThickness= 0;
								m_InnerLayerno= 0;
								m_InnertilesizeX= 0;
								m_InnertilesizeY= 0;
								Innerabsofirst= true;

								m_OuterAbsoThickness= 0;
								m_OuterGapThickness= 0;
								m_OuterLayerno= 0;
								m_OutertilesizeX= 0;
								m_OutertilesizeY= 0;
								Outerabsofirst= true;
								//m_airgapThickness=0.;
								m_calsizeXY= 0;
}
void B4ROOTEvent::Reset(Option_t * /*option*/)
{
								// Static function to reset all static objects for this event
								delete aHits; aHits = 0;
}




B4ROOTHit* B4ROOTEvent::AddHit(B4ROOTHit& cand) {
								//cout << "Adding Pi0 candidate after " << m_NPi0Candidates << endl;
								TClonesArray &aCand = *m_Hits;
								B4ROOTHit *piCand = new B4ROOTHit(cand);
								aCand[m_NHits++]=piCand;
								//std::cout<<"AddHitCand: "<<cand.PhotNr()<<" AddHitPiCand: "<<piCand->PhotNr()<<std::endl;
								return piCand;
}


B4ROOTHit::B4ROOTHit() : TObject() {
								m_X = 0;
								m_Y = 0;
								m_Z = 0;
								m_EnergyDeposit = 0;
								m_calorimeterPart.clear();

								m_radiusWeight=1;




}

B4ROOTHit::B4ROOTHit(const B4ROOTHit& orig) : TObject(orig) {
								m_X = orig.m_X;
								m_Y = orig.m_Y;
								m_Z = orig.m_Z;

								m_EnergyDeposit = orig.m_EnergyDeposit;
								m_calorimeterPart=orig.m_calorimeterPart;
								m_radiusWeight=orig.m_radiusWeight;

}

B4ROOTHit::~B4ROOTHit(){
}
