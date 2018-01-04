//---------------------------------------------------------------
//
// G4UserTrackingAction.hh
//
// Description:
// This class represents actions taken place by the user at
// the start/end point of processing one track.
//
//---------------------------------------------------------------

#include "G4TrackingManager.hh"
#include "B4cTrackInformation.hh"
///////////////////////////
class B4cTrackingAction : public G4UserTrackingAction
///////////////////////////
{

//--------
public:
//--------

// Constructor & Destructor
B4cTrackingAction();
virtual ~B4cTrackingAction();

// Member functions
virtual void PreUserTrackingAction(const G4Track*);

virtual void PostUserTrackingAction(const G4Track*);

//-----------
protected:
//-----------

private:

G4bool isdecayed=false;
G4bool isPion=false;

B4cTrackInformation * theInfo;

// Member data
//G4TrackingManager* fpTrackingManager;



};
