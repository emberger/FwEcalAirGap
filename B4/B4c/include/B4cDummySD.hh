#include "G4VSensitiveDetector.hh"

#include "B4cCalorHit.hh"

#include <vector>

class G4Step;
class G4HCofThisEvent;


class dummySD : public G4VSensitiveDetector{
	
public:
	
	dummySD(const G4String name);
	~dummySD();
	
	
	
	

};
