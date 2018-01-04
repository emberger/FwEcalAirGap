#include "G4VReadOutGeometry.hh"



#include "globals.hh"






class MyRO : public G4VReadOutGeometry {
	
	private:
		

	
	
	
	
	public:

		MyRO();
		MyRO(G4String aString);
		virtual ~MyRO();



	protected:

		virtual G4VPhysicalVolume* Build();
		
		
	
	
};
