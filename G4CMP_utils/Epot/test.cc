#include <string>
#include "Epot.hh"
#include "G4ElectricField.hh"
#include <iostream>
#include <fstream>

int main()
{
    G4String file("Epot_iZip4_test");

    G4ElectricField* myField = new DetectorField(file);

    G4double Point1[4] = {.00328*m,.0173*m,.0030*m};
    G4double Efield1[3] = {0,0,0};
    G4double Point2[4] = {0,0,1.27*cm,0};
    G4double Efield2[3] = {0,0,0};
    G4double Point3[4] = {-.00328*m,.0246*m,.0133*m,0};
    G4double Efield3[3] = {0,0,0};

    myField->GetFieldValue(Point1,Efield1);

    G4cout << "Field at point1 in xtal: " << Efield1[0]/volt*m << " " 
					     << Efield1[1]/volt*m << " "
					     << Efield1[2]/volt*m << G4endl;
    
    myField->GetFieldValue(Point2,Efield2);

    G4cout << "Field at point2 in xtal: " << Efield2[0]/volt*m << " "
					     << Efield2[1]/volt*m<< " "
					     << Efield2[2]/volt*m << G4endl;

    myField->GetFieldValue(Point3,Efield3);

    G4cout << "Field at point3 in xtal: " << Efield3[0]/volt*m << " "
					     << Efield3[1]/volt*m << " "
					     << Efield3[2]/volt*m << G4endl;
    return 0;
}
