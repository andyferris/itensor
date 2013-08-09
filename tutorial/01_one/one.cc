#include "core.h"

using namespace std;
using boost::format;

int
main(int argc, char* argv[])
    {

    //
    // Single-site wavefunction
    //
    
    Index s("s",2);

    ITensor psi(s); //default initialized to zero


    //
    // Initialize up spin
    //

    psi(s(1)) = 1;

    PrintDat(psi);
    
    //exit(0); //uncomment to exit here

    //
    // Operators 
    //

    ITensor Sz(s,primed(s)),
            Sx(s,primed(s));

    commaInit(Sz,s,primed(s)) << 0.5, 0.0,
                                 0.0,-0.5;

    commaInit(Sx,s,primed(s)) << 0.0, 0.5,
                                 0.5, 0.0;

    PrintDat(Sz);
    PrintDat(Sx);

    //exit(0); //uncomment to exit here

    //
    // Product Sx * phi 
    //

    ITensor phi = Sx * psi;

    phi.noprime();

    PrintDat(phi);

    //exit(0); //uncomment to exit here

    //
    // 45* angle spin
    //

    const Real theta = Pi/4;

    //Extra factors of two come from S=1/2 representation
    psi(s(1)) = cos(theta/2.);
    psi(s(2)) = sin(theta/2.);

    PrintDat(psi);

    //exit(0); //uncomment to exit here

    //
    // Expectation values
    //

    ITensor cpsi = conj(primed(psi));

    Real zz = (cpsi * Sz * psi).toReal();
    Real xx = (cpsi * Sx * psi).toReal();

    cout << format("<Sz> = %.5f") % zz << endl;
    cout << format("<Sx> = %.5f") % xx << endl;

    return 0;
    }
