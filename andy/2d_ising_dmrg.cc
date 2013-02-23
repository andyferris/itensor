#include "core.h"
#include "model/spinhalf.h"
#include "model/spinone.h"
#include "hams/ising.h"
using boost::format;
using namespace std;

int main(int argc, char* argv[])
{
    //Parse the input file
    if(argc != 2)
    {
        cout << "Usage: " << argv[0] << " inputfile_dmrg_table" << endl;
        return 0;
    }
    string infilename(argv[1]);
    InputFile infile(infilename);
    InputGroup basic(infile,"basic");
    
    //Read in individual parameters from the input file
    int Nx = 0;
    int Ny = 0;
    basic.GetIntM("Nx",Nx); //the 'M' stands for mandatory
    basic.GetIntM("Ny",Ny);
    int nsweeps = 0;
    basic.GetIntM("nsweeps",nsweeps);
    int quiet = 1;
    basic.GetYesNo("quiet",quiet);
    int cylinder = 0;
    basic.GetYesNo("cylinder",cylinder);
    Real h = 0;
    basic.GetRealM("h",h);
    
    
    // Read the sweeps parameters from a table.
    
    //Read in the sweeps table itself
    InputGroup table(basic,"sweeps");
    
    //Create the sweeps class & print
    Sweeps sweeps(nsweeps,table);
    cout << sweeps;
    
    // Derived quantities
    int N = Nx*Ny;
    
    //
    // Initialize the site degrees of freedom.
    //
    //SpinOne model(N);    // make a chain of N spin 1's
    SpinHalf model(N); // make a chain of N spin 1/2's

    //
    // Create the Hamiltonian matrix product operator (MPO)
    //
    //double h = 3.05; // Magnetic field
    MPO H = Ising(model,h,Ny,cylinder); // Not cylindrical b.c.

    //
    // Set the initial wavefunction matrix product state (MPS)
    // to be a Neel state.
    //
    InitState initState(N);
    for(int i = 1; i <= N; ++i) 
        initState(i) = (i%2==1 ? model.Up(i) : model.Dn(i));

    MPS psi(model,initState);

    //
    // psiHphi calculates matrix elements of MPO's with respect to MPS's
    // psiHphi(psi,H,psi) = <psi|H|psi>
    //
    cout << format("Initial energy = %.5f") % psiHphi(psi,H,psi) << endl;

    //
    // Set the parameters controlling the accuracy of the DMRG
    // calculation for each DMRG sweep. Here less than 20 maxm
    // values are provided, so all remaining sweeps will use the
    // last maxm (= 500).
    //
    //Sweeps sweeps(20);
    //sweeps.maxm() = 50,50,100,100,200,200,300,300,400,400,500;
    //sweeps.cutoff() = 1E-10;
    //cout << sweeps;

    //
    // Begin the DMRG calculation
    //
    Real En = dmrg(psi,H,sweeps,Quiet());

    //
    // Print the final energy reported by DMRG
    //
    cout << format("\nGround State Energy = %.10f")%En << endl;
    
    return 0;
}
