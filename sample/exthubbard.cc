#include "core.h"
#include "model/spinhalf.h"
#include "hams/ExtendedHubbard.h"
using boost::format;
using namespace std;

int main(int argc, char* argv[])
    {
    //Parse the input file
    if(argc != 2)
        {
        cout << "Usage: " << argv[0] << " inputfile." << endl;
        return 0;
        }
    string infilename(argv[1]);
    InputFile infile(infilename);
    InputGroup basic(infile,"basic");

    int N = 0;
    basic.GetIntM("N",N); //the 'M' stands for mandatory
    int Npart = N; //number of particles, default is N (half filling)
    basic.GetInt("Npart",Npart);

    int nsweeps = 0;
    basic.GetIntM("nsweeps",nsweeps);
    Real t1 = 0;
    basic.GetRealM("t1",t1);
    Real t2 = 0;
    basic.GetRealM("t2",t2);
    Real U = 0;
    basic.GetRealM("U",U);
    Real V1 = 0;
    basic.GetRealM("V1",V1);

    InputGroup table(basic,"sweeps");
    Sweeps sweeps(nsweeps,table);
    cout << sweeps;

    //
    // Initialize the site degrees of freedom.
    //
    Hubbard model(N);

    //
    // Create the Hamiltonian matrix product operator.
    // Here we use the IQMPO class which is an MPO of 
    // IQTensors, tensors whose indices are sorted
    // with respect to quantum numbers
    //
    IQMPO H = ExtendedHubbard(model,
                              Opt("U",U)
                              & Opt("t1",t1)
                              & Opt("t2",t2)
                              & Opt("V1",V1));

    //
    // Set the initial wavefunction matrix product state
    // to be a Neel state.
    //
    InitState initState(N);
    int p = Npart;
    for(int i = N; i >= 1; --i) 
        {
        if(p > i)
            {
            cout << "Doubly occupying site " << i << endl;
            initState(i) = model.UpDn(i);
            p -= 2;
            }
        else
        if(p > 0)
            {
            cout << "Singly occupying site " << i << endl;
            initState(i) = (i%2==1 ? model.Up(i) : model.Dn(i));
            p -= 1;
            }
        else
            {
            initState(i) = model.Emp(i);
            }
        }

    IQMPS psi(model,initState);

    cout << totalQN(psi) << endl;

    //
    // Begin the DMRG calculation
    //
    Real En = dmrg(psi,H,sweeps,Quiet());

    //
    // Measure spin densities
    //
    Vector upd(N),dnd(N);
    for(int j = 1; j <= N; ++j)
        {
        psi.position(j);
        upd(j) = Dot(conj(primed(psi.AA(j),Site)),model.Nup(j)*psi.AA(j));
        dnd(j) = Dot(conj(primed(psi.AA(j),Site)),model.Ndn(j)*psi.AA(j));
        }

    cout << "Up Density:" << endl;
    for(int j = 1; j <= N; ++j)
        cout << format("%d %.10f\n") % j % upd(j);
    cout << endl;

    cout << "Dn Density:" << endl;
    for(int j = 1; j <= N; ++j)
        cout << format("%d %.10f\n") % j % dnd(j);
    cout << endl;

    cout << "Total Density:" << endl;
    for(int j = 1; j <= N; ++j)
        cout << format("%d %.10f\n") % j % (upd(j)+dnd(j));
    cout << endl;

    //
    // Print the final energy reported by DMRG
    //
    cout << format("\nGround State Energy = %.10f\n")%En;

    return 0;
    }
