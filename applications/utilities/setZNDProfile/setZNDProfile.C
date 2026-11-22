/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | CDMGpuFoam: Chemical-Diffusive Model solver
   \\    /   O peration     | for OpenFOAM-13 (openfoam.org)
    \\  /    A nd           |
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of CDMGpuFoam, derived from OpenFOAM.
    GNU General Public License v3 or later; see COPYING.

Application
    setZNDProfile

Description
    Maps a 1D ZND profile (written by tools/CDM_ZND.py --write) onto the
    mesh to initialise p, T, U and Y.

    system/setZNDProfileDict:
    \verbatim
        file            "znd.dat";      // relative to the case directory
        direction       (1 0 0);        // propagation direction
        shockPosition   0.01;           // coordinate of the leading shock
        frame           lab;            // lab | shock

        // Optional sinusoidal perturbation of the shock position,
        // used to trigger transverse waves in 2D/3D
        perturbation
        {
            tangent     (0 1 0);
            amplitude   2e-4;           // m
            wavelength  5e-3;           // m
        }
    \endverbatim

    frame lab:   unburnt gas at rest, detonation moving along +direction.
    frame shock: detonation at rest, unburnt gas entering at -D_CJ*direction
                 (the configuration of Lu, Kaplan & Oran 2021, Fig. 3).

    Behind the end of the tabulated profile the final (CJ) state is used.

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "Time.H"
#include "fvMesh.H"
#include "volFields.H"
#include "IFstream.H"
#include "IStringStream.H"
#include "mathematicalConstants.H"

#include <algorithm>

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"

    IOdictionary dict
    (
        IOobject
        (
            "setZNDProfileDict",
            runTime.system(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

    fileName file(dict.lookup<fileName>("file"));
    file.expand();
    if (!file.isAbsolute())
    {
        file = runTime.path()/file;
    }

    vector dir(dict.lookup<vector>("direction"));
    dir /= mag(dir);

    const scalar xShock0 = dict.lookup<scalar>("shockPosition");

    const word frame(dict.lookupOrDefault<word>("frame", "lab"));
    if (frame != "lab" && frame != "shock")
    {
        FatalIOErrorInFunction(dict)
            << "frame must be lab or shock" << exit(FatalIOError);
    }

    const bool perturb = dict.isDict("perturbation");
    vector tangent(Zero);
    scalar amplitude = 0;
    scalar wavelength = 1;
    if (perturb)
    {
        const dictionary& pDict = dict.subDict("perturbation");
        tangent = pDict.lookup<vector>("tangent");
        tangent /= mag(tangent);
        amplitude = pDict.lookup<scalar>("amplitude");
        wavelength = pDict.lookup<scalar>("wavelength");
    }


    // Read the table
    // ~~~~~~~~~~~~~~

    Info<< "Reading ZND profile " << file << endl;

    IFstream is(file);
    if (!is.good())
    {
        FatalErrorInFunction
            << "Cannot open " << file << exit(FatalError);
    }

    scalar p0 = -1, T0 = -1, D = -1;
    DynamicList<scalar> s, pT, TT, uT, YT;

    string line;
    while (is.good() && !is.eof())
    {
        is.getLine(line);

        if (line.empty())
        {
            continue;
        }

        if (line[0] == '#')
        {
            // Header: "# p0 <v> T0 <v> D <v> xd <v>"
            if (line.find("p0") != string::npos)
            {
                IStringStream hs(line.substr(1));
                word key;
                scalar value;
                while (hs.good())
                {
                    hs >> key;
                    if (!hs.good()) break;
                    hs >> value;
                    if (key == "p0") p0 = value;
                    else if (key == "T0") T0 = value;
                    else if (key == "D") D = value;
                }
            }
            continue;
        }

        IStringStream ls(line);
        scalar si, pi, Ti, ui, Yi, rhoi;
        ls >> si >> pi >> Ti >> ui >> Yi >> rhoi;

        s.append(si);
        pT.append(pi);
        TT.append(Ti);
        uT.append(ui);
        YT.append(Yi);
    }

    if (p0 < 0 || T0 < 0 || D < 0 || s.size() < 2)
    {
        FatalErrorInFunction
            << "Invalid ZND table " << file
            << " (missing header or data)" << exit(FatalError);
    }

    Info<< "    " << s.size() << " points, p0 = " << p0 << ", T0 = " << T0
        << ", D = " << D << ", length = " << s.last() << nl << endl;


    // Read the fields
    // ~~~~~~~~~~~~~~~

    volScalarField p
    (
        IOobject("p", runTime.name(), mesh, IOobject::MUST_READ),
        mesh
    );
    volScalarField T
    (
        IOobject("T", runTime.name(), mesh, IOobject::MUST_READ),
        mesh
    );
    volScalarField Y
    (
        IOobject("Y", runTime.name(), mesh, IOobject::MUST_READ),
        mesh
    );
    volVectorField U
    (
        IOobject("U", runTime.name(), mesh, IOobject::MUST_READ),
        mesh
    );


    // Map
    // ~~~

    const vectorField& C = mesh.C().primitiveField();
    const scalar twoPi = constant::mathematical::twoPi;

    label nReacting = 0;

    forAll(C, celli)
    {
        scalar xs = xShock0;
        if (perturb)
        {
            xs += amplitude*Foam::sin(twoPi*(C[celli] & tangent)/wavelength);
        }

        const scalar dist = xs - (C[celli] & dir);

        scalar pi, Ti, ui, Yi;

        if (dist < 0)
        {
            // Unshocked gas
            pi = p0;
            Ti = T0;
            ui = D;
            Yi = 1;
        }
        else if (dist >= s.last())
        {
            pi = pT.last();
            Ti = TT.last();
            ui = uT.last();
            Yi = YT.last();
            nReacting++;
        }
        else
        {
            const label j =
                std::upper_bound(s.begin(), s.end(), dist) - s.begin();
            const label i0 = max(j - 1, 0);
            const label i1 = min(j, s.size() - 1);
            const scalar w =
                i1 == i0 ? 0 : (dist - s[i0])/(s[i1] - s[i0]);

            pi = (1 - w)*pT[i0] + w*pT[i1];
            Ti = (1 - w)*TT[i0] + w*TT[i1];
            ui = (1 - w)*uT[i0] + w*uT[i1];
            Yi = (1 - w)*YT[i0] + w*YT[i1];
            nReacting++;
        }

        p[celli] = pi;
        T[celli] = Ti;
        Y[celli] = Yi;

        // u is the speed relative to the shock, directed towards -dir
        U[celli] = frame == "lab" ? (D - ui)*dir : -ui*dir;
    }

    Info<< "Set " << returnReduce(nReacting, sumOp<label>())
        << " cells behind the shock out of "
        << returnReduce(C.size(), sumOp<label>()) << nl << endl;

    p.write();
    T.write();
    Y.write();
    U.write();

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
