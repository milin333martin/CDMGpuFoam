#!/usr/bin/env python3
"""
writeSourceCase.py: generate the wave-fixed (shock-frame) detonation case
that later initialises an arrestor or narrow-channel run.

Fresh gas enters the right-hand boundary at -D_CJ, so the detonation stands
still in the domain and only the cellular structure evolves.  A lab-frame
soot foil is reconstructed in situ (Ramachandran, Dushe & Yang, AIAA
2025-1386, Eq. 13), so the cells can be measured even though the front does
not move.

Driven by runSource.sh; every setting arrives as a command-line option.
"""

import argparse
import os
import shutil

HEAD = """/*--------------------------------*- C++ -*----------------------------------*\\
  =========                 |
  \\\\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\\\    /   O peration     | Version:  13
    \\\\  /    A nd           |
     \\\\/     M anipulation  |
\\*---------------------------------------------------------------------------*/
FoamFile
{
    format      ascii;
    class       %s;
    location    "%s";
    object      %s;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
"""


def w(case, path, cls, obj, body, location):
    p = os.path.join(case, path)
    os.makedirs(os.path.dirname(p), exist_ok=True)
    open(p, "w").write(HEAD % (cls, location, obj) + body.strip() + "\n\n"
                       "// *********************************************"
                       "******************************** //\n")


def main():
    a = argparse.ArgumentParser()
    g = a.add_argument
    g("--case", required=True)
    # geometry, mm
    g("--H", type=float, default=45.0, help="channel height")
    g("--L", type=float, default=90.0, help="domain length")
    g("--thickness", type=float, default=0.2)
    g("--shockPosition", type=float, default=70.0,
      help="where the front is placed, mm from the outlet")
    # mesh
    g("--dx", type=float, default=0.2, help="base cell, mm")
    g("--levels", type=int, default=2, help="AMR levels (2 -> 50 um at 0.2 mm)")
    g("--maxCells", type=int, default=4000000)
    # perturbation that seeds the cells
    g("--nWaves", type=float, default=3.0,
      help="transverse waves across the channel at t = 0")
    g("--amplitude", type=float, default=0.25,
      help="perturbation amplitude in units of x_d")
    # model
    g("--DCJ", type=float, required=True, help="m/s, from CDM_ZND.py")
    g("--xd", type=float, default=0.626, help="half-reaction length, mm")
    g("--A", default="1.01e12")
    g("--EaByRT0", default="39.98")
    g("--qMwByRT0", default="62.97")
    g("--kappa0", default="2.00e-5")
    g("--T0", default="298")
    g("--Mw", default="22.77")
    g("--Cv", default="2921.1990")
    g("--p0", default="101325")
    # run
    g("--endTime", type=float, default=1e-4)
    g("--writeInterval", type=float, default=5e-6)
    g("--maxCo", type=float, default=0.4)
    g("--fluxScheme", default="HLLC")
    g("--nProcs", type=int, default=8)
    g("--decompose", default="8 1")
    g("--distributor", default="zoltan")
    g("--fileHandler", default="uncollated")
    g("--foilDelta", type=float, default=0.05, help="soot foil pixel, mm")
    args = a.parse_args()

    mm = 1e-3
    H, L, th = args.H*mm, args.L*mm, args.thickness*mm
    case = args.case

    nx = max(1, int(round(L/(args.dx*mm))))
    ny = max(1, int(round(H/(args.dx*mm))))
    finest = args.dx*1000/2**args.levels          # um

    # the foil must cover the domain plus everything the front travels over
    foilMax = L + args.DCJ*args.endTime

    if os.path.isdir(case):
        shutil.rmtree(case)
    os.makedirs(case)

    w(case, "system/blockMeshDict", "dictionary", "blockMeshDict", """
// Wave-fixed (shock-frame) channel: fresh gas enters at x = L with -D_CJ,
// so the detonation stands still and only its structure evolves.
// %.4g x %.4g mm, base cell %.4g mm, %d AMR levels -> %.4g um at the front
// (x_d = %.4g mm, so the reaction zone has %.1f cells).
convertToMeters 1;

vertices
(
    (0 0 0) (%.9g 0 0) (%.9g %.9g 0) (0 %.9g 0)
    (0 0 %.9g) (%.9g 0 %.9g) (%.9g %.9g %.9g) (0 %.9g %.9g)
);

blocks
(
    hex (0 1 2 3 4 5 6 7) (%d %d 1) simpleGrading (1 1 1)
);

boundary
(
    internalFaces { type internal; faces (); }

    outlet  { type patch; faces ((0 4 7 3)); }
    inlet   { type patch; faces ((1 2 6 5)); }
    bottom  { type symmetryPlane; faces ((0 1 5 4)); }
    top     { type symmetryPlane; faces ((3 7 6 2)); }
    frontAndBack { type empty; faces ((0 3 2 1) (4 5 6 7)); }
);""" % (args.L, args.H, args.dx, args.levels, finest, args.xd,
         args.xd/(finest/1000),
         L, L, H, H, th, L, th, L, H, th, H, th, nx, ny), location="system")

    w(case, "constant/physicalProperties", "dictionary", "physicalProperties",
      """
thermoType
{
    type            hePsiThermo;
    mixture         pureMixture;
    transport       const;
    thermo          eConst;
    equationOfState perfectGas;
    specie          specie;
    energy          sensibleInternalEnergy;
}

mixture
{
    specie      { molWeight %s; }
    thermodynamics { Cv %s; hf 0; }
    transport   { mu 0; Pr 1; }
}""" % (args.Mw, args.Cv), location="constant")

    w(case, "constant/CDMProperties", "dictionary", "CDMProperties", """
// Chemical-diffusive model (Lu, Kaplan & Oran, Combust. Flame 230 (2021))
A_cgs       %s;
EaByRT0     %s;
qMwByRT0    %s;
kappa0_cgs  %s;
T0          %s;

n           0.7;
Pr          1;
Le          1;
cpuLoad     yes;""" % (args.A, args.EaByRT0, args.qMwByRT0, args.kappa0,
                       args.T0), location="constant")

    w(case, "constant/dynamicMeshDict", "dictionary", "dynamicMeshDict", """
topoChanger
{
    type            refiner2D;
    libs            ("libCDMDynamicMesh2D.so");

    refineInterval  2;

    field               refinementIndicator;
    lowerRefineLevel    0.05;
    upperRefineLevel    1e6;
    unrefineLevel       0.015;

    nBufferLayers   2;
    maxRefinement   %d;
    maxCells        %d;

    correctFluxes   ((phi none));
    dumpLevel       true;
}

distributor
{
    type            loadBalancer;
    libs            ("libfvMeshDistributors.so");

    redistributionInterval  20;
    maxImbalance    0.1;
}""" % (args.levels, args.maxCells), location="constant")

    fh = ("\nOptimisationSwitches\n{\n    fileHandler     collated;\n"
          "    maxThreadFileBufferSize 2e9;\n}\n"
          if args.fileHandler == "collated" else "")

    w(case, "system/controlDict", "dictionary", "controlDict", """
application     foamRun;
solver          CDMGpuFoam;

startFrom       latestTime;
startTime       0;
stopAt          endTime;
endTime         %g;

deltaT          1e-10;
writeControl    adjustableRunTime;
writeInterval   %g;

purgeWrite      0;
writeFormat     binary;
writePrecision  10;
writeCompression off;
timeFormat      general;
timePrecision   8;
runTimeModifiable yes;

adjustTimeStep  yes;
maxCo           %g;
maxDeltaT       1e-7;
%s""" % (args.endTime, args.writeInterval, args.maxCo, fh), location="system")

    w(case, "system/fvSchemes", "dictionary", "fvSchemes", """
fluxScheme      %s;         // HLLC | HLLCLM | HLLCP (experimental)

ddtSchemes      { default none; }
gradSchemes     { default Gauss linear; }
divSchemes      { default none; }
laplacianSchemes { default Gauss linear corrected; }
interpolationSchemes
{
    default         linear;
    reconstruct(rho) vanLeer;
    reconstruct(U)  vanLeerV;
    reconstruct(p)  vanLeer;
    reconstruct(Y)  vanLeer;
}
snGradSchemes   { default corrected; }""" % args.fluxScheme,
      location="system")

    w(case, "system/fvSolution", "dictionary", "fvSolution", """
solvers {}

PIMPLE
{
    nRKStages           2;
    viscous             yes;
    chemistry           yes;
    chemistrySplitting  Strang;

    // The front should hold station in this frame: x stays put and the
    // speed hovers about zero once the detonation is quasi-steady.
    frontTracking
    {
        direction       (1 0 0);
        pThreshold      5e5;
        YFlame          0.5;
        window          20;
        interval        20;
    }

    // Lab-frame soot foil reconstructed in situ from the shock-frame run
    // (Ramachandran, Dushe & Yang, AIAA 2025-1386, Eq. 13)
    sootFoil
    {
        frameVelocity   (%g 0 0);
        min             (0 0);
        max             (%.9g %.9g);
        delta           %.9g;
    }

    refinementIndicator
    {
        TReacting       500;
        YReacting       1e-3;
    }
}""" % (args.DCJ, foilMax, H, args.foilDelta*mm), location="system")

    nxp, nyp = args.decompose.split()
    w(case, "system/decomposeParDict", "dictionary", "decomposeParDict", """
numberOfSubdomains %d;

decomposer      hierarchical;
hierarchicalCoeffs
{
    n           (%s %s 1);
    order       xyz;
}

distributor     %s;
zoltanCoeffs
{
    lb_method   graph;
    lb_approach repartition;
}

libs            ("libCDMDynamicMesh2D.so"%s);

constraints
{
    refinementHistory
    {
        type        refinementHistory2D;
    }
}""" % (args.nProcs, nxp, nyp, args.distributor,
        ' "libzoltanDecomp.so"' if args.distributor == "zoltan" else
        (' "libptscotchDecomp.so"' if args.distributor == "ptscotch" else "")),
      location="system")

    w(case, "system/setZNDProfileDict", "dictionary", "setZNDProfileDict", """
file            "znd.dat";
direction       (1 0 0);
shockPosition   %.9g;
frame           shock;

// %.4g transverse waves across the channel seed the cellular instability
perturbation
{
    tangent         (0 1 0);
    amplitude       %.9g;
    wavelength      %.9g;
}""" % (args.shockPosition*mm, args.nWaves, args.amplitude*args.xd*mm,
        H/args.nWaves), location="system")

    def field(name, cls, dims, internal, inlet):
        w(case, "0/" + name, cls, name, """
dimensions      %s;

internalField   uniform %s;

boundaryField
{
    internalFaces   { type internal; }
    outlet          { type zeroGradient; }
    bottom          { type symmetryPlane; }
    top             { type symmetryPlane; }
    frontAndBack    { type empty; }
    inlet
    {
        %s
    }
}""" % (dims, internal, inlet), location="0")

    def fixed(v):
        return "type            fixedValue;\n        value           uniform %s;" % v

    field("p", "volScalarField", "[1 -1 -2 0 0 0 0]", args.p0,
          "type            zeroGradient;")
    field("T", "volScalarField", "[0 0 0 1 0 0 0]", args.T0, fixed(args.T0))
    field("Y", "volScalarField", "[0 0 0 0 0 0 0]", "1", fixed("1"))
    field("U", "volVectorField", "[0 1 -1 0 0 0 0]",
          "(%g 0 0)" % -args.DCJ, fixed("(%g 0 0)" % -args.DCJ))

    print("Wave-fixed source case: %s" % case)
    print("  domain     : %.4g x %.4g mm, front placed at x = %.4g mm"
          % (args.L, args.H, args.shockPosition))
    print("  mesh       : %d x %d = %d base cells of %.4g mm, %d levels "
          "-> %.4g um" % (nx, ny, nx*ny, args.dx, args.levels, finest))
    print("  resolution : %.1f cells per half-reaction length (x_d = %.4g mm)"
          % (args.xd/(finest/1000), args.xd))
    print("  frame      : fresh gas enters at -%.6g m/s; the foil covers "
          "%.4g mm of travel" % (args.DCJ, foilMax/mm))
    if args.xd/(finest/1000) < 10:
        print("  WARNING: fewer than 10 cells per x_d; Lu et al. (Appendix A) "
              "found the cell size converged from x_d/12")


if __name__ == "__main__":
    main()
