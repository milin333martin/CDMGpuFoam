#!/usr/bin/env python3
"""
writeCHTTest.py: a minimal two-region conjugate heat transfer test for
CDMGpuFoam.

Hot, quiescent burnt gas in a channel sits on a stainless steel plate.  No
flow, no reaction, no AMR: just conduction across the interface, which is the
one thing that has to work before CHT is put into an arrestor case.

Both bodies are thermally semi-infinite over the test duration, so the
interface temperature has a closed-form answer,

    Ti = (e_f Tf + e_s Ts)/(e_f + e_s),     e = sqrt(kappa rho cp)

which the script prints for comparison with the result.

Stages:
    --stage base    write the mesh and dictionaries (before splitMeshRegions)
    --stage fields  write 0/fluid and 0/solid using the split patch names
"""

import argparse
import math
import os
import re
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
    %sobject      %s;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
"""


def w(case, path, cls, obj, body, location=None):
    p = os.path.join(case, path)
    os.makedirs(os.path.dirname(p), exist_ok=True)
    loc = 'location    "%s";\n    ' % location if location else ""
    open(p, "w").write(HEAD % (cls, loc, obj) + body.strip() + "\n\n"
                       "// *********************************************"
                       "******************************** //\n")


def mapped_patch(case, region):
    """name of the region's mapped (interface) patch after splitMeshRegions"""
    b = os.path.join(case, "constant", region, "polyMesh", "boundary")
    txt = open(b).read()
    names = re.findall(r"^\s{4}(\w+)\s*$\s*\{([^}]*)\}", txt, re.M)
    for name, body in names:
        if "mapped" in body:
            return name
    raise SystemExit("no mapped patch found in %s" % b)


def main():
    a = argparse.ArgumentParser()
    a.add_argument("--case", default="chtTest")
    a.add_argument("--stage", default="base", choices=("base", "fields"))
    # geometry, mm
    a.add_argument("--L", type=float, default=4.0, help="length")
    a.add_argument("--hFluid", type=float, default=2.0)
    a.add_argument("--hSolid", type=float, default=1.0)
    a.add_argument("--dx", type=float, default=0.02, help="cell size")
    # states
    a.add_argument("--Tgas", type=float, default=2800.0)
    a.add_argument("--Tsolid", type=float, default=300.0)
    a.add_argument("--p0", type=float, default=1772635.0)
    # steel
    a.add_argument("--rhoS", type=float, default=7900.0)
    a.add_argument("--cvS", type=float, default=500.0)
    a.add_argument("--kappaS", type=float, default=16.0)
    # CDM (Table 6)
    a.add_argument("--kappa0", default="2.00e-5")
    a.add_argument("--Cv", type=float, default=2921.199)
    a.add_argument("--Mw", type=float, default=22.77)
    a.add_argument("--endTime", type=float, default=2e-5)
    args = a.parse_args()

    case, mm = args.case, 1e-3
    L, hf, hs = args.L*mm, args.hFluid*mm, args.hSolid*mm
    nx = max(1, int(round(L/(args.dx*mm))))
    nyf = max(1, int(round(hf/(args.dx*mm))))
    nys = max(1, int(round(hs/(args.dx*mm))))
    th = args.dx*mm

    if args.stage == "fields":
        fp = mapped_patch(case, "fluid")
        sp = mapped_patch(case, "solid")

        def field(region, name, cls, dims, internal, patches):
            w(case, "0/%s/%s" % (region, name), cls, name, """
dimensions      %s;

internalField   uniform %s;

boundaryField
{
%s
    frontAndBack    { type empty; }

    #includeEtc "caseDicts/setConstraintTypes"
}""" % (dims, internal, patches), location="0/" + region)

        wallsZG = """    walls           { type zeroGradient; }
    ends            { type zeroGradient; }
"""
        field("fluid", "p", "volScalarField", "[1 -1 -2 0 0 0 0]", args.p0,
              wallsZG + "    %s { type zeroGradient; }\n" % fp)
        field("fluid", "T", "volScalarField", "[0 0 0 1 0 0 0]", args.Tgas,
              wallsZG + """    %s
    {
        type            CDMCoupledTemperature;
        value           uniform %g;
    }
""" % (fp, args.Tgas))
        field("fluid", "Y", "volScalarField", "[0 0 0 0 0 0 0]", "0",
              wallsZG + "    %s { type zeroGradient; }\n" % fp)
        field("fluid", "U", "volVectorField", "[0 1 -1 0 0 0 0]", "(0 0 0)",
              """    walls           { type slip; }
    ends            { type slip; }
    %s              { type noSlip; }
""" % fp)
        field("solid", "T", "volScalarField", "[0 0 0 1 0 0 0]", args.Tsolid,
              """    walls           { type zeroGradient; }
    ends            { type zeroGradient; }
    %s
    {
        type            coupledTemperature;
        value           uniform %g;
    }
""" % (sp, args.Tsolid))

        # analytic contact temperature
        R = 8314.46/args.Mw
        cp = args.Cv + R
        # K = kappa0 cp T^0.7 with kappa0 in cgs g/(s cm K^n) -> SI kg/(m s K^n)
        kappa0 = float(args.kappa0)*0.1
        kf = kappa0*cp*args.Tgas**0.7
        rhof = args.p0/(R*args.Tgas)
        ef = math.sqrt(kf*rhof*cp)
        es = math.sqrt(args.kappaS*args.rhoS*args.cvS)
        Ti = (ef*args.Tgas + es*args.Tsolid)/(ef + es)
        alphaS = args.kappaS/(args.rhoS*args.cvS)
        print("  fluid at interface : k = %.4g W/m/K, rho = %.4g kg/m3, "
              "cp = %.5g J/kg/K" % (kf, rhof, cp))
        print("  effusivity         : gas %.4g, steel %.4g J/(m2 K s^0.5)"
              % (ef, es))
        print("  EXPECTED interface temperature: %.2f K" % Ti)
        print("  steel penetration depth at %.3g s: %.2f um"
              % (args.endTime, 1e6*math.sqrt(alphaS*args.endTime)))
        return

    # ---------------- base case ----------------
    if os.path.isdir(case):
        shutil.rmtree(case)
    os.makedirs(case)

    ys = [-hs, 0.0, hf]
    verts = []
    for z in (0.0, th):
        for y in ys:
            for x in (0.0, L):
                verts.append("    (%.9g %.9g %.9g)" % (x, y, z))

    def V(ix, iy, iz):
        return iz*6 + iy*2 + ix

    def hexb(iy, zone, ny):
        v = [V(0, iy, 0), V(1, iy, 0), V(1, iy+1, 0), V(0, iy+1, 0),
             V(0, iy, 1), V(1, iy, 1), V(1, iy+1, 1), V(0, iy+1, 1)]
        return ("    hex (%s) %s (%d %d 1) simpleGrading (1 1 1)"
                % (" ".join(str(i) for i in v), zone, nx, ny))

    def face(iy, side):
        q = {
            "x0": [V(0, iy, 0), V(0, iy, 1), V(0, iy+1, 1), V(0, iy+1, 0)],
            "x1": [V(1, iy, 0), V(1, iy+1, 0), V(1, iy+1, 1), V(1, iy, 1)],
            "y0": [V(0, iy, 0), V(1, iy, 0), V(1, iy, 1), V(0, iy, 1)],
            "y1": [V(0, iy+1, 0), V(0, iy+1, 1), V(1, iy+1, 1), V(1, iy+1, 0)],
            "z0": [V(0, iy, 0), V(0, iy+1, 0), V(1, iy+1, 0), V(1, iy, 0)],
            "z1": [V(0, iy, 1), V(1, iy, 1), V(1, iy+1, 1), V(0, iy+1, 1)],
        }[side]
        return "        (%s)" % " ".join(str(i) for i in q)

    w(case, "system/blockMeshDict", "dictionary", "blockMeshDict", """
// CHT test: %.4g mm of burnt gas on %.4g mm of steel, cells %.4g mm
convertToMeters 1;

vertices
(
%s
);

blocks
(
%s
%s
);

boundary
(
    walls
    {
        type wall;
        faces
        (
%s
%s
        );
    }
    ends
    {
        type wall;
        faces
        (
%s
        );
    }
    frontAndBack
    {
        type empty;
        faces
        (
%s
        );
    }
);""" % (args.hFluid, args.hSolid, args.dx, "\n".join(verts),
         hexb(0, "solid", nys), hexb(1, "fluid", nyf),
         face(0, "y0"), face(1, "y1"),
         "\n".join([face(0, "x0"), face(0, "x1"),
                    face(1, "x0"), face(1, "x1")]),
         "\n".join([face(0, "z0"), face(0, "z1"),
                    face(1, "z0"), face(1, "z1")])), location="system")

    w(case, "system/controlDict", "dictionary", "controlDict", """
application     foamMultiRun;

// the solid is constructed first: the fluid's coupled temperature condition
// needs the solid region to exist, while the solid reads its fields from
// file without evaluating coupled conditions
regionSolvers
{
    solid           solid;
    fluid           CDMGpuFoam;
}

startFrom       latestTime;
startTime       0;
stopAt          endTime;
endTime         %g;

deltaT          1e-10;
writeControl    adjustableRunTime;
writeInterval   %g;

purgeWrite      0;
writeFormat     ascii;
writePrecision  8;
writeCompression off;
timeFormat      general;
timePrecision   8;
runTimeModifiable yes;

adjustTimeStep  yes;
maxCo           0.4;
maxDeltaT       1e-8;""" % (args.endTime, args.endTime/4), location="system")

    # foamMultiRun reads the outer control from the case-level fvSolution
    w(case, "system/fvSolution", "dictionary", "fvSolution", """
PIMPLE
{
    nOuterCorrectors 1;
}""", location="system")
    w(case, "system/fvSchemes", "dictionary", "fvSchemes", """
ddtSchemes      { default none; }
gradSchemes     { default Gauss linear; }
divSchemes      { default none; }
laplacianSchemes { default Gauss linear corrected; }
interpolationSchemes { default linear; }
snGradSchemes   { default corrected; }""", location="system")

    # fluid region
    w(case, "system/fluid/fvSchemes", "dictionary", "fvSchemes", """
fluxScheme      HLLC;

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
snGradSchemes   { default corrected; }""", location="system/fluid")

    w(case, "system/fluid/fvSolution", "dictionary", "fvSolution", """
solvers {}

PIMPLE
{
    nRKStages           2;
    viscous             yes;
    chemistry           no;     // conduction only for this test
    CHT                 yes;    // publish K for CDMCoupledTemperature
}""", location="system/fluid")

    w(case, "constant/fluid/physicalProperties", "dictionary",
      "physicalProperties", """
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
    specie      { molWeight %g; }
    thermodynamics { Cv %g; hf 0; }
    transport   { mu 0; Pr 1; }
}""" % (args.Mw, args.Cv), location="constant/fluid")

    w(case, "constant/fluid/CDMProperties", "dictionary", "CDMProperties", """
A_cgs       1.01e12;
EaByRT0     39.98;
qMwByRT0    62.97;
kappa0_cgs  %s;
T0          298;

n           0.7;
Pr          1;
Le          1;""" % args.kappa0, location="constant/fluid")

    # solid region
    w(case, "system/solid/fvSchemes", "dictionary", "fvSchemes", """
ddtSchemes      { default Euler; }
gradSchemes     { default Gauss linear; }
divSchemes      { default none; }
laplacianSchemes { default Gauss linear corrected; }
interpolationSchemes { default linear; }
snGradSchemes   { default corrected; }""", location="system/solid")

    w(case, "system/solid/fvSolution", "dictionary", "fvSolution", """
solvers
{
    e
    {
        solver           PCG;
        preconditioner   DIC;
        tolerance        1e-08;
        relTol           0.1;
    }

    eFinal
    {
        $e;
        relTol           0;
    }
}

PIMPLE
{
    nNonOrthogonalCorrectors 0;
}""", location="system/solid")

    w(case, "constant/solid/physicalProperties", "dictionary",
      "physicalProperties", """
// Stainless steel
thermoType          constSolidThermo;

rho     { type uniform; value %g; }
Cv      { type uniform; value %g; }
kappa   { type uniform; value %g; }""" % (args.rhoS, args.cvS, args.kappaS),
      location="constant/solid")

    print("CHT test case: %s" % case)
    print("  fluid %.4g mm of burnt gas at %.4g K, solid %.4g mm of steel "
          "at %.4g K" % (args.hFluid, args.Tgas, args.hSolid, args.Tsolid))
    print("  mesh  %d x %d (fluid) + %d x %d (solid) cells of %.4g mm"
          % (nx, nyf, nx, nys, args.dx))


if __name__ == "__main__":
    main()
