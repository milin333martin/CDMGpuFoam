#!/usr/bin/env python3
"""
writeArrestorCase.py: generate a CDMGpuFoam case for a flame arrestor,
a stack of N narrow channels between parallel plates, in a wide channel.

Driven by runArrestor.sh; every setting comes in as a command-line
option so the script itself needs no editing.

Geometry (2D, lengths in mm on the command line, metres in the case):

    |<--- Lup --->|<- Lb ->|<--- Ldown --->|
    +-------------+========+---------------+   y = D     (wall)
    |             |  plate |               |
    |             +--------+               |   channel N
    |             |  plate |               |
    |    wide     +--------+    wide       |   ...
    |             |  plate |               |
    |             +--------+               |   channel 1
    |             |  plate |               |
    +-------------+========+---------------+   y = 0     (wall)

N channels of width d, N+1 plates of thickness (D - N d)/(N+1), so the
blockage ratio is BR = 1 - N d / D.
"""

import argparse
import math
import os
import shutil

HEAD = """/*--------------------------------*- C++ -*----------------------------------*\\
  =========                 |
  \\\\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\\\    /   O peration     | Website:  https://openfoam.org
    \\\\  /    A nd           | Version:  13
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
    with open(p, "w") as f:
        f.write(HEAD % (cls, loc, obj))
        f.write(body.strip() + "\n\n// ******************************"
                "******************************************* //\n")


def thermalBC_outer(args):
    spec = args.outerWallT
    if spec.lower() in ("adiabatic", "zerogradient"):
        return "zeroGradient"
    return "fixedValue; value uniform " + spec


def velBC(spec):
    return "noSlip" if spec.lower() == "noslip" else "slip"


def mappedPatch(case, region):
    """the region's interface patch, after splitMeshRegions"""
    import re
    b = os.path.join(case, "constant", region, "polyMesh", "boundary")
    for name, body in re.findall(r"^\s{4}(\w+)\s*$\s*\{([^}]*)\}",
                                 open(b).read(), re.M):
        if "mapped" in body:
            return name
    raise SystemExit("no mapped patch in " + b)


def main():
    a = argparse.ArgumentParser()
    g = a.add_argument
    g("--case", required=True)
    # geometry, mm
    g("--D", type=float, default=45.0, help="wide channel height")
    g("--d", type=float, default=0.4, help="channel (gap) width")
    g("--BR", type=float, default=None, help="blockage ratio, sets N")
    g("--N", type=int, default=None, help="number of channels (overrides BR)")
    g("--Lb", type=float, default=10.0, help="arrestor length")
    g("--Lup", type=float, default=120.0, help="length upstream of it")
    g("--Ldown", type=float, default=100.0, help="length downstream of it")
    g("--thickness", type=float, default=0.4, help="2D depth")
    # mesh
    g("--dx", type=float, default=0.4, help="base cell size")
    g("--levelsWide", type=int, default=3)
    g("--levelsPassage", type=int, default=5)
    g("--maxCells", type=int, default=6000000)
    # source snapshot and placement
    g("--sourceCase", default="../run_T6")
    g("--sourceTime", default="1e-04")
    g("--DCJ", type=float, default=1970.06)
    g("--sourceFront", type=float, default=74.5, help="front x in the source, mm")
    g("--frontX", type=float, default=40.0, help="where to put it, mm")
    g("--burntLength", type=float, default=8.0,
      help="mapped products kept behind the front, mm")
    g("--preRefineAhead", type=float, default=6.0, help="mm")
    g("--mirrorY", default="no")
    # run
    g("--endTime", type=float, default=1.4e-4)
    g("--writeInterval", type=float, default=5e-6)
    g("--maxCo", type=float, default=0.4)
    g("--fluxScheme", default="HLLC")
    g("--nProcs", type=int, default=24)
    g("--decompose", default="12 2")
    g("--distributor", default="zoltan")
    g("--fileHandler", default="uncollated")
    g("--foilDelta", type=float, default=0.05, help="soot foil pixel, mm")
    # CDM model
    g("--A", default="1.01e12")
    g("--EaByRT0", default="39.98")
    g("--qMwByRT0", default="62.97")
    g("--kappa0", default="2.00e-5")
    g("--T0", default="298")
    g("--gamma", default="1.125")
    g("--Mw", default="22.77")
    g("--Cv", default="2921.1990")
    g("--p0", default="101325")
    # walls
    g("--arrestorWallT", default="adiabatic",
      help="adiabatic, or a temperature in K for isothermal channel walls")
    g("--arrestorWallU", default="slip", help="slip | noSlip")
    g("--outerWallT", default="adiabatic",
      help="adiabatic, or a temperature in K for the outer channel walls")
    g("--outerWallU", default="slip", help="slip | noSlip")
    # conjugate heat transfer
    g("--cht", default="no", help="yes: plates as a solid region (CHT)")
    g("--stage", default="all", choices=("all", "fields"),
      help="fields: write 0/ after splitMeshRegions")
    g("--solidRho", type=float, default=7900.0)
    g("--solidCv", type=float, default=500.0)
    g("--solidKappa", type=float, default=16.0)
    g("--solidT", type=float, default=298.0)
    args = a.parse_args()

    CHT = args.cht.lower() in ("yes", "true", "1")
    case = args.case

    # ---------------- CHT: write 0/ after splitMeshRegions ----------------
    if args.stage == "fields":
        fp = mappedPatch(case, "fluid")
        sp = mappedPatch(case, "plates")

        def rfield(region, name, cls, dims, internal, patches, loc):
            w(case, "0/%s/%s" % (region, name), cls, name, """
dimensions      %s;

internalField   uniform %s;

boundaryField
{
%s
    frontAndBack    { type empty; }

    #includeEtc "caseDicts/setConstraintTypes"
}""" % (dims, internal, patches), location=loc)

        outer = """    internalFaces   { type internal; }
    outlet          { type zeroGradient; }
    inlet           { type zeroGradient; }
    walls           { type %s; }
"""
        zgOuter = outer % "zeroGradient"
        rfield("fluid", "p", "volScalarField", "[1 -1 -2 0 0 0 0]", args.p0,
               zgOuter + "    %s { type zeroGradient; }\n" % fp, "0/fluid")
        rfield("fluid", "T", "volScalarField", "[0 0 0 1 0 0 0]", args.T0,
               (outer % thermalBC_outer(args)) + """    %s
    {
        type            CDMCoupledTemperature;
        value           uniform %s;
    }
""" % (fp, args.T0), "0/fluid")
        rfield("fluid", "Y", "volScalarField", "[0 0 0 0 0 0 0]", "1",
               zgOuter + "    %s { type zeroGradient; }\n" % fp, "0/fluid")
        rfield("fluid", "U", "volVectorField", "[0 1 -1 0 0 0 0]", "(0 0 0)",
               """    internalFaces   { type internal; }
    outlet          { type zeroGradient; }
    inlet           { type zeroGradient; }
    walls           { type %s; }
    %s              { type %s; }
""" % (velBC(args.outerWallU), fp, velBC(args.arrestorWallU)), "0/fluid")
        rfield("plates", "T", "volScalarField", "[0 0 0 1 0 0 0]",
               args.solidT, """    walls           { type zeroGradient; }
    %s
    {
        type            coupledTemperature;
        value           uniform %g;
    }
""" % (sp, args.solidT), "0/plates")

        print("CHT fields written: fluid interface %s, solid interface %s"
              % (fp, sp))
        return

    mm = 1e-3
    D, d = args.D*mm, args.d*mm
    Lb, Lup, Ldown = args.Lb*mm, args.Lup*mm, args.Ldown*mm
    th = args.thickness*mm
    L = Lup + Lb + Ldown

    # number of channels
    if args.N is not None:
        N = args.N
    elif args.BR is not None:
        N = int(round((1 - args.BR)*D/d))
    else:
        raise SystemExit("give --N or --BR")
    if N < 1:
        raise SystemExit("N < 1: the arrestor would be solid")
    if N*d >= D:
        raise SystemExit("N*d >= D: the arrestor has no plates")

    t = (D - N*d)/(N + 1)          # plate thickness
    BR = 1 - N*d/D

    # y levels: plate, gap, plate, gap, ..., plate
    ys, kinds = [0.0], []
    for i in range(N):
        ys.append(ys[-1] + t); kinds.append("plate")
        ys.append(ys[-1] + d); kinds.append("gap")
    ys.append(D); kinds.append("plate")
    ys[-1] = D                      # snap
    nBands = len(kinds)

    xs = [0.0, Lup, Lup + Lb, L]

    def nc(length):
        return max(1, int(round(length/(args.dx*mm))))

    nxs = [nc(xs[1] - xs[0]), nc(xs[2] - xs[1]), nc(xs[3] - xs[2])]
    nys = [nc(ys[i + 1] - ys[i]) for i in range(nBands)]

    case = args.case
    if os.path.isdir(case):
        shutil.rmtree(case)
    os.makedirs(case)

    # ---------------- blockMeshDict ----------------
    nxv, nyv = 4, nBands + 1

    def V(ix, iy, iz):
        return iz*nxv*nyv + iy*nxv + ix

    verts = []
    for z in (0.0, th):
        for y in ys:
            for x in xs:
                verts.append("    (%.9g %.9g %.9g)" % (x, y, z))

    def hexb(ix, iy, zone, nx, ny):
        v = [V(ix, iy, 0), V(ix+1, iy, 0), V(ix+1, iy+1, 0), V(ix, iy+1, 0),
             V(ix, iy, 1), V(ix+1, iy, 1), V(ix+1, iy+1, 1), V(ix, iy+1, 1)]
        return ("    hex (%s) %s (%d %d 1) simpleGrading (1 1 1)"
                % (" ".join(str(i) for i in v), zone, nx, ny))

    blocks, present = [], []           # present: (ix, iy) blocks that exist
    for ix in (0, 1, 2):
        for iy in range(nBands):
            plate = (ix == 1 and kinds[iy] == "plate")
            if plate and not CHT:
                continue               # solid, not meshed
            if CHT:
                # one zone per region; the refinement zones (wide, passage)
                # are made inside the fluid region after splitMeshRegions
                zone = "plates" if plate else "fluid"
            else:
                zone = "passage" if (ix == 1) else "wide"
            blocks.append(hexb(ix, iy, zone, nxs[ix], nys[iy]))
            present.append((ix, iy))

    def face(ix, iy, side):
        q = {
            "x0": [V(ix, iy, 0), V(ix, iy, 1), V(ix, iy+1, 1), V(ix, iy+1, 0)],
            "x1": [V(ix+1, iy, 0), V(ix+1, iy+1, 0), V(ix+1, iy+1, 1), V(ix+1, iy, 1)],
            "y0": [V(ix, iy, 0), V(ix+1, iy, 0), V(ix+1, iy, 1), V(ix, iy, 1)],
            "y1": [V(ix, iy+1, 0), V(ix, iy+1, 1), V(ix+1, iy+1, 1), V(ix+1, iy+1, 0)],
            "z0": [V(ix, iy, 0), V(ix, iy+1, 0), V(ix+1, iy+1, 0), V(ix+1, iy, 0)],
            "z1": [V(ix, iy, 1), V(ix+1, iy, 1), V(ix+1, iy+1, 1), V(ix, iy+1, 1)],
        }[side]
        return "        (%s)" % " ".join(str(i) for i in q)

    outlet = [face(0, iy, "x0") for iy in range(nBands)]
    inlet = [face(2, iy, "x1") for iy in range(nBands)]
    walls = [face(0, 0, "y0"), face(2, 0, "y0"),
             face(0, nBands-1, "y1"), face(2, nBands-1, "y1")]
    arrestor = []
    if CHT:
        # the arrestor section is meshed, so its outer faces are walls too
        walls += [face(1, 0, "y0"), face(1, nBands-1, "y1")]
    else:
        for iy in range(nBands):
            if kinds[iy] == "plate":   # upstream/downstream faces of a plate
                arrestor.append(face(0, iy, "x1"))
                arrestor.append(face(2, iy, "x0"))
            else:                      # walls of a channel
                arrestor.append(face(1, iy, "y0"))
                arrestor.append(face(1, iy, "y1"))
    empty = []
    for (ix, iy) in present:
        empty += [face(ix, iy, "z0"), face(ix, iy, "z1")]

    w(case, "system/blockMeshDict", "dictionary", "blockMeshDict", """
// Flame arrestor: %d channels of %.3g mm in a %.4g mm channel,
// plates %.4g mm thick, blockage ratio %.4f, arrestor %.4g mm long.
// Base cell %.4g mm; AMR: %d levels in the wide zone, %d in the passages.
convertToMeters 1;

vertices
(
%s
);

blocks
(
%s
);

boundary
(
    internalFaces { type internal; faces (); }

    outlet
    {
        type patch;
        faces
        (
%s
        );
    }
    inlet
    {
        type patch;
        faces
        (
%s
        );
    }
    walls
    {
        type wall;
        faces
        (
%s
        );
    }
%s
    frontAndBack
    {
        type empty;
        faces
        (
%s
        );
    }
);""" % (N, args.d, args.D, t/mm, BR, args.Lb, args.dx,
         args.levelsWide, args.levelsPassage,
         "\n".join(verts), "\n".join(blocks), "\n".join(outlet),
         "\n".join(inlet), "\n".join(walls),
         ("" if CHT else """    arrestor
    {
        type wall;
        faces
        (
%s
        );
    }""" % "\n".join(arrestor)),
         "\n".join(empty)), location="system")

    # ---------------- constant ----------------
    rdir = "constant/fluid" if CHT else "constant"
    sdir = "system/fluid" if CHT else "system"

    if CHT:
        # solid region: stainless steel
        w(case, "constant/plates/physicalProperties", "dictionary",
          "physicalProperties", """
// Stainless steel plates
thermoType          constSolidThermo;

rho     { type uniform; value %g; }
Cv      { type uniform; value %g; }
kappa   { type uniform; value %g; }""" % (args.solidRho, args.solidCv,
                                          args.solidKappa),
          location="constant/plates")

        w(case, "system/plates/fvSchemes", "dictionary", "fvSchemes", """
ddtSchemes      { default Euler; }
gradSchemes     { default Gauss linear; }
divSchemes      { default none; }
laplacianSchemes { default Gauss linear corrected; }
interpolationSchemes { default linear; }
snGradSchemes   { default corrected; }""", location="system/plates")

        w(case, "system/plates/fvSolution", "dictionary", "fvSolution", """
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
}""", location="system/plates")

        # the refinement zones, made inside the fluid region after splitting
        w(case, "system/fluid/topoSetDict", "dictionary", "topoSetDict", """
actions
(
    {
        name    passage;
        type    cellZoneSet;
        action  new;
        source  boxToCell;
        box     (%.9g -1 -1) (%.9g 1 1);
    }
    {
        name    wide;
        type    cellZoneSet;
        action  new;
        source  boxToCell;
        box     (-1 -1 -1) (%.9g 1 1);
    }
    {
        name    wide;
        type    cellZoneSet;
        action  add;
        source  boxToCell;
        box     (%.9g -1 -1) (1 1 1);
    }
);""" % (xs[1], xs[2], xs[1], xs[2]), location="system/fluid")

        # case-level control for foamMultiRun
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

    w(case, rdir + "/physicalProperties", "dictionary", "physicalProperties", """
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
}""" % (args.Mw, args.Cv), location=rdir)

    w(case, rdir + "/CDMProperties", "dictionary", "CDMProperties", """
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
                       args.T0), location=rdir)

    w(case, rdir + "/dynamicMeshDict", "dictionary", "dynamicMeshDict", """
// Levels count from the base cell (%.4g mm):
//   wide     %d -> %.4g um
//   passage  %d -> %.4g um (%d cells across a %.3g mm channel)
topoChanger
{
    type            refiner2D;
    libs            ("libCDMDynamicMesh2D.so");

    refineInterval  2;

    nBufferLayers   2;
    maxCells        %d;

    correctFluxes   ((phi none));
    dumpLevel       true;

    refinementRegions
    {
        wide
        {
            cellZone         wide;
            field            refinementIndicator;
            lowerRefineLevel 0.05;
            upperRefineLevel 1e6;
            unrefineLevel    0.015;
            maxRefinement    %d;
        }
        passage
        {
            cellZone         passage;
            field            refinementIndicator;
            lowerRefineLevel 0.05;
            upperRefineLevel 1e6;
            unrefineLevel    0.015;
            maxRefinement    %d;
        }
    }
}

distributor
{
    type            loadBalancer;
    libs            ("libfvMeshDistributors.so");

    redistributionInterval  20;
    maxImbalance    0.1;
}""" % (args.dx, args.levelsWide, args.dx*1000/2**args.levelsWide,
        args.levelsPassage, args.d*1000/max(1, nys[1])/2**args.levelsPassage,
        nys[1]*2**args.levelsPassage, args.d,
        args.maxCells, args.levelsWide, args.levelsPassage), location=rdir)

    # ---------------- system ----------------
    fh = ("\nOptimisationSwitches\n{\n    fileHandler     collated;\n"
          "    maxThreadFileBufferSize 2e9;\n}\n"
          if args.fileHandler == "collated" else "")

    w(case, "system/controlDict", "dictionary", "controlDict", """
%s

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
%s""" % (("""application     foamMultiRun;

// the solid is constructed first: the fluid's coupled temperature condition
// needs the solid region to exist when it is built
regionSolvers
{
    plates          solid;
    fluid           CDMGpuFoam;
}""" if CHT else """application     foamRun;
solver          CDMGpuFoam;"""),
        args.endTime, args.writeInterval, args.maxCo, fh), location="system")

    w(case, sdir + "/fvSchemes", "dictionary", "fvSchemes", """
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
snGradSchemes   { default corrected; }""" % args.fluxScheme, location="system")

    bands, names = [], []
    for iy in range(nBands):
        if kinds[iy] == "gap":
            bands.append("            (%.9g %.9g)" % (ys[iy], ys[iy + 1]))
            names.append("ch%d" % (len(names) + 1))

    w(case, sdir + "/fvSolution", "dictionary", "fvSolution", """
solvers {}

PIMPLE
{
    nRKStages           2;
    viscous             yes;
    chemistry           yes;
    chemistrySplitting  Strang;
%s

    // Leading shock and leading flame, over the whole domain and inside
    // each channel of the arrestor
    frontTracking
    {
        direction       (1 0 0);
        pThreshold      5e5;
        YFlame          0.5;
        window          20;
        interval        20;
        bands
        (
%s
        );
        bandNames       (%s);
    }

    sootFoil
    {
        frameVelocity   (0 0 0);
        min             (0 0);
        max             (%.9g %.9g);
        delta           %.9g;
    }

    refinementIndicator
    {
        TReacting       500;
        YReacting       1e-3;
    }
}""" % ("    CHT                 yes;\n" if CHT else "",
        "\n".join(bands), " ".join(names), L, D, args.foilDelta*mm),
      location=sdir)

    nx, ny = args.decompose.split()
    decomposeParDict = """
numberOfSubdomains %d;

// decomposePar runs in serial, so the initial split is coordinate based
// (zoltan and ptscotch are parallel only)
decomposer      hierarchical;
hierarchicalCoeffs
{
    n           (%s %s 1);
    order       xyz;
}

// run-time rebalancing: zoltan is a parallel graph partitioner and finds
// splits that satisfy the refinement-history constraint, so it avoids the
// "unable to find desired decomposition split" warnings of hierarchical
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
}""" % (args.nProcs, nx, ny, args.distributor,
        ' "libzoltanDecomp.so"' if args.distributor == "zoltan" else
        (' "libptscotchDecomp.so"' if args.distributor == "ptscotch" else ""))

    w(case, sdir + "/decomposeParDict", "dictionary", "decomposeParDict",
      decomposeParDict, location=sdir)

    if CHT:
        # decomposePar -allRegions also reads the case-level dictionary
        w(case, "system/decomposeParDict", "dictionary", "decomposeParDict",
          decomposeParDict, location="system")

    shift = (args.frontX - args.sourceFront)*mm
    w(case, sdir + "/mapDetonationDict", "dictionary", "mapDetonationDict", """
sourceCase      "%s";
sourceTime      %s;
frameVelocity   (%g 0 0);          // shock frame -> lab frame
shift           (%.9g 0 0);        // source front -> x = %.9g m
mirrorY         %s;
mapMethod       nearest;
sourceRegion    region0;

// products further behind the front are replaced by a uniform burnt state
behindX         %.9g;

fresh
{
    p           %s;
    T           %s;
    U           (0 0 0);
    Y           1;
}

preRefine
{
    min         (%.9g 0);
    max         (%.9g %.9g);
    levels      %d;
}""" % (args.sourceCase, args.sourceTime, args.DCJ, shift, args.frontX*mm,
        args.mirrorY, (args.frontX - args.burntLength)*mm, args.p0, args.T0,
        (args.frontX - args.burntLength)*mm,
        (args.frontX + args.preRefineAhead)*mm, D, args.levelsWide),
      location="system")

    if CHT:
        w(case, "system/plates/decomposeParDict", "dictionary",
          "decomposeParDict", """
numberOfSubdomains %d;

decomposer      hierarchical;
hierarchicalCoeffs
{
    n           (%s %s 1);
    order       xyz;
}

distributor     hierarchical;""" % (args.nProcs, nx, ny),
          location="system/plates")

    # ---------------- 0 ----------------
    def thermalBC(spec):
        if spec.lower() in ("adiabatic", "zerogradient"):
            return "type            zeroGradient;"
        return ("type            fixedValue;\n        value           uniform %s;"
                % spec)

    def velocityBC(spec):
        if spec.lower() == "noslip":
            return "type            noSlip;"
        return "type            slip;"

    def field(name, cls, dims, internal, walls, arrestor):
        w(case, "0/" + name, cls, name, """
dimensions      %s;

internalField   uniform %s;

boundaryField
{
    internalFaces   { type internal; }
    outlet          { type zeroGradient; }
    inlet           { type zeroGradient; }
    walls
    {
        %s
    }
    arrestor
    {
        %s
    }
    frontAndBack    { type empty; }
}""" % (dims, internal, walls, arrestor), location="0")

    if CHT:
        nBase = sum(nxs[ix]*nys[iy] for (ix, iy) in present)
        print("Arrestor case (CHT): %s" % case)
        print("  channels   : N = %d of %.4g mm, plates %.4g mm, BR = %.4f"
              % (N, args.d, t/mm, BR))
        print("  regions    : fluid + plates (stainless steel: rho %g, "
              "Cv %g, kappa %g)" % (args.solidRho, args.solidCv,
                                    args.solidKappa))
        print("  base mesh  : %d cells (%.4g mm) over both regions" %
              (nBase, args.dx))
        print("  next       : splitMeshRegions, topoSet, then --stage fields")
        return

    zg = "type            zeroGradient;"
    field("p", "volScalarField", "[1 -1 -2 0 0 0 0]", args.p0, zg, zg)
    field("T", "volScalarField", "[0 0 0 1 0 0 0]", args.T0,
          thermalBC(args.outerWallT), thermalBC(args.arrestorWallT))
    field("Y", "volScalarField", "[0 0 0 0 0 0 0]", "1", zg, zg)
    field("U", "volVectorField", "[0 1 -1 0 0 0 0]", "(0 0 0)",
          velocityBC(args.outerWallU), velocityBC(args.arrestorWallU))

    nBase = sum(nxs[ix]*nys[iy] for (ix, iy) in present)
    print("Arrestor case: %s" % case)
    print("  channels   : N = %d of %.4g mm, plates %.4g mm, BR = %.4f"
          % (N, args.d, t/mm, BR))
    print("  domain     : %.4g x %.4g mm, arrestor at x = %.4g to %.4g mm"
          % (L/mm, D/mm, Lup/mm, (Lup + Lb)/mm))
    print("  base mesh  : %d cells (%.4g mm), finest %.4g um (wide), "
          "%.4g um (passage)"
          % (nBase, args.dx, args.dx*1000/2**args.levelsWide,
             (args.d/nys[1])*1000/2**args.levelsPassage))
    print("  channel    : %d base cells across, %d at the finest level"
          % (nys[1], nys[1]*2**args.levelsPassage))
    print("  walls      : channels %s / %s, outer %s / %s"
          % (args.arrestorWallT, args.arrestorWallU,
             args.outerWallT, args.outerWallU))
    print("  detonation : front at x = %.4g mm, run-in %.4g mm to the arrestor"
          % (args.frontX, Lup/mm - args.frontX))
    if Lup/mm - args.frontX < 50:
        print("  WARNING: less than 50 mm of run-in; the transplanted front "
              "may not have settled before it reaches the arrestor")


if __name__ == "__main__":
    main()
