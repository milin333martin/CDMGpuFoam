/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | CDMGpuFoam: Chemical-Diffusive Model solver
   \\    /   O peration     | for OpenFOAM-13 (openfoam.org)
    \\  /    A nd           |
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of CDMGpuFoam, derived from OpenFOAM (shockFluid).
    GNU General Public License v3 or later; see COPYING.

\*---------------------------------------------------------------------------*/

#include "cdmFoam.H"
#include "linear.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace solvers
{
    defineTypeNameAndDebug(cdmFoam, 0);
    addToRunTimeSelectionTable(solver, cdmFoam, fvMesh);
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::solvers::cdmFoam::readControls()
{
    const dictionary& dict = pimple.dict();

    nRKStages_ = dict.lookupOrDefault<label>("nRKStages", 2);

    if (nRKStages_ < 1 || nRKStages_ > 3)
    {
        FatalIOErrorInFunction(dict)
            << "nRKStages = " << nRKStages_ << " is not supported;"
            << " choose 1 (Euler), 2 (SSP-RK2) or 3 (SSP-RK3)"
            << exit(FatalIOError);
    }

    viscous_ = dict.lookupOrDefault<bool>("viscous", true);
    chemistry_ = dict.lookupOrDefault<bool>("chemistry", true);

    const word splitting
    (
        dict.lookupOrDefault<word>("chemistrySplitting", "Strang")
    );

    if (splitting == "Strang")
    {
        strang_ = true;
    }
    else if (splitting == "Godunov")
    {
        strang_ = false;
    }
    else
    {
        FatalIOErrorInFunction(dict)
            << "Unknown chemistrySplitting " << splitting
            << "; valid choices: Strang, Godunov" << exit(FatalIOError);
    }

    {
        const dictionary rDict
        (
            dict.subOrEmptyDict("refinementIndicator")
        );
        TReacting_ = rDict.lookupOrDefault<scalar>("TReacting", 500);
        YReacting_ = rDict.lookupOrDefault<scalar>("YReacting", 1e-3);
    }

    if (dict.isDict("sootFoil") && !sootFoil_.valid())
    {
        sootFoil_.reset(new CDM::sootFoil(mesh, dict.subDict("sootFoil")));
    }

    {
        const dictionary pDict(dict.subOrEmptyDict("positivityLimiter"));
        positivityLimiter_ = pDict.lookupOrDefault<bool>("active", true);
        if (pDict.found("rhoMin"))
        {
            rhoFloor_ = pDict.lookup<scalar>("rhoMin");
        }
        if (pDict.found("pMin"))
        {
            pFloor_ = pDict.lookup<scalar>("pMin");
        }
    }

    {
        const dictionary dDict(dict.subOrEmptyDict("diagnostics"));
        TWarn_ = dDict.lookupOrDefault<scalar>("TWarn", 6000);
        TFail_ = dDict.lookupOrDefault<scalar>("TFail", 30000);
    }

    CHT_ = dict.lookupOrDefault<bool>("CHT", false);

    frontTracking_ = dict.isDict("frontTracking");

    if (frontTracking_)
    {
        const dictionary& fDict = dict.subDict("frontTracking");

        frontDirection_ = fDict.lookupOrDefault<vector>("direction", vector(1, 0, 0));
        frontDirection_ /= mag(frontDirection_) + vSmall;
        frontPThreshold_ = fDict.lookupOrDefault<scalar>("pThreshold", 5e5);
        YFlame_ = fDict.lookupOrDefault<scalar>("YFlame", 0.5);
        frontWindow_ = max(fDict.lookupOrDefault<label>("window", 20), 2);
        frontInterval_ = max(fDict.lookupOrDefault<label>("interval", 20), 1);

        frontBands_.clear();
        frontNames_.clear();

        if (fDict.found("bands"))
        {
            const List<Pair<scalar>> bands(fDict.lookup("bands"));
            frontBands_ = bands;
        }

        if (fDict.found("bandNames"))
        {
            frontNames_ = wordList(fDict.lookup("bandNames"));
        }

        if (frontNames_.size() != frontBands_.size())
        {
            frontNames_.setSize(frontBands_.size());
            forAll(frontNames_, i)
            {
                frontNames_[i] = "band" + Foam::name(i);
            }
        }
    }

    shockTracking_ = dict.isDict("shockTracking");

    if (shockTracking_)
    {
        const dictionary& sDict = dict.subDict("shockTracking");
        shockDirection_ = sDict.lookup<vector>("direction");
        shockDirection_ /= mag(shockDirection_);
        shockPThreshold_ = sDict.lookup<scalar>("pThreshold");
        shockWindow_ = max(sDict.lookupOrDefault<label>("window", 100), 2);
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::solvers::cdmFoam::cdmFoam(fvMesh& mesh)
:
    fluidSolver(mesh),

    thermoPtr_(psiThermo::New(mesh)),

    thermo_(thermoPtr_()),

    p_(thermo_.p()),

    rho_
    (
        IOobject
        (
            "rho",
            runTime.name(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        thermo_.renameRho()
    ),

    U_
    (
        IOobject
        (
            "U",
            runTime.name(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    ),

    Y_
    (
        IOobject
        (
            "Y",
            runTime.name(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    ),

    rhoU_
    (
        IOobject
        (
            "rhoU",
            runTime.name(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        rho_*U_
    ),

    rhoE_
    (
        IOobject
        (
            "rhoE",
            runTime.name(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        rho_*(thermo_.he() + 0.5*magSqr(U_))
    ),

    rhoY_
    (
        IOobject
        (
            "rhoY",
            runTime.name(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        rho_*Y_
    ),

    phi_
    (
        IOobject
        (
            "phi",
            runTime.name(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        linearInterpolate(rho_*U_) & mesh.Sf()
    ),

    cdm_
    (
        mesh,
        gAverage(thermo_.Cp().primitiveField())
      - gAverage(thermo_.Cv().primitiveField())
    ),

    fluxScheme_(CDM::fluxScheme::New(mesh)),

    Qdot_
    (
        IOobject
        (
            "Qdot",
            runTime.name(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar(dimPower/dimVolume, 0)
    ),

    pMax_
    (
        IOobject
        (
            "pMax",
            runTime.name(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        p_
    ),

    refinementIndicator_
    (
        IOobject
        (
            "refinementIndicator",
            runTime.name(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar(dimless, 0)
    ),

    nRKStages_(2),
    viscous_(true),
    chemistry_(true),
    strang_(true),

    CHT_(false),
    KPtr_(nullptr),

    frontTracking_(false),
    frontDirection_(1, 0, 0),
    frontPThreshold_(5e5),
    YFlame_(0.5),
    frontWindow_(20),
    frontInterval_(20),
    frontBands_(),
    frontNames_(),
    frontStartName_(mesh.time().name()),
    frontFilePtr_(nullptr),

    shockTracking_(false),
    shockDirection_(1, 0, 0),
    shockPThreshold_(great),
    shockWindow_(100),

    TReacting_(500),
    YReacting_(1e-3),

    gamma_(1.4),
    eOffset_(0),

    positivityLimiter_(true),
    rhoFloor_(-1),
    pFloor_(-1),
    nLimitedFaces_(0),

    TWarn_(6000),
    TFail_(30000),

    tRHS_(0),
    tStage_(0),
    tChem_(0),
    tTail_(0),
    tDiag_(0),
    tSoot_(0),
    tCourant_(0),
    tDec_{0, 0, 0, 0, 0, 0, 0, 0},

    thermo(thermo_),
    p(p_),
    rho(rho_),
    U(U_),
    Y(Y_),
    phi(phi_)
{
    thermo.validate(type(), "e");

    if (!transient())
    {
        FatalErrorInFunction
            << type() << " is an explicit time-accurate solver;"
            << " steady-state and LTS ddtSchemes are not supported"
            << exit(FatalError);
    }

    readControls();

    // Equation of state constants (uniform for the CDM gas)
    gamma_ =
        gAverage(thermo_.Cp().primitiveField())
       /gAverage(thermo_.Cv().primitiveField());

    eOffset_ = gAverage
    (
        thermo_.he().primitiveField()
      - p_.primitiveField()/(rho_.primitiveField()*(gamma_ - 1))
    );

    Info<< "CDM gas: gamma = " << gamma_ << ", e = p/(rho (gamma - 1)) + "
        << eOffset_ << " J/kg" << nl << endl;

    // Make the primitive and conservative states consistent: clips Y, applies
    // fvConstraints, evaluates BCs. With CHT the boundary evaluation is
    // deferred to the first time step, when every region exists.
    decode(!CHT_);

    if (CHT_)
    {
        // The coupled temperature condition needs K as soon as it is first
        // evaluated, which is before the first call to advance()
        updateK(cdm_.kappaByCp(thermo_.T()));
    }

    if (rhoFloor_ < 0)
    {
        rhoFloor_ = 1e-3*gMin(rho_.primitiveField());
    }
    if (pFloor_ < 0)
    {
        pFloor_ = 1e-3*gMin(p_.primitiveField());
    }

    Info<< "Positivity limiter: " << (positivityLimiter_ ? "on" : "off")
        << ", rhoMin = " << rhoFloor_ << ", pMin = " << pFloor_ << nl
        << endl;

    updateDiagnostics();

    if (sootFoil_.valid())
    {
        sootFoil_->update(p_);
    }

    correctCoNum(waveSpeedSf());
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::solvers::cdmFoam::~cdmFoam()
{
    reportTimeSplit();
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::solvers::cdmFoam::preSolve()
{
    readControls();

    fvModels().preUpdateMesh();

    // Update the mesh for topology change and redistribution. All registered
    // volFields, including the conservative variables, are mapped here.
    mesh_.update();

    if (mesh.dynamic())
    {
        // Rebuild the primitive variables from the mapped conservative ones
        decode();
    }

    // Time-step limit evaluated on the (possibly refined) current mesh
    const double t0 = CDMWallTime();

    correctCoNum(waveSpeedSf());

    tCourant_ += CDMWallTime() - t0;
}


void Foam::solvers::cdmFoam::moveMesh()
{
    if (pimple.firstIter() || pimple.moveMeshOuterCorrectors())
    {
        mesh_.move();
    }
}


void Foam::solvers::cdmFoam::motionCorrector()
{
    if (pimple.firstIter() || pimple.moveMeshOuterCorrectors())
    {
        if (mesh.changing())
        {
            meshCourantNo();
        }
    }
}


void Foam::solvers::cdmFoam::prePredictor()
{
    // The explicit update is performed once per time step
    if (pimple.firstIter())
    {
        advance();
    }
}


void Foam::solvers::cdmFoam::momentumTransportPredictor()
{}


void Foam::solvers::cdmFoam::thermophysicalTransportPredictor()
{}


void Foam::solvers::cdmFoam::momentumPredictor()
{}


void Foam::solvers::cdmFoam::thermophysicalPredictor()
{}


void Foam::solvers::cdmFoam::pressureCorrector()
{}


void Foam::solvers::cdmFoam::momentumTransportCorrector()
{}


void Foam::solvers::cdmFoam::thermophysicalTransportCorrector()
{}


void Foam::solvers::cdmFoam::postSolve()
{
    const double t0 = CDMWallTime();

    updateDiagnostics();

    const double t1 = CDMWallTime();
    tDiag_ += t1 - t0;

    if (sootFoil_.valid())
    {
        sootFoil_->update(p_);

        if (runTime.writeTime())
        {
            sootFoil_->write();
        }
    }

    tSoot_ += CDMWallTime() - t1;
}


// ************************************************************************* //
