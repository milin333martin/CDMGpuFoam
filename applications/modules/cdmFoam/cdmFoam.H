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

Class
    Foam::solvers::cdmFoam

Description
    Density-based solver module for unsteady, multidimensional, reactive,
    compressible Navier-Stokes flow with the chemical-diffusive model (CDM)
    of Kaplan et al. (2019) / Lu, Kaplan & Oran (2021).

    Run with
    \verbatim
        foamRun -solver cdmFoam      (or: solver cdmFoam; in controlDict)
    \endverbatim

    Numerics
      - Conservative variables rho, rhoU, rhoE, rhoY (Y: reactant fraction)
      - MUSCL reconstruction + approximate Riemann solver (HLLC, HLLCP)
      - Explicit viscous stress, heat conduction and reactant diffusion
        with the CDM power-law transport (K = kappa0 Cp T^n)
      - Explicit SSP Runge-Kutta, 1-3 stages
      - Reaction Strang- (or Godunov-) split, integrated per cell with an
        exponential sub-stepped update (see CDMModel)

    Adaptive mesh refinement / load balancing
      Conservative variables are registered fields and are mapped by the
      mesh topology changer; primitives are rebuilt from them after every
      mesh update, and the time-step limit is evaluated on the updated mesh.
      A dimensionless refinement indicator field "refinementIndicator"
      (maximum relative density jump across the faces of a cell) is provided
      for the refiner.

    Controls: system/fvSolution
    \verbatim
    PIMPLE
    {
        nRKStages           2;          // 1: Euler, 2: SSP-RK2, 3: SSP-RK3
        viscous             yes;        // no: reactive Euler equations
        chemistry           yes;
        chemistrySplitting  Strang;     // Strang | Godunov

        CHT                 yes;        // optional: publish K for a
                                        // CDMCoupledTemperature patch

        frontTracking                   // optional; replaces shockTracking
        {
            direction       (1 0 0);
            pThreshold      5e5;        // Pa, leading shock
            YFlame          0.5;        // leading flame: tip of Y < YFlame
            window          20;         // samples averaged for the speeds
            interval        20;         // steps between samples
            bands                       // optional: tracked separately
            (
                (0.0223 0.0227)
            );
            bandNames       (ch0);      // optional
        }

        shockTracking                   // optional, superseded by frontTracking
        {
            direction       (1 0 0);
            pThreshold      5e5;        // Pa
            window          100;        // time steps for speed average
        }

        sootFoil                        // optional: see CDM::sootFoil
        {
            frameVelocity   (1967.15 0 0);  // shock frame; (0 0 0) lab frame
            min             (0 0);
            max             (0.135 0.01);
            delta           2e-5;
            interval        1;              // steps between samples
        }

        positivityLimiter               // optional, defaults shown
        {
            active          yes;
            rhoMin          <1e-3 min(rho(t=0))>;
            pMin            <1e-3 min(p(t=0))>;
        }

        diagnostics                     // optional, defaults shown
        {
            TWarn           6000;       // report hot-spot location above
            TFail           30000;      // stop with cell description above
        }

        refinementIndicator             // optional, defaults shown
        {
            TReacting       500;        // K
            YReacting       1e-3;
        }
    }
    \endverbatim

    refinementIndicator (for the AMR refiner, dimensionless):
        max over the faces of a cell of
            |rho_N - rho_P|/min(rho_N, rho_P),
            |p_N - p_P|/min(p_N, p_P),
            |Y_N - Y_P|
        and 1 in reacting gas (T > TReacting and Y > YReacting), which keeps
        the induction zone refined where the gradients are small.
        Typical refiner settings: lowerRefineLevel 0.05; upperRefineLevel 1e6;

    \verbatim
    \endverbatim

    Additional written fields: Y, Qdot (heat-release rate W/m^3),
    pMax (maximum pressure history: numerical smoke foil),
    refinementIndicator.

SourceFiles
    cdmFoam.C
    advance.C
    decode.C
    diagnostics.C

\*---------------------------------------------------------------------------*/

#ifndef cdmFoam_H
#define cdmFoam_H

#include "fluidSolver.H"
#include "psiThermo.H"
#include "CDMModel.H"
#include "fluxScheme.H"
#include "OFstream.H"
#include "Pair.H"
#include "sootFoil.H"

#include <chrono>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace solvers
{

//- Monotonic wall-clock seconds since the first call, for the solver's own
//  per-step time split.  Not a substitute for a profiler; it is here so that
//  the split between the flux evaluation, the chemistry and the diagnostics
//  is a measurement rather than an assumption.
inline double CDMWallTime()
{
    using cdmClock = std::chrono::steady_clock;
    static const cdmClock::time_point t0 = cdmClock::now();
    return std::chrono::duration<double>(cdmClock::now() - t0).count();
}


/*---------------------------------------------------------------------------*\
                          Class cdmFoam Declaration
\*---------------------------------------------------------------------------*/

class cdmFoam
:
    public fluidSolver
{

protected:

    // Thermophysical properties

        //- Pointer to the fluid thermophysical properties
        autoPtr<psiThermo> thermoPtr_;

        //- Reference to the fluid thermophysical properties
        psiThermo& thermo_;

        //- Reference to the pressure field
        volScalarField& p_;

        //- Density
        volScalarField rho_;


    // Primitive variables

        //- Velocity
        volVectorField U_;

        //- Reactant mass fraction
        volScalarField Y_;


    // Conservative variables (registered so that they are mapped by AMR)

        volVectorField rhoU_;

        volScalarField rhoE_;

        volScalarField rhoY_;

        //- Mass flux (for flux-dependent boundary conditions)
        surfaceScalarField phi_;


    // Models

        //- Chemical-diffusive model
        CDM::CDMModel cdm_;

        //- Riemann solver
        autoPtr<CDM::fluxScheme> fluxScheme_;

        //- Optional in-situ lab-frame soot foil
        autoPtr<CDM::sootFoil> sootFoil_;


    // Diagnostic fields

        //- Heat-release rate [W/m^3]
        volScalarField Qdot_;

        //- Maximum pressure history (numerical smoke foil)
        volScalarField pMax_;

        //- Refinement indicator (max relative density jump)
        volScalarField refinementIndicator_;


    // Controls

        label nRKStages_;

        bool viscous_;

        bool chemistry_;

        bool strang_;


    // Shock tracking

        bool shockTracking_;

        vector shockDirection_;

        scalar shockPThreshold_;

        label shockWindow_;

        DynamicList<scalar> shockTimes_;

        DynamicList<scalar> shockPositions_;


    // Refinement indicator controls

        scalar TReacting_;

        scalar YReacting_;


    // Equation of state constants of the calorically perfect CDM gas
    // e = p/(rho (gamma - 1)) + eOffset

        scalar gamma_;

        scalar eOffset_;


    // Conjugate heat transfer

        //- Publish the CDM thermal conductivity K for the coupled
        //  temperature boundary condition of an adjacent solid region
        bool CHT_;

        autoPtr<volScalarField> KPtr_;


    // Shock and flame front tracking

        bool frontTracking_;

        vector frontDirection_;

        scalar frontPThreshold_;

        scalar YFlame_;

        //- Samples used for the front speeds, and their spacing in steps
        label frontWindow_;

        label frontInterval_;

        //- Bands (y ranges) tracked separately, e.g. arrestor channels
        List<Pair<scalar>> frontBands_;

        wordList frontNames_;

        DynamicList<scalar> frontTimes_;

        DynamicList<scalarField> frontShock_;

        DynamicList<scalarField> frontFlame_;

        word frontStartName_;

        autoPtr<OFstream> frontFilePtr_;


    // Positivity-preserving flux limiter (Hu, Adams & Shu 2013)

        bool positivityLimiter_;

        //- Density and pressure floors (default 1e-3 of the initial minima)
        scalar rhoFloor_;

        scalar pFloor_;

        //- Number of limited faces in the current time step
        label nLimitedFaces_;


    // Diagnostics controls

        //- Report the location of the maximum temperature above TWarn
        scalar TWarn_;

        //- Stop with a description of the cell if T exceeds TFail or the
        //  density/temperature becomes non-positive
        scalar TFail_;


    // Wall-clock accumulators in seconds for the per-step time split, reported
    // once at the end of the run, at a cost of a dozen clock reads per step.
    // The balance between the flux evaluation, the chemistry and the
    // diagnostics is what decides whether moving work to a device is worth it.

        double tRHS_;

        double tStage_;

        double tChem_;

        double tTail_;

        double tDiag_;

        double tSoot_;

        double tCourant_;

        //- The eight segments of decode(), in the order they run.  decode()
        //  is the single largest piece of host-side work in a GPU run and
        //  which of its parts has to move to the device is a question for
        //  the clock, not for intuition.
        double tDec_[8];


protected:

    // Protected Member Functions. These were private until the GPU back end
    // needed to replace the right-hand side; see applications/modules/cdmGpuFoam.

        //- Read the solution controls from PIMPLE
        void readControls();

        //- Primitive -> conservative (including boundaries)
        void encode();

        //- Conservative -> primitive, thermo update, BCs, re-encode
        void decode(const bool evaluateBCs = true);

        //- Speed of sound
        tmp<volScalarField> soundSpeed() const;

        //- Maximum wave speed times face area of the current state
        tmp<surfaceScalarField> waveSpeedSf() const;

        //- Update the (convective and diffusive) Courant number
        void correctCoNum(const surfaceScalarField& amaxSf);

        //- Evaluate the right-hand sides of the conservation laws
        //  (per unit volume) for the current state and a forward-Euler
        //  stage of length deltaT (used by the positivity limiter)
        virtual void computeRHS
        (
            const scalar deltaT,
            scalarField& dRho,
            vectorField& dRhoU,
            scalarField& dRhoE,
            scalarField& dRhoY
        );

        //- Update the published conductivity field (CHT)
        void updateK(const volScalarField& kByCp);

        //- Track the leading shock and flame, globally and per band
        void trackFronts();

        //- Positivity-preserving limiting of the face fluxes
        void limitPositivity
        (
            const scalar deltaT,
            surfaceScalarField& rhoPhi,
            surfaceVectorField& rhoUPhi,
            surfaceScalarField& rhoEPhi,
            surfaceScalarField& rhoYPhi
        );

        //- Integrate the reaction over deltaT and decode
        virtual void react(const scalar deltaT);

        //- Advance the flow over one time step.
        //  Virtual so that a back end which keeps the state on a device can
        //  run the whole step there and hand the host a consistent state at
        //  the step boundary; see applications/modules/cdmGpuFoam.
        virtual void advance();

        //- Update pMax, the refinement indicator and shock tracking
        void updateDiagnostics();

        //- Check the conservative state before the thermo update; stop
        //  with a description of the first bad cells
        void checkState() const;

        //- Report the location of the hottest cell above TWarn
        void reportHotSpot() const;

        //- Write the per-step wall-clock split at the end of the run
        void reportTimeSplit() const;


public:

    // Public Data

        //- Reference to the fluid thermophysical properties
        const psiThermo& thermo;

        //- Reference to the pressure field
        const volScalarField& p;

        //- Reference to the density field
        const volScalarField& rho;

        //- Reference to the velocity field
        const volVectorField& U;

        //- Reference to the reactant mass fraction
        const volScalarField& Y;

        //- Reference to the mass-flux field
        const surfaceScalarField& phi;


    //- Runtime type information
    TypeName("cdmFoam");


    // Constructors

        //- Construct from region mesh
        cdmFoam(fvMesh& mesh);

        //- Disallow default bitwise copy construction
        cdmFoam(const cdmFoam&) = delete;


    //- Destructor
    virtual ~cdmFoam();


    // Member Functions

        //- Called at the start of the time-step, before the PIMPLE loop
        virtual void preSolve();

        //- Called at the start of the PIMPLE loop to move the mesh
        virtual void moveMesh();

        //- Corrections that follow mesh motion
        virtual void motionCorrector();

        //- Called at the start of the PIMPLE loop: advances the solution
        virtual void prePredictor();

        //- Not used
        virtual void momentumTransportPredictor();

        //- Not used
        virtual void thermophysicalTransportPredictor();

        //- Not used (explicit update in prePredictor)
        virtual void momentumPredictor();

        //- Not used (explicit update in prePredictor)
        virtual void thermophysicalPredictor();

        //- Not used (explicit update in prePredictor)
        virtual void pressureCorrector();

        //- Not used
        virtual void momentumTransportCorrector();

        //- Not used
        virtual void thermophysicalTransportCorrector();

        //- Called after the PIMPLE loop at the end of the time-step
        virtual void postSolve();


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const cdmFoam&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace solvers
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
