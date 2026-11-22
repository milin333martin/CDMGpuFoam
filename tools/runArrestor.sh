#!/bin/sh
#------------------------------------------------------------------------------
# runArrestor.sh <config> [options]
#
# Builds and runs one arrestor case from a configuration file (see
# arrestor.cfg).  Every step is done here: blockMesh, the fluid/solid split
# and interface patches (CHT), the refinement zones, the initial fields, the
# transplanted quasi-steady detonation, the decomposition and the run.
#
#     ./runArrestor.sh myCase.cfg           build and run
#     ./runArrestor.sh myCase.cfg -g        build only, do not run
#     ./runArrestor.sh myCase.cfg -f        overwrite an existing case
#     ./runArrestor.sh myCase.cfg -r        restart: just run what is there
#     ./runArrestor.sh myCase.cfg -n other  another case name
#
# A copy of the configuration is written into the case as case.cfg.
#------------------------------------------------------------------------------
set -e

SCRIPT_DIR=$(cd "$(dirname "$0")" && pwd)

usage() { sed -n '3,16p' "$0"; exit ${1:-0}; }

[ $# -ge 1 ] || usage 1
CONFIG=$1; shift
[ -f "$CONFIG" ] || { echo "no such configuration: $CONFIG" >&2; exit 1; }
CONFIG=$(cd "$(dirname "$CONFIG")" && pwd)/$(basename "$CONFIG")

BUILD_ONLY=no; FORCE=no; RESTART=no; NAME=
while getopts "gfrn:h" opt; do
    case $opt in
        g) BUILD_ONLY=yes ;;
        f) FORCE=yes ;;
        r) RESTART=yes ;;
        n) NAME=$OPTARG ;;
        h) usage ;;
        *) usage 1 ;;
    esac
done

# ---------------------------------------------------------------- settings --
. "$CONFIG"
[ -n "$NAME" ] && CASE=$NAME

# ---------------------------------------------------------------- checks ----
fail() { echo "ERROR: $*" >&2; exit 1; }

command -v foamRun > /dev/null || fail "OpenFOAM is not in the environment"
command -v mapDetonation > /dev/null \
    || fail "mapDetonation not found: build CDMGpuFoam first (./Allwmake)"

if [ "$RESTART" = no ]; then
    case "$SOURCE_CASE" in
        /*) ;;
        *) SOURCE_CASE=$(cd "$(dirname "$CONFIG")" && cd "$SOURCE_CASE" \
                         2>/dev/null && pwd) || true ;;
    esac
    [ -d "$SOURCE_CASE" ] || fail "SOURCE_CASE does not exist: $SOURCE_CASE"
    [ -n "$(ls -d "$SOURCE_CASE"/[0-9]* 2>/dev/null)" ] \
        || fail "SOURCE_CASE has no time directories (reconstructPar first)"
fi

NX=$(echo "$DECOMPOSE" | cut -d' ' -f1)
NY=$(echo "$DECOMPOSE" | cut -d' ' -f2)
[ $((NX*NY)) -eq "$NPROCS" ] \
    || fail "DECOMPOSE ($DECOMPOSE) does not multiply to NPROCS ($NPROCS)"

[ -n "$N" ] && NOPT="--N $N" || NOPT="--BR $BR"

cd "$(dirname "$CONFIG")"

if [ "$RESTART" = no ] && [ -d "$CASE" ]; then
    [ "$FORCE" = yes ] || fail "$CASE exists (use -f to overwrite, -r to restart)"
    rm -rf "$CASE"
fi

# ---------------------------------------------------------------- build -----
gen()
{
    python3 "$SCRIPT_DIR/writeArrestorCase.py" --case "$CASE" "$@" \
        --D "$D" --d "$d" $NOPT --Lb "$LB" --Lup "$LUP" --Ldown "$LDOWN" \
        --thickness "$THICKNESS" \
        --dx "$DX" --levelsWide "$LEVELS_WIDE" \
        --levelsPassage "$LEVELS_PASSAGE" --maxCells "$MAXCELLS" \
        --sourceCase "$SOURCE_CASE" --sourceTime "$SOURCE_TIME" \
        --sourceFront "$SOURCE_FRONT" --DCJ "$DCJ" --frontX "$FRONT_X" \
        --burntLength "$BURNT_LENGTH" --mirrorY "$MIRROR_Y" \
        --endTime "$END_TIME" --writeInterval "$WRITE_INTERVAL" \
        --maxCo "$MAX_CO" --fluxScheme "$FLUX_SCHEME" --nProcs "$NPROCS" \
        --decompose "$DECOMPOSE" --distributor "$DISTRIBUTOR" \
        --fileHandler "$FILE_HANDLER" --foilDelta "$FOIL_DELTA" \
        --arrestorWallT "$ARRESTOR_WALL_T" --arrestorWallU "$ARRESTOR_WALL_U" \
        --outerWallT "$OUTER_WALL_T" --outerWallU "$OUTER_WALL_U" \
        --cht "$CHT" --solidRho "$SOLID_RHO" --solidCv "$SOLID_CV" \
        --solidKappa "$SOLID_KAPPA" --solidT "$SOLID_T" \
        --A "$CDM_A" --EaByRT0 "$CDM_EaByRT0" --qMwByRT0 "$CDM_qMwByRT0" \
        --kappa0 "$CDM_kappa0" --T0 "$CDM_T0" --gamma "$CDM_GAMMA" \
        --Mw "$CDM_MW" --Cv "$CDM_CV" --p0 "$P0"
}

run()
{
    echo "  $*"
    if ! "$@" > "log.$1" 2>&1; then
        echo "FAILED: $* (see $CASE/log.$1)" >&2
        tail -20 "log.$1" >&2
        exit 1
    fi
}

if [ "$RESTART" = no ]; then
    gen
    cp "$CONFIG" "$CASE/case.cfg"
    cd "$CASE"

    run blockMesh
    grep -E "^  nCells" log.blockMesh || true

    if [ "$CHT" = yes ]; then
        # one region per cell zone; -noFields because 0/ is written below
        run splitMeshRegions -cellZonesOnly -overwrite -noFields
        grep "added patches" -A2 log.splitMeshRegions || true

        # the split deletes empty patches, including the one the load
        # balancer needs
        python3 "$SCRIPT_DIR/addInternalPatchToMesh.py" . fluid
        rm -f 0/cellToRegion 0/*/cellToRegion

        # the refinement zones live inside the fluid region
        run topoSet -region fluid

        # 0/ needs the interface patch names, so it comes after the split
        cd ..
        gen --stage fields
        cd "$CASE"

        run mapDetonation -region fluid
    else
        run mapDetonation
    fi
    grep -E "Method|uniform burnt" log.mapDetonation || true

    if [ "$NPROCS" -gt 1 ]; then
        # the decomposition is always uncollated: the history utility and
        # addInternalPatchToMesh read the processor directories
        if [ "$CHT" = yes ]; then
            run decomposePar -allRegions -fileHandler uncollated
            run decomposeRefinementHistory2D -region fluid \
                -fileHandler uncollated
        else
            run decomposePar -fileHandler uncollated
            run decomposeRefinementHistory2D -fileHandler uncollated
        fi
        grep "Decomposed" log.decomposeRefinementHistory2D || true
    fi
else
    cd "$CASE"
fi

if [ "$BUILD_ONLY" = yes ]; then
    echo "built $CASE (not run: -g)"
    exit 0
fi

# ---------------------------------------------------------------- run -------
[ "$CHT" = yes ] && APP=foamMultiRun || APP=foamRun

if [ "$NPROCS" -gt 1 ]; then
    echo "  $APP on $NPROCS processors"
    nohup mpirun -np "$NPROCS" --bind-to core "$APP" -parallel \
        -fileHandler "$FILE_HANDLER" > "log.$APP" 2>&1 &
else
    nohup "$APP" -fileHandler "$FILE_HANDLER" > "log.$APP" 2>&1 &
fi

cat << EOF

started $CASE

  progress   tail -f $CASE/log.$APP
  fronts     $CASE/postProcessing/fronts/0/fronts.dat
  channels   $SCRIPT_DIR/channels.sh $CASE
  soot foil  $CASE/postProcessing/sootFoil/<time>/sootFoil.vtk
  settings   $CASE/case.cfg
EOF
