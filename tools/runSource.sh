#!/bin/sh
#------------------------------------------------------------------------------
# runSource.sh <config> [options]
#
# Runs the wave-fixed (shock-frame) detonation that later initialises an
# arrestor or narrow-channel case, and hands its result on.
#
#     ./runSource.sh src.cfg              build and run
#     ./runSource.sh src.cfg -g           build only
#     ./runSource.sh src.cfg -f           overwrite an existing case
#     ./runSource.sh src.cfg -r           restart what is already built
#     ./runSource.sh src.cfg -s           progress: front, speed, cells, cells/foil
#     ./runSource.sh src.cfg -F           finish: reconstruct the latest time,
#                                         measure the front, report the snapshot
#     ./runSource.sh src.cfg -F -t 8e-05  finish at a particular time
#     ./runSource.sh src.cfg -F -a a.cfg  the same, and write the measured
#                                         snapshot into the arrestor config
#                                         a.cfg (SOURCE_CASE, SOURCE_TIME,
#                                         SOURCE_FRONT, DCJ and the CDM block)
#
# The CJ speed and the half-reaction length come from CDM_ZND.py, so the
# configuration carries the model parameters, not derived numbers.
#------------------------------------------------------------------------------
set -e
SCRIPT_DIR=$(cd "$(dirname "$0")" && pwd)

usage() { sed -n '3,25p' "$0"; exit ${1:-0}; }
fail()  { echo "ERROR: $*" >&2; exit 1; }

[ $# -ge 1 ] || usage 1
CONFIG=$1; shift
[ -f "$CONFIG" ] || fail "no such configuration: $CONFIG"
CONFIG=$(cd "$(dirname "$CONFIG")" && pwd)/$(basename "$CONFIG")

BUILD_ONLY=no; FORCE=no; RESTART=no; STATUS=no; FINISH=no; FTIME=; ARR=; NAME=
while getopts "gfrsFt:a:n:h" opt; do
    case $opt in
        g) BUILD_ONLY=yes ;;
        f) FORCE=yes ;;
        r) RESTART=yes ;;
        s) STATUS=yes ;;
        F) FINISH=yes ;;
        t) FTIME=$OPTARG ;;
        a) ARR=$OPTARG ;;
        n) NAME=$OPTARG ;;
        h) usage ;;
        *) usage 1 ;;
    esac
done

. "$CONFIG"
[ -n "$NAME" ] && CASE=$NAME
cd "$(dirname "$CONFIG")"

# ---------------------------------------------------------------- status ---
if [ "$STATUS" = yes ]; then
    [ -d "$CASE" ] || fail "$CASE has not been built"
    cd "$CASE"
    echo "$CASE"
    grep -E "^Time =" log.foamRun 2>/dev/null | tail -1 | sed 's/^/  /'
    grep "^Totals" log.foamRun 2>/dev/null | tail -1 | sed 's/^/  /'
    F=$(find . -name fronts.dat 2>/dev/null | sort | tail -1)
    [ -n "$F" ] && tail -1 "$F" | awk '{
        printf "  front x = %.4f m   speed in this frame = %+.1f m/s   (0 means it is holding station)\n", $2, $3 }'
    ls postProcessing/sootFoil 2>/dev/null | sort -g | tail -1 \
        | sed 's|^|  latest soot foil: postProcessing/sootFoil/|'
    exit 0
fi

# ---------------------------------------------------------------- finish ---
if [ "$FINISH" = yes ]; then
    [ -d "$CASE" ] || fail "$CASE has not been built"
    CASEDIR=$(cd "$CASE" && pwd)
    cd "$CASE"

    if [ -z "$FTIME" ]; then
        FTIME=$(ls -d processor0/[0-9]* processors*/[0-9]* 2>/dev/null \
                | sed 's|.*/||' | sort -g | tail -1)
        [ -n "$FTIME" ] || FTIME=$(ls -d [0-9]*e-0* 2>/dev/null | sort -g | tail -1)
    fi
    [ -n "$FTIME" ] || fail "no written times found in $CASE"

    echo "reconstructing $CASE at t = $FTIME"
    reconstructPar -time "$FTIME" -fileHandler "$FILE_HANDLER" \
        > log.reconstruct 2>&1 || {
            tail -10 log.reconstruct >&2; fail "reconstructPar failed"; }

    # the front position at that time, from the tracking the solver wrote
    F=$(find . -name fronts.dat 2>/dev/null | sort | tail -1)
    [ -n "$F" ] || fail "no fronts.dat: cannot measure the front"
    FRONT_MM=$(awk -v t="$FTIME" 'BEGIN { best = 1e30 }
        !/^#/ && NF > 2 { d = $1 - t; if (d < 0) d = -d
                          if (d < best) { best = d; x = $2 } }
        END { printf "%.4f", x*1000 }' "$F")

    DCJ_USED=$(grep -m1 "D_CJ" log.CDM_ZND 2>/dev/null | awk '{print $3}')
    [ -n "$DCJ_USED" ] || DCJ_USED=$DCJ

    cat << EOF

the snapshot to transplant

  SOURCE_CASE=$CASEDIR
  SOURCE_TIME=$FTIME
  SOURCE_FRONT=$FRONT_MM
  DCJ=$DCJ_USED

check it is quasi-steady before using it: the cells should have stopped
growing and the frame speed should hover about zero
  ./runSource.sh $(basename "$CONFIG") -s
  paraview $CASEDIR/postProcessing/sootFoil/$FTIME/sootFoil.vtk
EOF

    if [ -n "$ARR" ]; then
        cd ..
        [ -f "$ARR" ] || fail "no such arrestor configuration: $ARR"
        cp "$ARR" "$ARR.bak"
        sed -i \
            -e "s|^SOURCE_CASE=.*|SOURCE_CASE=$CASEDIR|" \
            -e "s|^SOURCE_TIME=.*|SOURCE_TIME=$FTIME|" \
            -e "s|^SOURCE_FRONT=.*|SOURCE_FRONT=$FRONT_MM|" \
            -e "s|^DCJ=.*|DCJ=$DCJ_USED|" \
            -e "s|^CDM_A=.*|CDM_A=$CDM_A|" \
            -e "s|^CDM_EaByRT0=.*|CDM_EaByRT0=$CDM_EaByRT0|" \
            -e "s|^CDM_qMwByRT0=.*|CDM_qMwByRT0=$CDM_qMwByRT0|" \
            -e "s|^CDM_kappa0=.*|CDM_kappa0=$CDM_kappa0|" \
            -e "s|^CDM_T0=.*|CDM_T0=$CDM_T0|" \
            -e "s|^CDM_GAMMA=.*|CDM_GAMMA=$CDM_GAMMA|" \
            -e "s|^CDM_MW=.*|CDM_MW=$CDM_MW|" \
            -e "s|^CDM_CV=.*|CDM_CV=$CDM_CV|" \
            -e "s|^P0=.*|P0=$P0|" \
            "$ARR"
        echo
        echo "$ARR updated (previous copy: $ARR.bak).  Run the geometry with"
        echo "    $SCRIPT_DIR/runArrestor.sh $ARR"
    fi
    exit 0
fi

# ---------------------------------------------------------------- checks ---
command -v foamRun > /dev/null || fail "OpenFOAM is not in the environment"
command -v setZNDProfile > /dev/null \
    || fail "setZNDProfile not found: build CDMGpuFoam first (./Allwmake)"

NX=$(echo "$DECOMPOSE" | cut -d' ' -f1)
NY=$(echo "$DECOMPOSE" | cut -d' ' -f2)
[ $((NX*NY)) -eq "$NPROCS" ] \
    || fail "DECOMPOSE ($DECOMPOSE) does not multiply to NPROCS ($NPROCS)"

if [ "$RESTART" = no ] && [ -d "$CASE" ]; then
    [ "$FORCE" = yes ] || fail "$CASE exists (use -f to overwrite, -r to restart)"
    rm -rf "$CASE"
fi

# ---------------------------------------------------------------- build ----
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
    # the CJ state and the ZND profile follow from the model parameters
    TMP=$(mktemp -d)
    python3 "$SCRIPT_DIR/CDM_ZND.py" --paper --A "$CDM_A" \
        --EaByRT0 "$CDM_EaByRT0" --qMwByRT0 "$CDM_qMwByRT0" \
        --kappa0 "$CDM_kappa0" --gamma "$CDM_GAMMA" --Mw "$CDM_MW" \
        --T0 "$CDM_T0" --p0 "$P0" --write "$TMP/znd.dat" \
        > "$TMP/log.CDM_ZND" 2>&1 \
        || { tail -20 "$TMP/log.CDM_ZND" >&2; rm -rf "$TMP"
             fail "CDM_ZND.py failed"; }

    DCJ=$(awk '/D_CJ/ { print $3; exit }' "$TMP/log.CDM_ZND")
    XD=$(awk '/x_d  / { printf "%.6g", $3*10; exit }' "$TMP/log.CDM_ZND")
    [ -n "$DCJ" ] || { rm -rf "$TMP"; fail "no D_CJ in the CDM_ZND output"; }
    [ -n "$XD" ] && XDOPT="--xd $XD" || XDOPT=
    echo "  CDM_ZND: D_CJ = $DCJ m/s, x_d = ${XD:-unknown} mm"

    python3 "$SCRIPT_DIR/writeSourceCase.py" --case "$CASE" \
        --H "$H" --L "$L" --thickness "$THICKNESS" \
        --shockPosition "$SHOCK_POSITION" \
        --dx "$DX" --levels "$LEVELS" --maxCells "$MAXCELLS" \
        --nWaves "$N_WAVES" --amplitude "$AMPLITUDE" \
        --DCJ "$DCJ" $XDOPT \
        --A "$CDM_A" --EaByRT0 "$CDM_EaByRT0" --qMwByRT0 "$CDM_qMwByRT0" \
        --kappa0 "$CDM_kappa0" --T0 "$CDM_T0" --Mw "$CDM_MW" \
        --Cv "$CDM_CV" --p0 "$P0" \
        --endTime "$END_TIME" --writeInterval "$WRITE_INTERVAL" \
        --maxCo "$MAX_CO" --fluxScheme "$FLUX_SCHEME" --nProcs "$NPROCS" \
        --decompose "$DECOMPOSE" --distributor "$DISTRIBUTOR" \
        --fileHandler "$FILE_HANDLER" --foilDelta "$FOIL_DELTA"

    cp "$TMP/znd.dat" "$TMP/log.CDM_ZND" "$CASE/"
    rm -rf "$TMP"
    cp "$CONFIG" "$CASE/case.cfg"
    cd "$CASE"

    run blockMesh
    grep -E "^  nCells" log.blockMesh || true
    run setZNDProfile
    [ "$NPROCS" -gt 1 ] && run decomposePar -fileHandler uncollated
else
    cd "$CASE"
fi

if [ "$BUILD_ONLY" = yes ]; then
    echo "built $CASE (not run: -g)"
    exit 0
fi

# ---------------------------------------------------------------- run ------
if [ "$NPROCS" -gt 1 ]; then
    echo "  foamRun on $NPROCS processors"
    nohup mpirun -np "$NPROCS" --bind-to core foamRun -parallel \
        -fileHandler "$FILE_HANDLER" > log.foamRun 2>&1 &
else
    nohup foamRun -fileHandler "$FILE_HANDLER" > log.foamRun 2>&1 &
fi

cat << EOF

started $CASE (wave-fixed detonation)

  progress    $SCRIPT_DIR/runSource.sh $(basename "$CONFIG") -s
  soot foil   $CASE/postProcessing/sootFoil/<time>/sootFoil.vtk
  when quasi-steady, hand it to a geometry:
              $SCRIPT_DIR/runSource.sh $(basename "$CONFIG") -F -a <arrestor.cfg>
EOF
