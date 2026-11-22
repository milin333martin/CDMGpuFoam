#!/bin/sh
#------------------------------------------------------------------------------
# runPipeline.sh <source.cfg> <geometry.cfg> [geometry2.cfg ...]
#
# The whole study, start to finish, one case at a time:
#
#   1. run the wave-fixed (shock-frame) detonation until the cells settle
#   2. reconstruct its final time and measure the front
#   3. write that snapshot into each geometry configuration
#   4. build and run each geometry case in turn (arrestor, narrow channel)
#
# Each stage waits for the previous one, so the cores are never shared.  The
# pipeline can be stopped and started again: a source case that has already
# finished is reused rather than rerun.
#
#     ./runPipeline.sh source-T6.cfg arrestor-BR09.cfg arrestor-BR09-CHT.cfg
#     ./runPipeline.sh source-T6.cfg arrestor-BR09.cfg -t 8e-05   # snapshot time
#------------------------------------------------------------------------------
set -e
SCRIPT_DIR=$(cd "$(dirname "$0")" && pwd)

usage() { sed -n '3,18p' "$0"; exit ${1:-0}; }
fail()  { echo "ERROR: $*" >&2; exit 1; }

[ $# -ge 2 ] || usage 1

SRC=$1; shift
GEOM=""
TIME_OPT=""
while [ $# -gt 0 ]; do
    case $1 in
        -t) TIME_OPT="-t $2"; shift 2 ;;
        -h) usage ;;
        *)  GEOM="$GEOM $1"; shift ;;
    esac
done
[ -n "$GEOM" ] || fail "no geometry configuration given"

for f in "$SRC" $GEOM; do
    [ -f "$f" ] || fail "no such configuration: $f"
done

wait_for_run()      # wait until no solver of ours is running
{
    while pgrep -f "foamRun -parallel|foamMultiRun -parallel" > /dev/null \
       || pgrep -x foamRun > /dev/null || pgrep -x foamMultiRun > /dev/null
    do
        sleep 30
    done
}

finished()          # finished <config>: has this case reached its endTime?
{
    c=$(. "./$1" > /dev/null 2>&1; echo "$CASE")
    e=$(. "./$1" > /dev/null 2>&1; echo "$END_TIME")
    [ -d "$c" ] || return 1
    grep -q "^End$" "$c"/log.foam*Run 2>/dev/null && return 0
    t=$(grep -hE "^ *Time = " "$c"/log.foam*Run 2>/dev/null | tail -1 \
        | sed 's/.*= *//; s/s$//')
    [ -n "$t" ] || return 1
    awk -v t="$t" -v e="$e" 'BEGIN { exit !(t >= 0.999*e) }'
}

echo "==== 1. wave-fixed detonation: $SRC"
if finished "$SRC"; then
    echo "     already finished, reusing it"
else
    "$SCRIPT_DIR/runSource.sh" "$SRC" || fail "the source run did not start"
    echo "     waiting for it to finish (progress: runSource.sh $SRC -s)"
    sleep 30
    wait_for_run
    finished "$SRC" || fail "the source run stopped early: check its log"
fi

echo
echo "==== 2. snapshot"
for g in $GEOM; do
    "$SCRIPT_DIR/runSource.sh" "$SRC" -F $TIME_OPT -a "$g" | sed 's/^/     /'
done

for g in $GEOM; do
    echo
    echo "==== 3. geometry: $g"
    if finished "$g"; then
        echo "     already finished, skipping"
        continue
    fi
    "$SCRIPT_DIR/runArrestor.sh" "$g" -f || fail "$g did not start"
    echo "     waiting for it to finish"
    sleep 30
    wait_for_run
done

echo
echo "==== done"
for g in $GEOM; do
    c=$(. "./$g" > /dev/null 2>&1; echo "$CASE")
    echo "  $c"
    "$SCRIPT_DIR/channels.sh" "$c" 2>/dev/null | sed 's/^/      /' || true
done
echo
echo "  compare them:  $SCRIPT_DIR/plotFronts.py <case> <case> --save compare.png"
