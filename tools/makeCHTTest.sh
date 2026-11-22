#!/bin/sh
#------------------------------------------------------------------------------
# makeCHTTest.sh: build and run the two-region conjugate heat transfer test.
#
#     ./makeCHTTest.sh            generate, split the regions and run
#     ./makeCHTTest.sh -g         generate only
#
# Hot burnt gas on a steel plate, no flow and no reaction, so that the
# interface temperature can be compared with the analytic contact value the
# script prints.
#------------------------------------------------------------------------------
set -e
SCRIPT_DIR=$(cd "$(dirname "$0")" && pwd)
CASE=chtTest
GEN_ONLY=no

while getopts "gn:" opt; do
    case $opt in
        g) GEN_ONLY=yes ;;
        n) CASE=$OPTARG ;;
        *) exit 1 ;;
    esac
done

python3 "$SCRIPT_DIR/writeCHTTest.py" --case "$CASE" --stage base

cd "$CASE"

run()
{
    echo "$1"
    if ! "$@" > "log.$1" 2>&1; then
        echo "FAILED: $* (see $CASE/log.$1)" >&2
        tail -20 "log.$1" >&2
        exit 1
    fi
}

run blockMesh
run splitMeshRegions -cellZones -overwrite

echo "regions:"
ls constant | grep -v polyMesh | sed 's/^/    /'

cd ..
python3 "$SCRIPT_DIR/writeCHTTest.py" --case "$CASE" --stage fields
cd "$CASE"

[ "$GEN_ONLY" = yes ] && { echo "generated only (-g)"; exit 0; }

run foamMultiRun

echo
echo "interface temperature reached:"
LAST=$(ls -d [0-9]* 0.* 2>/dev/null | sort -g | tail -1)
echo "    last time $LAST"
echo "compare with the EXPECTED value printed above"
