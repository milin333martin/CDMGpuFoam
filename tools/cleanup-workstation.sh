#!/bin/sh
#------------------------------------------------------------------------------
# cleanup-workstation.sh
#
# Tidies an OpenFOAM workstation after a CDMGpuFoam study:
#
#   1. rescues the wave-fixed detonation from source-T6, reducing it to a
#      transplant snapshot (reconstructed final time + soot foils + the
#      configuration) and reporting the numbers a future transplant needs;
#   2. moves that snapshot OUT of the package tree, because clean-workspace.sh
#      deletes any directory holding a case.cfg;
#   3. cleans the package tree: build products, logs, VTK, run output.
#
# Nothing is deleted until you pass --apply.  Run it without that first and
# read what it says it would do.
#
#     ./cleanup-workstation.sh                 # dry run: report only
#     ./cleanup-workstation.sh --apply         # actually do it
#     ./cleanup-workstation.sh --apply --no-snapshot   # skip step 1 and 2
#------------------------------------------------------------------------------
set -e

TREE=${CDMREACTFOAM_TREE:-$HOME/OpenFOAM/CDMGpuFoam-1.4.0}
CASE=${SOURCE_CASE_DIR:-$TREE/examples/source-T6}
KEEP=${SNAPSHOT_DIR:-$HOME/OpenFOAM/snapshots}

APPLY=no
SNAPSHOT=yes
for a in "$@"
do
    case $a in
        --apply)       APPLY=yes ;;
        --no-snapshot) SNAPSHOT=no ;;
        -h|--help)     sed -n '3,21p' "$0"; exit 0 ;;
        *) echo "unknown option: $a" >&2; exit 1 ;;
    esac
done

say()  { printf '%s\n' "$*"; }
step() { printf '\n== %s\n' "$*"; }
size() { du -sh "$1" 2>/dev/null | cut -f1; }

run()
{
    if [ "$APPLY" = yes ]
    then
        say "   + $*"
        "$@"
    else
        say "   would run: $*"
    fi
}

# A guard: never operate on something that is not a CDMGpuFoam tree.
[ -f "$TREE/Allwmake" ] && [ -d "$TREE/src/CDMModels" ] || {
    echo "ERROR: $TREE does not look like a CDMGpuFoam tree." >&2
    echo "       Set CDMREACTFOAM_TREE to the right path." >&2
    exit 1; }

[ -n "$WM_PROJECT_DIR" ] || {
    echo "ERROR: source the OpenFOAM-13 environment first." >&2; exit 1; }

say "tree      $TREE          ($(size "$TREE"))"
say "case      $CASE          ($(size "$CASE"))"
say "snapshots $KEEP"
[ "$APPLY" = yes ] || say "
DRY RUN - nothing will be deleted.  Re-run with --apply when the plan below
looks right."

#--- 1. rescue the detonation -------------------------------------------------
if [ "$SNAPSHOT" = yes ] && [ -d "$CASE" ]
then
    step "1. reduce $CASE to a transplant snapshot"

    LATEST=$(ls -d "$CASE"/processor0/[0-9]* "$CASE"/processors*/[0-9]* \
             2>/dev/null | sed 's|.*/||' | sort -g | tail -1)

    if [ -z "$LATEST" ]
    then
        say "   no written times found in the processor directories;"
        say "   nothing to rescue, skipping to the clean."
        SNAPSHOT=no
    else
        say "   latest written time: $LATEST"
        say "   processor data:      $(du -sh "$CASE"/processor* 2>/dev/null \
                                       | tail -1 | cut -f1) total"

        # reconstructPar needs the handler the run actually used
        if ls -d "$CASE"/processors* > /dev/null 2>&1
        then HANDLER=collated
        else HANDLER=uncollated
        fi
        say "   file handler:        $HANDLER"

        if [ "$APPLY" = yes ]
        then
            ( cd "$CASE" && \
              reconstructPar -time "$LATEST" -fileHandler "$HANDLER" \
                  > log.reconstructSnapshot 2>&1 ) || {
                echo "   reconstructPar FAILED - see $CASE/log.reconstructSnapshot" >&2
                echo "   stopping here; nothing has been deleted." >&2
                exit 1; }

            # the fields must actually be there, not just a polyMesh
            if [ ! -f "$CASE/$LATEST/p" ]
            then
                echo "   $CASE/$LATEST holds no fields - reconstruction did not" >&2
                echo "   work.  Stopping; nothing has been deleted." >&2
                exit 1
            fi
            say "   reconstructed $LATEST ($(size "$CASE/$LATEST"))"
        else
            say "   would run: reconstructPar -time $LATEST -fileHandler $HANDLER"
        fi

        # the numbers a future transplant needs
        step "   the snapshot's transplant parameters"
        if [ "$APPLY" = yes ]
        then
            FRONTS=$(find "$CASE" -name fronts.dat 2>/dev/null | sort | tail -1)
            FRONT_MM=""
            [ -n "$FRONTS" ] && FRONT_MM=$(awk -v t="$LATEST" 'BEGIN { best=1e30 }
                !/^#/ && NF > 2 { d=$1-t; if (d<0) d=-d
                                  if (d<best) { best=d; x=$2 } }
                END { if (x != "") printf "%.4f", x*1000 }' "$FRONTS")
            DCJ=$(grep -m1 "D_CJ" "$CASE/log.CDM_ZND" 2>/dev/null | awk '{print $3}')

            {
                echo "# transplant parameters for $CASE"
                echo "# written by cleanup-workstation.sh on $(date)"
                echo "SOURCE_CASE=$KEEP/$(basename "$CASE")"
                echo "SOURCE_TIME=$LATEST"
                echo "SOURCE_FRONT=${FRONT_MM:-UNKNOWN}"
                echo "DCJ=${DCJ:-UNKNOWN}"
            } > "$CASE/SNAPSHOT.txt"
            cat "$CASE/SNAPSHOT.txt" | sed 's/^/   /'
            [ -n "$FRONT_MM" ] || say "
   NOTE: no fronts.dat, so the front position could not be measured.  You
   will have to measure it by hand before using this as a transplant source."
        else
            say "   would measure the front from fronts.dat and write SNAPSHOT.txt"
        fi

        #--- 2. move it out of the package tree ------------------------------
        step "2. move the snapshot out of the package tree"
        say "   (clean-workspace.sh deletes any directory holding a case.cfg,"
        say "    so the snapshot cannot stay under $TREE)"
        run mkdir -p "$KEEP"
        run mv "$CASE" "$KEEP/"

        TARGET=$KEEP/$(basename "$CASE")
        say "   dropping the processor data from the snapshot:"
        if [ "$APPLY" = yes ]
        then
            rm -rf "$TARGET"/processor[0-9]* "$TARGET"/processors[0-9]*
            say "   snapshot is now $(size "$TARGET")"
        else
            say "   would run: rm -rf $TARGET/processor*"
        fi
    fi
fi

#--- 3. clean the package tree ------------------------------------------------
step "3. clean the package tree"

say "   build products (Allwclean):"
if [ "$APPLY" = yes ]
then
    ( cd "$TREE" && ./Allwclean > /dev/null 2>&1 ) || true
    say "   + done"
else
    say "   would run: (cd $TREE && ./Allwclean)"
fi

say "   run output, logs, VTK, leftover cases (clean-workspace.sh):"
if [ -x "$TREE/clean-workspace.sh" ]
then
    if [ "$APPLY" = yes ]
    then
        "$TREE/clean-workspace.sh" -y "$TREE"
    else
        "$TREE/clean-workspace.sh" "$TREE" || true
    fi
else
    say "   clean-workspace.sh not found; cleaning by hand"
    run find "$TREE" -type d \( -name "processor[0-9]*" -o -name "processors[0-9]*" \
         -o -name postProcessing -o -name VTK -o -name Alltest.tmp \) -prune -exec rm -rf {} +
    run find "$TREE" -type f \( -name "log.*" -o -name "*.bak" -o -name "nohup.out" \) -delete
fi

#--- 4. report ----------------------------------------------------------------
step "4. result"
say "   package tree  $TREE  -> $(size "$TREE")"
[ -d "$KEEP" ] && say "   snapshots     $KEEP  -> $(size "$KEEP")"
df -h "$HOME" | tail -1 | sed 's/^/   /'

if [ "$APPLY" != yes ]
then
    say "
Nothing was changed.  Re-run with --apply to carry the plan out."
else
    say "
Done.  To share the code, send the archive that ./package.sh makes rather
than repackaging this tree by hand: the release carries the manual, the
build-environment check and the tutorial fixes, and nothing from a run."
fi
