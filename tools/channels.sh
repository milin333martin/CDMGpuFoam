#!/bin/sh
#------------------------------------------------------------------------------
# channels.sh <case> [x0] [x1]
#
# Penetration into each channel of an arrestor at the latest sample of the
# front tracking.  x0 and x1 are the arrestor's start and end in metres
# (default 0.12 and 0.13); a channel whose shock has passed x1 is marked
# THROUGH.
#------------------------------------------------------------------------------
[ $# -ge 1 ] || { sed -n '3,9p' "$0"; exit 1; }

F=$(find "$1" -name fronts.dat 2>/dev/null | sort | tail -1)
[ -n "$F" ] || { echo "no fronts.dat under $1" >&2; exit 1; }

X0=${2:-0.12}
X1=${3:-0.13}

N=$(( ($(tail -1 "$F" | wc -w) - 1) / 5 - 1 ))

tail -1 "$F" | awk -v n="$N" -v x0="$X0" -v x1="$X1" '{
  printf "t = %.1f us   global: shock %.4f  flame %.4f  separation %.3f mm\n",
         $1*1e6, $2, $4, $6*1000
  inside = 0
  for (i = 1; i <= n; i++) {
    xs = $(2+5*i); xf = $(4+5*i)
    if (xs > x0) {
      inside++
      printf "  ch%-3d shock %.4f (%5.2f mm in)  flame %.4f  sep %.3f mm%s\n",
             i, xs, (xs-x0)*1000, xf, (xs-xf)*1000, (xs > x1 ? "   THROUGH" : "")
    }
  }
  if (!inside) print "  no channel reached yet"
}'
