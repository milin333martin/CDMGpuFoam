#!/usr/bin/env python3
"""
addInternalPatchToMesh.py <case> [region]

Adds the empty 'internalFaces' patch of type 'internal' to a mesh that has
already been written, which splitMeshRegions removes along with the other
empty patches ("Deleting empty patches").  OpenFOAM-13's mesh redistribution
(fvMeshDistribute, used by the load balancer) requires it, and without it a
parallel run stops with

    Cannot find a internal patch in ...

Safe to run twice: it does nothing if the patch is already there.
"""

import os
import re
import sys

case = sys.argv[1] if len(sys.argv) > 1 else "."
region = sys.argv[2] if len(sys.argv) > 2 else ""

b = os.path.join(case, "constant", region, "polyMesh", "boundary")
if not os.path.isfile(b):
    raise SystemExit("no boundary file: " + b)

s = open(b).read()

if re.search(r"type\s+internal\s*;", s):
    print("already present:", b)
    raise SystemExit(0)

# the end of the face list: the largest startFace + nFaces of any patch
ends = [int(sf) + int(nf) for nf, sf in
        re.findall(r"nFaces\s+(\d+);\s*startFace\s+(\d+);", s)]
if not ends:
    raise SystemExit("no patches found in " + b)
end = max(ends)

entry = """    internalFaces
    {
        type            internal;
        nFaces          0;
        startFace       %d;
    }
""" % end

# insert before the final closing bracket of the patch list
i = s.rindex(")")
s = s[:i] + entry + s[i:]

# The patch count is the FIRST "<n>" on its own line followed by "(".
# A mapped patch writes its inGroups the same way, so a later match would
# be that list rather than the patch count.
m = re.search(r"^[ \t]*(\d+)[ \t]*\r?\n[ \t]*\(", s, re.M)
if not m:
    raise SystemExit("could not find the patch count in " + b)
s = s[:m.start(1)] + str(int(m.group(1)) + 1) + s[m.end(1):]

open(b, "w").write(s)
print("added internalFaces to", b, "(startFace %d)" % end)
