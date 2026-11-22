#!/usr/bin/env python3
"""Add the empty 'internal' patch required by OpenFOAM-13 mesh
redistribution (fvMeshDistribute) to a case: blockMeshDict and 0/ fields."""
import sys, os, re
case = sys.argv[1] if len(sys.argv) > 1 else "."
bm = os.path.join(case, "system/blockMeshDict")
s = open(bm).read()
if "type internal;" not in s:
    s = re.sub(r"(\nboundary\s*\n\()",
               r"\1\n    // Empty patch required by OpenFOAM-13 for dynamic load balancing"
               r"\n    internalFaces { type internal; faces (); }\n", s, count=1)
    open(bm, "w").write(s)
    print("patched", bm)
for f in sorted(os.listdir(os.path.join(case, "0"))):
    p = os.path.join(case, "0", f)
    t = open(p).read()
    if "boundaryField" in t and "type            internal;" not in t:
        t = re.sub(r"(boundaryField\s*\n\{\n)",
                   r"\1    internalFaces\n    {\n        type            internal;\n    }\n", t, count=1)
        open(p, "w").write(t)
        print("patched", p)
