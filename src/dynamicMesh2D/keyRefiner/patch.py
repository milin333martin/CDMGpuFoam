import io, sys
C = 'keyRefiner_fvMeshTopoChanger.C'
H = 'keyRefiner_fvMeshTopoChanger.H'

if 'verifyCentres' in io.open(C).read():
    sys.exit('ALREADY APPLIED')

E = [
 (C,
  '''    label nx = dict_.lookupOrDefault<label>("nx", 0);
    label ny = dict_.lookupOrDefault<label>("ny", 0);

    if (nx <= 0 || ny <= 0)
    {''',
  '''    label nx = dict_.lookupOrDefault<label>("nx", 0);
    label ny = dict_.lookupOrDefault<label>("ny", 0);

    // A refined mesh cannot be counted back to its base, so on a restart
    // nx and ny come from the key file's own header.
    const bool restarting = isFile(keysPath(mesh().time().name()));

    if ((nx <= 0 || ny <= 0) && restarting)
    {
        peekKeysHeader(nx, ny);
    }

    if (nx <= 0 || ny <= 0)
    {'''),

 (C,
  '    if (nx <= 0 || ny <= 0 || nx*ny != m.nCells())',
  '    if (nx <= 0 || ny <= 0 || (!restarting && nx*ny != m.nCells()))'),

 (C,
  '''    scalar worst = 0;
    label worstCell = -1;
    const scalar h = min(Lx/nx, Ly/ny);''',
  '''    // On a restart the case mesh is refined, so the base mesh cannot be
    // compared with it cell for cell here.  verifyCentres(), called once
    // the keys are read, compares the engine's refined mesh against it
    // instead, which is the stronger statement.
    if (restarting)
    {
        Info<< "keyRefiner: base mesh " << nx << " x " << ny
            << ", restarting from a key file" << endl;
        return;
    }

    scalar worst = 0;
    label worstCell = -1;
    const scalar h = min(Lx/nx, Ly/ny);'''),

 (C,
  '''    readKeys();

    seedNames();''',
  '''    readKeys();

    // After the keys, so that a restart compares the engine's refined
    // mesh with the case's refined mesh.
    verifyCentres();

    seedNames();'''),

 (C,
  '            << 1.0e6*total/(nChanges_*max(amr_.nCells(), 1))',
  '            << 1.0e6*total/(scalar(nChanges_)*max(amr_.nCells(), 1))'),

 (H,
  '        void fitBaseMesh();',
  '''        void fitBaseMesh();

        //- nx and ny from the key file's header, without reading the keys
        bool peekKeysHeader(label& nx, label& ny) const;

        //- Check the engine's mesh is the case's mesh cell for cell
        void verifyCentres() const;'''),
]

NEWFUNCS = r'''

bool Foam::fvMeshTopoChangers::keyRefiner::peekKeysHeader
(
    label& nx,
    label& ny
) const
{
    const fileName f(keysPath(mesh().time().name()));

    if (!isFile(f)) { return false; }

    IFstream is(f);
    label n = 0;
    is >> nx >> ny >> n;

    return nx > 0 && ny > 0;
}


void Foam::fvMeshTopoChangers::keyRefiner::verifyCentres() const
{
    // The assumption the whole scheme rests on: the engine's cell i is
    // the fvMesh's cell i.  If it were false every field would be
    // permuted on the first rebuild while the mesh stayed perfectly
    // valid, and nothing downstream could detect it.
    const fvMesh& m = mesh();

    if (amr_.nCells() != m.nCells())
    {
        FatalErrorInFunction
            << "the engine has " << amr_.nCells() << " cells but the mesh"
            << " has " << m.nCells() << "." << nl
            << "  On a restart that means the key file and the mesh in"
            << " the same time directory are not the same mesh."
            << exit(FatalError);
    }

    const scalar h = min(amr_.Lx/amr_.nx, amr_.Ly/amr_.ny);

    scalar worst = 0;
    label worstCell = -1;

    forAll(m.C(), celli)
    {
        const scalar d =
            mag(m.C()[celli].x() - amr_.C[2*celli])
          + mag(m.C()[celli].y() - amr_.C[2*celli + 1]);

        if (d > worst)
        {
            worst = d;
            worstCell = celli;
        }
    }

    if (worst > 1e-8*h)
    {
        FatalErrorInFunction
            << "the key engine's mesh is not the case's mesh cell for"
            << " cell." << nl
            << "  Worst at cell " << worstCell << ": the mesh has it at "
            << m.C()[worstCell] << ", the engine at ("
            << amr_.C[2*worstCell] << " " << amr_.C[2*worstCell + 1]
            << "), off by " << worst << " against a base cell size of "
            << h << "." << exit(FatalError);
    }

    Info<< "    cell centres agree with the engine to " << worst
        << " over " << m.nCells() << " cells (base cell " << h << ")"
        << endl;
}


// ************************************************************************* //
'''

bad = 0
for path, old, new in E:
    s = io.open(path).read()
    if s.count(old) != 1:
        print('FAIL %s x%d: %s' % (path, s.count(old), old.strip()[:50]))
        bad += 1
        continue
    io.open(path, 'w').write(s.replace(old, new, 1))
    print('ok   %s  %s' % (path, old.strip()[:50]))

if not bad:
    io.open(C, 'a').write(NEWFUNCS)
    print('ok   appended peekKeysHeader + verifyCentres')

sys.exit(1 if bad else 0)
