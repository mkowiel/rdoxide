"""rdMolStandardize bindings (the behaviour is checked by the upstream
testMolStandardize.py and the `standardize` differential)."""
import pytest
from rdoxide import Chem
from rdoxide.Chem.MolStandardize import ReorderTautomers
from rdoxide.Chem.MolStandardize import rdMolStandardize as ms


def smi(m):
    return Chem.MolToSmiles(m)


def test_module_and_classes():
    assert ms.__name__ == "rdoxide.Chem.MolStandardize.rdMolStandardize"
    for cls in (ms.CleanupParameters, ms.TautomerEnumerator, ms.Uncharger, ms.Pipeline):
        assert cls.__module__ == "rdoxide.Chem.MolStandardize.rdMolStandardize"
    assert issubclass(ms.RDKitValidation, ms.ValidationMethod)
    assert ms.TautomerEnumerator.tautomerScoreVersion == "1.0.0"


def test_cleanup_parameters():
    p = ms.CleanupParameters()
    assert (p.maxRestarts, p.maxTautomers, p.preferOrganic, p.doCanonical) == (200, 1000, False, True)
    ms.UpdateParamsFromJSON(p, '{"maxTautomers": 7, "tautomerRemoveSp3Stereo": false}')
    assert p.maxTautomers == 7 and not p.tautomerRemoveSp3Stereo
    with pytest.raises(AttributeError, match="Cannot set unknown attribute 'foo'"):
        p.foo = 1


def test_functions_and_none():
    m = Chem.MolFromSmiles("CCC(=O)O[Na]")
    before = smi(m)
    assert smi(ms.Cleanup(m)) == "CCC(=O)[O-].[Na+]"
    assert smi(m) == before  # the input is not modified
    with pytest.raises(ValueError, match="Molecule is None"):
        ms.Cleanup(None)
    assert smi(ms.ChargeParent(Chem.MolFromSmiles("C[NH+](C)C.[Cl-]"))) == "CN(C)C"
    assert smi(ms.FragmentParent(m, skipStandardize=True)) == before
    assert ms.StandardizeSmiles("CCC(=O)O[Na]") == "CCC(=O)[O-].[Na+]"
    with pytest.raises(ValueError):
        ms.StandardizeSmiles("C1CC")


def test_in_place_overloads():
    m = Chem.MolFromSmiles("CCC(=O)O[Na]")
    ms.CleanupInPlace(m)
    assert smi(m) == "CCC(=O)[O-].[Na+]"
    mols = [Chem.MolFromSmiles(s) for s in ("CCC(=O)O[Na]", "C[NH+](C)C.[Cl-]", "CN(=O)=O")]
    ms.ChargeParentInPlace(mols, 2)
    assert [smi(x) for x in mols] == ["CCC(=O)O", "CN(C)C", "C[N+](=O)[O-]"]
    ms.NormalizeInPlace(mols=mols, numThreads=1)
    with pytest.raises(ValueError, match="duplicate molecule"):
        ms.CleanupInPlace([mols[0], mols[0]], 1)


def test_components():
    assert smi(ms.Uncharger(canonicalOrder=False).uncharge(Chem.MolFromSmiles("[O-]c1ccccc1"))) == "Oc1ccccc1"
    assert smi(ms.LargestFragmentChooser(preferOrganic=True).choose(Chem.MolFromSmiles("CC.[Na+]"))) == "CC"
    assert smi(ms.FragmentRemoverFromData("sodium\t[Na+]\n").remove(Chem.MolFromSmiles("CC.[Na+]"))) == "CC"
    r = ms.ReionizerFromData("carboxylic acid\tC(=O)[OH]\tC(=O)[O-]\n", [ms.ChargeCorrection("Na", "[Na;X0+0]", 1)])
    assert smi(r.reionize(Chem.MolFromSmiles("CC(=O)O.[Na]"))) == "CC(=O)[O-].[Na+]"
    assert [c.Name for c in ms.CHARGE_CORRECTIONS()] == ["[Li,Na,K]", "[Mg,Ca]", "[Cl]"]
    md = ms.MetalDisconnector()
    assert md.MetalNof.startswith("[#3,#11,#19,")
    opts = ms.MetalDisconnectorOptions()
    opts.splitGrignards = True
    assert smi(ms.MetalDisconnector(opts).Disconnect(Chem.MolFromSmiles("[Mg](Cl)CC"))) == "[CH2-]C.[Cl-].[Mg+2]"


def test_validation():
    assert ms.ValidateSmiles("ClCCCl.c1ccccc1O") == ["INFO: [FragmentValidation] 1,2-dichloroethane is present"]
    v = ms.AllowedAtomsValidation([Chem.Atom(6), Chem.Atom(1)])
    assert v.validate(Chem.MolFromSmiles("CO")) == ["INFO: [AllowedAtomsValidation] Atom O is not in allowedAtoms list"]

    class OnlyCarbon(ms.ValidationMethod):
        def validate(self, mol, reportAllFailures=False):
            return [] if all(a.GetAtomicNum() == 6 for a in mol.GetAtoms()) else ["not carbon"]

    composite = ms.MolVSValidation([ms.NeutralValidation(), OnlyCarbon()])
    assert composite.validate(Chem.MolFromSmiles("CO")) == ["not carbon"]
    assert composite.validate(Chem.MolFromSmiles("C[N+](C)(C)C"), True) == [
        "INFO: [NeutralValidation] Not an overall neutral system (+1)", "not carbon"]


def test_tautomers():
    te = ms.TautomerEnumerator()
    res = te.Enumerate(Chem.MolFromSmiles("Oc1ccccn1"))
    assert list(res.smiles) == ["O=C1CC=CC=N1", "O=c1cccc[nH]1", "Oc1ccccn1"]
    assert len(res) == 3 and smi(res[-1]) == "Oc1ccccn1" and res.modifiedAtoms == (0, 2, 6)
    assert res.status == ms.TautomerEnumeratorStatus.Completed
    assert list(res.smilesTautomerMap.keys()) == list(res.smiles)
    assert smi(res.smilesTautomerMap.values()[0].kekulized) == "O=C1CC=CC=N1"
    assert smi(te.PickCanonical(res)) == "O=c1cccc[nH]1"
    assert smi(te.PickCanonical(set(res()))) == "O=c1cccc[nH]1"
    with pytest.raises(TypeError):
        te.PickCanonical([0, 1])
    with pytest.raises(TypeError):
        te.Canonicalize(Chem.MolFromSmiles("Oc1ccccn1"), lambda m: "fail")
    assert smi(te.Canonicalize(Chem.MolFromSmiles("Oc1ccccn1"), lambda m: m.GetNumAtoms())) in res.smiles
    assert smi(ReorderTautomers(Chem.MolFromSmiles("Oc1ccccn1"))[0]) == "O=c1cccc[nH]1"
    terms = ms.SubstructTermVector()
    terms.append(ms.SubstructTerm("C=O", "[#6]=,:[#8]", 1000))
    assert ms.ScoreSubstructs(Chem.MolFromSmiles("CC(=O)C"), terms) == 1000
    assert len(ms.GetDefaultTautomerScoreSubstructs()) == 12


def test_tautomer_callback():
    class Stop(ms.TautomerEnumeratorCallback):
        def __call__(self, mol, res):
            assert isinstance(res, ms.TautomerEnumeratorResult)
            return len(res) < 2

    class Broken(ms.TautomerEnumeratorCallback):
        pass

    class Raises(ms.TautomerEnumeratorCallback):
        def __call__(self, mol, res):
            raise KeyError("boom")

    te = ms.TautomerEnumerator()
    with pytest.raises(AttributeError):
        te.SetCallback(Broken())
    cb = Stop()
    te.SetCallback(cb)
    assert te.GetCallback() is cb
    res = te.Enumerate(Chem.MolFromSmiles("CC(=O)CC(=O)C"))
    assert res.status == ms.TautomerEnumeratorStatus.Canceled and len(res) == 2
    te.SetCallback(Raises())
    with pytest.raises(KeyError, match="boom"):
        te.Enumerate(Chem.MolFromSmiles("CC(=O)CC(=O)C"))


def test_pipeline():
    m = Chem.MolFromSmiles("CC[O-].[Na+]")
    from rdoxide.Chem import rdDepictor
    rdDepictor.Compute2DCoords(m)
    res = ms.Pipeline().run(Chem.MolToMolBlock(m))
    assert res.stage == ms.PipelineStage.COMPLETED
    assert res.status & ms.PipelineStatus.FRAGMENTS_REMOVED
    assert any(e.detail.startswith("One or more disconnected fragments") for e in res.log)
    bad = ms.Pipeline().run("garbage")
    assert bad.status & ms.PipelineStatus.INPUT_ERROR


def test_vector_suites():
    v = ms.GetDefaultTautomerScoreSubstructs()
    assert [t.name for t in v[1:3]] == ["oxim", "C=O"] and type(v[1:3]) is ms.SubstructTermVector
    assert ms.SubstructTerm("C=O", "[#6]=,:[#8]", 2) in v and 5 not in v
    with pytest.raises(IndexError, match="Index out of range"):
        v[100]
    with pytest.raises(IndexError, match="slice step size not supported."):
        v[::2]
    with pytest.raises(TypeError, match="Attempting to append an invalid type"):
        v.append(1)
    del v[0:2]
    assert len(v) == 10
    for cls in (ms.PipelineLog, ms.PipelineResult, ms.TautomerEnumeratorResult):
        with pytest.raises(RuntimeError, match="cannot be instantiated"):
            cls()


def test_pipeline_log_is_shared():
    m = Chem.MolFromSmiles("CC[O-].[Na+]")
    from rdoxide.Chem import rdDepictor
    rdDepictor.Compute2DCoords(m)
    res = ms.Pipeline().run(Chem.MolToMolBlock(m))
    n = len(res.log)
    res.log.append(res.log[0])
    assert len(res.log) == n + 1


def _without_rings(smi):
    src = Chem.MolFromSmiles(smi)
    rw = Chem.RWMol()
    for a in src.GetAtoms():
        rw.AddAtom(Chem.Atom(a))
    for b in src.GetBonds():
        rw.AddBond(b.GetBeginAtomIdx(), b.GetEndAtomIdx(), b.GetBondType())
    m = rw.GetMol()
    m.UpdatePropertyCache(False)
    return m


def test_scoring_leaves_rings_on_the_input():
    m = _without_rings("c1ccccc1C1CCCC1")
    with pytest.raises(RuntimeError):
        m.GetRingInfo().NumRings()
    assert ms.ScoreRings(m) == 250
    assert m.GetRingInfo().NumRings() == 2


def test_catalog_files_are_cached(tmp_path):
    path = tmp_path / "acidbase.txt"
    path.write_text("acid\tC(=O)[OH]\tC(=O)[O-]\n")
    p = ms.CleanupParameters()
    p.acidbaseFile = str(path)
    first = smi(ms.Reionize(Chem.MolFromSmiles("OC(=O)CC(=O)[O-]"), p))
    path.write_text("not a valid line\n")
    assert smi(ms.Reionize(Chem.MolFromSmiles("OC(=O)CC(=O)[O-]"), p)) == first
