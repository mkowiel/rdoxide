import pytest
from rdkit import DataStructs
from rdkit.Chem import rdFingerprintGenerator, rdMolProcessing


def _write(tmp_path, name, text):
    path = tmp_path / name
    path.write_text(text)
    return str(path)


def test_results_do_not_depend_on_thread_count(tmp_path):
    smis = ["CCO", "c1ccccc1", "CC(=O)O", "CCN", "C1CC1"] * 300
    path = _write(tmp_path, "mols.smi", "".join(f"{s}\tm{i}\n" for i, s in enumerate(smis)))
    opts = rdMolProcessing.SupplierOptions()
    opts.titleLine = False
    ref = None
    for n in (1, 2, 4):
        opts.numThreads = n
        fps = rdMolProcessing.GetFingerprintsForMolsInFile(path, options=opts)
        assert len(fps) == len(smis)
        bits = [tuple(fp.GetOnBits()) for fp in fps]
        if ref is None:
            ref = bits
        assert bits == ref


def test_unreadable_molecules_and_generator(tmp_path):
    path = _write(tmp_path, "mols.smi", "CCO a\nnot_a_smiles b\nCCN c\nbad d\n")
    opts = rdMolProcessing.SupplierOptions()
    opts.titleLine = False
    opts.delimiter = " "
    fpg = rdFingerprintGenerator.GetRDKitFPGenerator(fpSize=512)
    fps = rdMolProcessing.GetFingerprintsForMolsInFile(path, generator=fpg, options=opts)
    assert len(fps) == 3  # trailing unreadable entries are dropped
    assert fps[1] is None
    assert fps[0].GetNumBits() == 512
    assert DataStructs.TanimotoSimilarity(fps[0], fps[2]) > 0


def test_bad_files_and_options(tmp_path):
    with pytest.raises(OSError):
        rdMolProcessing.GetFingerprintsForMolsInFile(str(tmp_path / "missing.smi"))
    with pytest.raises(OSError):
        rdMolProcessing.GetFingerprintsForMolsInFile(_write(tmp_path, "x.xyz", "CCO\n"))
    with pytest.raises(OSError):
        rdMolProcessing.GetFingerprintsForMolsInFile(_write(tmp_path, "empty.smi", ""))
    opts = rdMolProcessing.SupplierOptions()
    with pytest.raises(AttributeError):
        opts.notAnOption = 1
    opts.numThreads = 3
    assert opts.numThreads == 3
