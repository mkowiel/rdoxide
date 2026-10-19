import pickle

import pytest
from rdoxide import Chem
from rdoxide.Chem.Scaffolds import rdScaffoldNetwork


def test_brics_scaffolds_and_options():
    m = Chem.MolFromSmiles("c1ccccc1C(=O)NC1NC(=O)CCC1")
    ps = rdScaffoldNetwork.BRICSScaffoldParams()
    ps.includeScaffoldsWithoutAttachments = False
    ps.includeGenericScaffolds = False
    net = rdScaffoldNetwork.CreateScaffoldNetwork([m], ps)
    assert (len(net.nodes), len(net.edges)) == (10, 20)
    assert len(net.molCounts) == len(net.nodes)
    with pytest.raises(AttributeError):
        ps.notAnOption = True


def test_bad_inputs():
    ps = rdScaffoldNetwork.ScaffoldNetworkParams()
    with pytest.raises(ValueError):
        rdScaffoldNetwork.CreateScaffoldNetwork([Chem.MolFromSmiles("CC"), None], ps)
    with pytest.raises(Exception):
        rdScaffoldNetwork.ScaffoldNetworkParams(["[C:1]>>[C:1]"])
    net = rdScaffoldNetwork.CreateScaffoldNetwork([Chem.MolFromSmiles("c1ccccc1C")], ps)
    assert list(pickle.loads(pickle.dumps(net)).nodes) == list(net.nodes)
