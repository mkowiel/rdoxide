import pickle

from rdoxide import Chem
from rdoxide.Chem.Scaffolds import rdScaffoldNetwork


def test_pickle_round_trip_and_empty():
    ms = [Chem.MolFromSmiles(s) for s in ("c1ccccc1CC1NC(=O)CCC1", "c1cccnc1CC1NC(=O)CCC1")]
    net = rdScaffoldNetwork.CreateScaffoldNetwork(ms, rdScaffoldNetwork.ScaffoldNetworkParams())
    net2 = pickle.loads(pickle.dumps(net))
    assert list(net2.nodes) == list(net.nodes)
    assert list(net2.counts) == list(net.counts)
    assert list(net2.molCounts) == list(net.molCounts)
    assert [str(e) for e in net2.edges] == [str(e) for e in net.edges]
    empty = pickle.loads(pickle.dumps(rdScaffoldNetwork.ScaffoldNetwork()))
    assert len(empty.nodes) == 0 and len(empty.edges) == 0
