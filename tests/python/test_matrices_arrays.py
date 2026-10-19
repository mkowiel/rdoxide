import numpy as np

from rdoxide import Chem


def test_matrix_arrays_are_writable_with_expected_dtypes():
    m = Chem.MolFromSmiles("CCOC")
    d = Chem.GetDistanceMatrix(m)
    assert d.dtype == np.float64 and d.shape == (4, 4) and d.flags.writeable
    assert d[0, 3] == 3.0
    a = Chem.GetAdjacencyMatrix(m)
    assert a.dtype == np.intc and a.sum() == 6
    d[0, 0] = 5.0  # writing must not corrupt the cached matrix
    assert Chem.GetDistanceMatrix(m)[0, 0] == 0.0
    assert Chem.GetAdjacencyMatrix(m, useBO=True).dtype == np.float64
