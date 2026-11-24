import os

import numpy as np
import pytest

from scicone import SCICoNE
from scicone.scicone import find_binary


def test_binaries_are_found():
    """The four installed executables must resolve to something runnable."""
    sci = SCICoNE()
    for binary in (sci.simulation_binary, sci.bp_binary,
                   sci.inference_binary, sci.score_binary):
        assert binary is not None
        assert os.access(binary, os.X_OK)


def test_explicit_binary_path_is_not_silently_overridden():
    """An explicit path is used as given, rather than falling back to $PATH."""
    assert find_binary("inference", binary_path="/nonexistent-scicone-bin") is None
    with pytest.raises(RuntimeError, match="Could not find the SCICoNE executables"):
        SCICoNE(binary_path="/nonexistent-scicone-bin")


def _three_groups():
    """120 cells in three groups whose counts differ two- and three-fold."""
    rng = np.random.RandomState(0)
    centers = rng.poisson(50, size=(3, 20)) * np.array([[1], [2], [3]])
    return np.vstack([rng.poisson(c, size=(40, 20)) for c in centers]).astype(float)


def test_clustering_recovers_groups_reproducibly():
    """learn_tree clusters through PhenoGraph; that step shipped broken once (#135)."""
    sci = SCICoNE()
    data = _three_groups()
    _, sizes, first, _ = sci.condense_segmented_clusters(data, seed=7)
    _, _, second, _ = sci.condense_segmented_clusters(data, seed=7)
    assert sorted(sizes.astype(int)) == [40, 40, 40]
    assert np.array_equal(first, second)


def test_unknown_clustering_algo_is_rejected():
    with pytest.raises(ValueError, match="clustering_algo"):
        SCICoNE().condense_segmented_clusters(_three_groups(), clustering_algo="kmeans")
