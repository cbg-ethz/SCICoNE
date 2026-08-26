import os

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
