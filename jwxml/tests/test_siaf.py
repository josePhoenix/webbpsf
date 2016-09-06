import os.path
import pytest
from ..constants import DATA_ROOT
from ..siaf import SIAF

# The tests are written against a particular version of the SIAF. As long as the format
# and semantics don't change, it's easier to refer back to the earlier SIAF to run
# the tests:
TEST_PRD_VERSION = 'PRDDEVSOC-D-012'  # updated 2016-04-13
TEST_DATA_ROOT = os.path.join(DATA_ROOT, TEST_PRD_VERSION)

def pair_almost_equal(left, right, atol=1e-6):
    left_1, left_2 = left
    right_1, right_2 = right
    assert abs(left_1 - right_1) < atol, "first elements differ by more than {}".format(atol)
    assert abs(left_2 - right_2) < atol, "second elements differ by more than {}".format(atol)
    return True

@pytest.fixture
def nircam_siaf():
    return SIAF(instr='NIRCam', basepath=TEST_DATA_ROOT)

def test_transform_in_to_out(nircam_siaf):
    startx = 1023
    starty = 1024
    nca = nircam_siaf['NIRCAM A']

    assert pair_almost_equal(nca.Det2Sci(startx, starty), (1020., 1020.))
    assert pair_almost_equal(nca.Det2Idl(startx, starty), (0.0, 0.0))
    assert pair_almost_equal(nca.Det2Tel(startx, starty), (87.50, -497.10))

def test_transform_out_to_in(nircam_siaf):
    startV2 = 87.50
    startV3 = -497.10
    nca = nircam_siaf['NIRCAM A']

    assert pair_almost_equal(nca.Sci2Det(1020., 1020), (1023., 1024.))
    assert pair_almost_equal(nca.Tel2Idl(startV2, startV3), (0.0, 0.0))
    assert pair_almost_equal(nca.Tel2Sci(startV2, startV3), (1020., 1020.))
    assert pair_almost_equal(nca.Tel2Det(startV2, startV3), (1023., 1024.))

def test_inverse_transforms(nircam_siaf):
    nca = nircam_siaf['NIRCAM A']

    assert pair_almost_equal(nca.Det2Sci(*nca.Sci2Det(1020., 1020)), (1020., 1020))
    assert pair_almost_equal(nca.Sci2Det(*nca.Det2Sci(1020., 1020)), (1020., 1020))

    assert pair_almost_equal(nca.Tel2Idl(*nca.Idl2Tel(10., 10)), (10., 10))
    assert pair_almost_equal(nca.Idl2Tel(*nca.Tel2Idl(10., 10)), (10., 10))

    assert pair_almost_equal(nca.Tel2Sci(*nca.Sci2Tel(10., 10)), (10., 10))
    assert pair_almost_equal(nca.Sci2Tel(*nca.Tel2Sci(10., 10)), (10., 10))
