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
    """Compare whether two 2-tuples (pairs) of numeric values have an
    absolute difference of < `atol`"""
    left_1, left_2 = left
    right_1, right_2 = right
    assert abs(left_1 - right_1) < atol, "first elements differ by more than {}".format(atol)
    assert abs(left_2 - right_2) < atol, "second elements differ by more than {}".format(atol)
    return True

@pytest.fixture
def nircam_siaf():
    return SIAF(instr='NIRCam', basepath=TEST_DATA_ROOT)

def test_transform_in_to_out(nircam_siaf):
    nca = nircam_siaf['NRCA4_FULL']
    x_det_ref, y_det_ref = nca.XDetRef, nca.YDetRef

    # Test fixed points
    # -- det and sci ref points are == on NRC
    assert pair_almost_equal(nca.Det2Sci(x_det_ref, y_det_ref), (x_det_ref, y_det_ref))
    # -- {X,Y}DetRef is the Idl frame 0, 0 point
    assert pair_almost_equal(nca.Det2Idl(x_det_ref, y_det_ref), (0.0, 0.0))
    # -- The Tel frame ref point is {V2,V3}Ref
    assert pair_almost_equal(nca.Det2Tel(x_det_ref, y_det_ref), (nca.V2Ref, nca.V3Ref))

def test_transform_out_to_in(nircam_siaf):
    nca = nircam_siaf['NRCA4_FULL']
    v2_ref, v3_ref = nca.V2Ref, nca.V3Ref

    # Test fixed points
    # -- {X,Y}DetRef is the Idl frame 0, 0 point
    assert pair_almost_equal(nca.Tel2Idl(v2_ref, v3_ref), (0.0, 0.0))
    # -- det and sci ref points are == on NRC
    assert pair_almost_equal(nca.Tel2Sci(v2_ref, v3_ref), (nca.XDetRef, nca.YDetRef))
    # -- The Tel frame ref point is {V2,V3}Ref
    assert pair_almost_equal(nca.Tel2Det(v2_ref, v3_ref), (nca.XDetRef, nca.YDetRef))

def test_inverse_transforms(nircam_siaf):
    nca = nircam_siaf['NRCA4_FULL']

    assert pair_almost_equal(nca.Det2Sci(*nca.Sci2Det(1020., 1020)), (1020., 1020))
    assert pair_almost_equal(nca.Sci2Det(*nca.Det2Sci(1020., 1020)), (1020., 1020))

    assert pair_almost_equal(nca.Tel2Idl(*nca.Idl2Tel(10., 10)), (10., 10))
    assert pair_almost_equal(nca.Idl2Tel(*nca.Tel2Idl(10., 10)), (10., 10))

    assert pair_almost_equal(nca.Tel2Sci(*nca.Sci2Tel(10., 10)), (10., 10))
    assert pair_almost_equal(nca.Sci2Tel(*nca.Tel2Sci(10., 10)), (10., 10))
