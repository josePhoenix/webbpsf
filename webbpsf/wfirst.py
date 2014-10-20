import os.path
import poppy

import webbpsf_core

class WFIRSTInstrument(webbpsf_core.SpaceTelescopeInstrument):
    """
    WFIRSTInstrument contains data and functionality common to WFIRST
    instruments, such as setting the pupil shape
    """
    def __init__(self, *args, **kwargs):
        super(WFIRSTInstrument, self).__init__(*args, **kwargs)
        # TODO:jlong: partly duplicates some logic in SpaceTelescopeInstrument
        self._WebbPSF_basepath = os.path.join(self._WebbPSF_basepath, 'WFIRST')
        self._datapath = os.path.join(self._WebbPSF_basepath, self.name)

        #TODO:jlong: Upon receipt of a real pupil for WFIRST, replace with
        #            a FITS pupil image.
        primary_radius = 2.4  # meters
        secondary_radius = 0.3 * primary_radius
        self.pupil = poppy.CompoundAnalyticalOptic((
            poppy.CircularAperture(radius=2.4), # meters
            poppy.CentralObscuration(secondary_radius=secondary_radius)
        ), name='WFIRST Pupil')
        self.pupilopd = None  # until we have some OPD maps and a FITS pupil of the right shape

class WFIRSTImager(webbpsf_core.WFIRSTInstrument):
    """
    WFIRSTImager represents to the to-be-named wide field imager
    for the WFIRST mission
    """
    def __init__(self):
        scale = 110e-3  # arcsec/px, WFIRST-AFTA SDT report v2 (p. 58)
        super(WFIRSTImager, self).__init__("WFIRSTImager", pixelscale=scale)
        self.pupil = os.path.join(self._WebbPSF_basepath, 'WFIRST')

    def _load_filters(self):
        pass
