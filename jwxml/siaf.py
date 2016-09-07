import numpy as np
import matplotlib.pyplot as plt

from .constants import PRD_DATA_ROOT
from .utils import iterchildren, etree

import logging
import os
import os.path

_log = logging.getLogger('jwxml')

try:
    from webbpsf import utils
except ImportError:
    utils = None

FRAMES = ('Det', 'Sci', 'Idl', 'Tel')

if not os.path.isdir(PRD_DATA_ROOT):
    raise RuntimeError("Could not find SIAF data "
                       "in {}".format(PRD_DATA_ROOT))

# ---------------------------------------------------------------------------------
#  SIAF related classes

class Aperture(object):
    """ An Aperture, as parsed from the XML.
    All XML nodes are converted into object attributes.

    See JWST-STScI-001550 for the reference on which this implementation was based.

    4 Coordinate systems:
        * Detector:  pixels, in raw detector read out axes orientation ("Det")
        * Science:   pixels, in conventional DMS axes orientation ("Sci")
        * Ideal:     arcsecs relative to aperture reference location. ("Idl")
        * Telescope: arcsecs V2,V3 ("Tel")


    Example
    ========

    ap = some_siaf['desired_aperture_name']     # extract one aperture from a SIAF

    ap.Det2Tel(1024, 512)                       # convert pixel coordinates to sky Tel coords.
                                                # takes pixel coords, returns arcsec
    ap.Idl2Sci( 10, 3)                          # convert Idl coords to Sci pixels
                                                # takes arcsec, returns pixels

    # there exist functions for all of the possible
    # {Tel,Idl,Sci,Det}2{Tel,Idl,Sci,Det} combinations.

    # you can also specify frames by string:
    ap.convert(1024, 512, frame_from='Det', frame_to='Tel')  # same as first example above

    ap.corners('Tel')                           # Get detector corners in Tel frame
    ap.center('Tel')                            # Get the reference point defined in the SIAF
                                                # this is typically the center of this region

    ap.plot('Idl', annotate=True, title=True)   # plot coordinates in Idl frame
    ap.plotDetectorChannels()                   # color in the readout channels
    """

    def __init__(self, xmlnode, instrument=None):

        self.instrument = instrument
        convfactors = {'RADIANS': 1, 'DEGREES': np.pi / 180, 'ARCSECS': np.pi / 180 / 60 / 60}

        for node in iterchildren(xmlnode):
            tag = node.tag.replace('{http://www.stsci.edu/SIAF}', '')
            if len(node.getchildren()) == 0:
                # if doens't have children,
                try:
                    value = float(node.text)  # do we care about ints vs floats?
                except (ValueError, TypeError):
                    value = node.text
                self.__dict__[tag] = value
            else:
                # if does have children:
                if '{http://www.stsci.edu/SIAF}units' in [c.tag for c in node.getchildren()]:
                    # this will be an angle/units pair. units are either in arcsec or degrees.
                    # Convert to radians in any case for internal use.
                    unit = node.find('{http://www.stsci.edu/SIAF}units').text
                    value = float(node.find('{http://www.stsci.edu/SIAF}value').text)
                    value *= convfactors[unit]
                    self.__dict__[tag] = value
                elif '{http://www.stsci.edu/SIAF}elt' in [c.tag for c in node.getchildren()]:
                    #  an array of values which should go to an NDarray
                    elts = [float(c.text) for c in
                            iterchildren(node, '{http://www.stsci.edu/SIAF}elt')]
                    self.__dict__[tag] = np.asarray(elts)

                else:
                    raise NotImplemented("Not sure how to parse that node.")

        # pack things into NDarrays for convenient access
        # first the vertices
        self.XIdlVert = np.asarray((self.XIdlVert1, self.XIdlVert2, self.XIdlVert3, self.XIdlVert4))
        self.YIdlVert = np.asarray((self.YIdlVert1, self.YIdlVert2, self.YIdlVert3, self.YIdlVert4))

        # then the transformation coefficients
        if self.Sci2IdlDeg is not None:
            self.Sci2IdlDeg = int(self.Sci2IdlDeg)
            self.Sci2IdlCoeffs_X = np.zeros((self.Sci2IdlDeg + 1, self.Sci2IdlDeg + 1))
            self.Sci2IdlCoeffs_Y = np.zeros((self.Sci2IdlDeg + 1, self.Sci2IdlDeg + 1))
            self.Idl2SciCoeffs_X = np.zeros((self.Sci2IdlDeg + 1, self.Sci2IdlDeg + 1))
            self.Idl2SciCoeffs_Y = np.zeros((self.Sci2IdlDeg + 1, self.Sci2IdlDeg + 1))
            for i in range(1, self.Sci2IdlDeg + 1):
                for j in range(0, i + 1):
                    self.Sci2IdlCoeffs_X[i, j] = self.__dict__['Sci2IdlX{0:1d}{1:1d}'.format(i, j)]
                    self.Sci2IdlCoeffs_Y[i, j] = self.__dict__['Sci2IdlY{0:1d}{1:1d}'.format(i, j)]
                    self.Idl2SciCoeffs_X[i, j] = self.__dict__['Idl2SciX{0:1d}{1:1d}'.format(i, j)]
                    self.Idl2SciCoeffs_Y[i, j] = self.__dict__['Idl2SciY{0:1d}{1:1d}'.format(i, j)]

    def __repr__(self):
        return "<jwxml.Aperture object AperName={0} >".format(self.AperName)

    # --- the actual fundamental transformation code follows in these next routines:
    def Det2Sci(self, XDet, YDet):
        """ Detector to Science, following Section 4.1 of JWST-STScI-001550"""
        XDet = np.asarray(XDet, dtype=float)
        YDet = np.asarray(YDet, dtype=float)
        ang = np.deg2rad(self.DetSciYAngle)
        XSci = (
            self.XSciRef +
            self.DetSciParity * (
                (XDet - self.XDetRef) * np.cos(ang) +
                (YDet - self.YDetRef) * np.sin(ang)
            )
        )
        YSci = (
            self.YSciRef -
            (XDet - self.XDetRef) * np.sin(ang) +
            (YDet - self.YDetRef) * np.cos(ang)
        )
        return XSci, YSci

    def Sci2Det(self, XSci, YSci):
        """ Science to Detector, following Section 4.1 of JWST-STScI-001550"""
        XSci = np.asarray(XSci, dtype=float)
        YSci = np.asarray(YSci, dtype=float)

        ang = np.deg2rad(self.DetSciYAngle)
        XDet = (
            self.XDetRef +
            self.DetSciParity * (XSci - self.XSciRef) * np.cos(ang) -
            (YSci - self.YSciRef) * np.sin(ang)
        )
        YDet = (
            self.YDetRef +
            self.DetSciParity * (XSci - self.XSciRef) * np.sin(ang) +
            (YSci - self.YSciRef) * np.cos(ang)
        )
        return XDet, YDet

    def Sci2Idl(self, XSci, YSci):
        """ Convert Sci to Idl
        input in pixel, output in arcsec """
        dX = np.asarray(XSci, dtype=float) - self.XSciRef
        dY = np.asarray(YSci, dtype=float) - self.YSciRef

        degree = self.Sci2IdlDeg
        XIdl = np.zeros_like(np.asarray(XSci), dtype=float)
        YIdl = np.zeros_like(np.asarray(YSci), dtype=float)

        for i in range(1, degree + 1):
            for j in range(0, i + 1):
                XIdl += self.Sci2IdlCoeffs_X[i, j] * dX ** (i - j) * dY ** j
                YIdl += self.Sci2IdlCoeffs_Y[i, j] * dX ** (i - j) * dY ** j

        return XIdl, YIdl

    def Idl2Sci(self, XIdl, YIdl):
        """ Convert Idl to  Sci
        input in arcsec, output in pixels """
        XIdl = np.asarray(XIdl, dtype=float)
        YIdl = np.asarray(YIdl, dtype=float)

        degree = self.Sci2IdlDeg
        XSci = np.zeros_like(np.asarray(XIdl), dtype=float)
        YSci = np.zeros_like(np.asarray(YIdl), dtype=float)

        for i in range(1, degree + 1):
            for j in range(0, i + 1):
                XSci += self.Idl2SciCoeffs_X[i, j] * XIdl ** (i - j) * YIdl ** j
                YSci += self.Idl2SciCoeffs_Y[i, j] * XIdl ** (i - j) * YIdl ** j

        return XSci + self.XSciRef, YSci + self.YSciRef

    def Idl2Tel(self, XIdl, YIdl):
        """ Convert Idl to  Tel

        input in arcsec, output in arcsec

        WARNING
        --------
        This is an implementation of the planar approximation, which is adequate for most
        purposes but may not be for all. Error is about 1.7 mas at 10 arcminutes from the tangent
        point. See JWST-STScI-1550 for more details.
        """
        XIdl = np.asarray(XIdl, dtype=float)
        YIdl = np.asarray(YIdl, dtype=float)

        # V2Ref and V3Ref are now in arcseconds in the XML file
        ang = np.deg2rad(self.V3IdlYAngle)
        V2 = self.V2Ref + self.VIdlParity * XIdl * np.cos(ang) + YIdl * np.sin(ang)
        V3 = self.V3Ref - self.VIdlParity * XIdl * np.sin(ang) + YIdl * np.cos(ang)
        return V2, V3

    def Tel2Idl(self, V2, V3):
        """ Convert Tel to Idl

        input in arcsec, output in arcsec

        This transformation involves going from global V2,V3 to local angles with respect to some
        reference point, and possibly rotating the axes and/or flipping the parity of the X axis.


        WARNING
        --------
        This is an implementation of the planar approximation, which is adequate for most
        purposes but may not be for all. Error is about 1.7 mas at 10 arcminutes from the tangent
        point. See JWST-STScI-1550 for more details.
        """

        dV2 = np.asarray(V2, dtype=float) - self.V2Ref
        dV3 = np.asarray(V3, dtype=float) - self.V3Ref
        ang = np.deg2rad(self.V3IdlYAngle)

        XIdl = self.VIdlParity * (dV2 * np.cos(ang) - dV3 * np.sin(ang))
        YIdl = dV2 * np.sin(ang) + dV3 * np.cos(ang)
        return XIdl, YIdl

    # --- and now some compound transformations that are less fundamental. This just nests calls
    # to the above.

    def Det2Idl(self, *args):
        return self.Sci2Idl(*self.Det2Sci(*args))

    def Det2Tel(self, *args):
        return self.Idl2Tel(*self.Sci2Idl(*self.Det2Sci(*args)))

    def Sci2Tel(self, *args):
        return self.Idl2Tel(*self.Sci2Idl(*args))

    def Idl2Det(self, *args):
        return self.Sci2Det(*self.Idl2Sci(*args))

    def Tel2Sci(self, *args):
        return self.Idl2Sci(*self.Tel2Idl(*args))

    def Tel2Det(self, *args):
        return self.Sci2Det(*self.Idl2Sci(*self.Tel2Idl(*args)))

    # --- now, functions other than direct coordinate transformations
    def convert(self, X, Y, frame_from, frame_to):
        """ Generic conversion routine, that calls one of the
        specific conversion routines based on the provided frame names as strings. """

        if frame_from not in FRAMES or frame_to not in FRAMES:
            raise ValueError("You must specify a frame_from value: {}".format(', '.join(FRAMES)))

        if frame_from == frame_to:
            return X, Y  # null transformation

        # With valid frame_from and frame_to, this method must exist:
        conversion_method = getattr(self, '{}2{}'.format(frame_from, frame_to))

        return conversion_method(X, Y)

    def corners(self, frame='Idl'):
        """ Return coordinates of the aperture outline"""
        return self.convert(self.XIdlVert, self.YIdlVert, 'Idl', frame)

    def center(self, frame='Tel'):
        """ Return the defining center point of the aperture"""
        return self.convert(self.V2Ref, self.V3Ref, 'Tel', frame)

    def plot(self, frame='Idl', label=True, ax=None, title=True, units='arcsec', annotate=False,
             color=None):
        """ Plot this one aperture

        Parameters
        -----------
        frame : str
            Which coordinate system to plot in: 'Tel', 'Idl', 'Sci', 'Det'
        label : bool
            Add text label stating aperture name
        ax : matplotlib.Axes
            Desired destination axes to plot into (If None, current
            axes are inferred from pyplot.)
        units : str
            one of 'arcsec', 'arcmin', 'deg'
        annotate : bool
            Add annotations for detector (0,0) pixels
        title : str
            If set, add a label to the plot indicating which frame was plotted.
        color : matplotlib-compatible color
            Color specification for this aperture's outline,
            passed through to `matplotlib.Axes.plot`
        """
        if units is None:
            units = 'arcsec'

        if ax is None:
            ax = plt.gca()
            ax.set_aspect('equal')
            if frame == 'Idl' or frame == 'Tel':
                ax.set_xlabel('V2 [{0}]'.format(units))
                ax.set_ylabel('V3 [{0}]'.format(units))

            elif frame == 'Sci' or frame == 'Det':
                ax.set_xlabel('X pixels [{0}]'.format(frame))
                ax.set_ylabel('Y pixels [{0}]'.format(frame))

        x, y = self.corners(frame=frame)

        if units.lower() == 'arcsec':
            scale = 1
        elif units.lower() == 'arcmin':
            scale = 1. / 60
        elif units.lower() == 'deg':
            scale = 1. / 60 / 60
        else:
            raise ValueError("Unknown units: " + units)

        x2 = np.concatenate([x, [x[0]]])  # close the box
        y2 = np.concatenate([y, [y[0]]])

        # convert arcsec to arcmin and plot
        if color is not None:
            ax.plot(x2 * scale, y2 * scale, color=color)
        else:
            ax.plot(x2 * scale, y2 * scale)

        if label:
            # partially mitigate overlapping NIRCam labels
            rotation = 30 if self.AperName.startswith('NRC') else 0
            ax.text(
                x.mean() * scale, y.mean() * scale, self.AperName,
                verticalalignment='center',
                horizontalalignment='center',
                rotation=rotation,
                color=ax.lines[-1].get_color()
            )
        if title:
            ax.set_title("{0} frame".format(frame))
        if annotate:
            self.plotDetectorOrigin(frame=frame)

    def plotDetectorOrigin(self, frame='Idl', which='both'):
        """ Draw red and blue squares to indicate the raw detector
        readout and science frame readout, respectively

        Parameters
        -----------
        which : str
            Which detector origin to plot: 'both', 'Det', 'Sci'
        frame : str
            Which coordinate system to plot in: 'Tel', 'Idl', 'Sci', 'Det'
        """

        # raw detector frame
        if which.lower() == 'det' or which.lower() == 'both':
            c1, c2 = self.convert(0, 0, 'Det', frame)
            plt.plot(c1, c2, color='red', marker='s', markersize=9)

        # science frame
        if which.lower() == 'sci' or which.lower() == 'both':
            c1, c2 = self.convert(0, 0, 'Sci', frame)
            plt.plot(c1, c2, color='blue', marker='s')

    def plotDetectorChannels(self, frame='Idl', color='0.5', alpha=0.3, evenoddratio=0.5):
        """ Mark on the plot the various detector readout channels

        These are depicted as alternating light/dark bars to show the
        regions read out by each of the output amps.

        Parameters
        ----------
        frame : str
            Which coordinate system to plot in: 'Tel', 'Idl', 'Sci', 'Det'
        color : matplotlib-compatible color
            Color specification for the amplifier shaded region,
            passed through to `matplotlib.patches.Polygon` as `facecolor`
        alpha : float
            Opacity of odd-numbered amplifier region overlays
            (for even, see `evenoddratio`)
        evenoddratio : float
            Ratio of opacity between even and odd amplifier region
            overlays
       """

        import matplotlib
        if self.instrument == 'MIRI':
            npixels = 1024
        else:
            npixels = 2048
        ch = npixels / 4

        ax = plt.gca()
        pts = ((0, 0), (ch, 0), (ch, npixels), (0, npixels))
        for chan in range(4):
            plotpoints = np.zeros((4, 2))
            for i, xy in enumerate(pts):
                plotpoints[i] = self.convert(xy[0] + chan * ch, xy[1], 'Det', frame)
            chan_alpha = alpha if chan % 2 == 1 else alpha * evenoddratio
            rect = matplotlib.patches.Polygon(
                plotpoints,
                closed=True,
                alpha=chan_alpha,
                facecolor=color,
                edgecolor='none',
                lw=0
            )
            ax.add_patch(rect)


class SIAF(object):
    """ Science Instrument Aperture File

    This is a class interface to SIAF information stored in an XML file.
    It lets you read (only) the SIAF information, retrieve apertures,
    plot them, and transform coordinates accordingly.

    This class is basically just a container. See the Aperture class for
    the detailed implementation of the transformations.

    Briefly, this class acts like a dict containing Aperture objects, accessible
    using their names defined in the SIAF

    Examples
    ---------

    fgs_siaf = SIAF('FGS')
    fgs_siaf.apernames                # returns a list of aperture names
    ap = fgs_siaf['FGS1_FULL']   # returns an aperture object
    ap.plot(frame='Tel')              # plot one aperture
    fgs_siaf.plot()                   # plot all apertures in this file

    """

    def __init__(self, instr='NIRISS', filename=None, basepath=None):
        """ Read a SIAF from disk

        Parameters
        -----------
        instr : string
            one of 'NIRCam', 'NIRSpec', 'NIRISS', 'MIRI', 'FGS'; case sensitive.
        basepath : string
            Directory to look in for SIAF files
        filename : string, optional
            Alternative method to specify a specific SIAF XML file.
        """

        if instr not in ['NIRCam', 'NIRSpec', 'NIRISS', 'MIRI', 'FGS']:
            raise ValueError(
                "Invalid instrument name: {0}. Note that this is case sensitive.".format(instr))

        self.instrument = instr

        if filename is None:
            if basepath is None:
                basepath = PRD_DATA_ROOT
            self.filename = os.path.join(basepath, instr + '_SIAF.xml')
        else:
            self.filename = filename

        self.apertures = {}

        self._tree = etree.parse(self.filename)

        self._last_plot_frame = None

        for entry in self._tree.getroot().iter('SiafEntry'):
            aperture = Aperture(entry, instrument=self.instrument)
            self.apertures[aperture.AperName] = aperture

    def __getitem__(self, key):
        return self.apertures[key]

    def __len__(self):
        return len(self.apertures)

    @property
    def apernames(self):
        """ List of aperture names defined in this SIAF"""
        return self.apertures.keys()

    def _getFullApertures(self):
        """ Return whichever subset of apertures correspond to the entire detectors. This is a
        helper function for the various plotting routines following"""
        fullaps = []
        if self.instrument == 'NIRCam':
            fullaps.append(self.apertures['NRCA5_FULL'])
            fullaps.append(self.apertures['NRCB5_FULL'])
        elif self.instrument == 'NIRSpec':
            fullaps.append(self.apertures['NRS_FULL_MSA1'])
            fullaps.append(self.apertures['NRS_FULL_MSA2'])
            fullaps.append(self.apertures['NRS_FULL_MSA3'])
            fullaps.append(self.apertures['NRS_FULL_MSA4'])
        elif self.instrument == 'NIRISS':
            fullaps.append(self.apertures['NIS-CEN'])
        elif self.instrument == 'MIRI':
            fullaps.append(self.apertures['MIRIM_FULL'])
        elif self.instrument == 'FGS':
            fullaps.append(self.apertures['FGS1_FULL'])
            fullaps.append(self.apertures['FGS2_FULL'])
        return fullaps

    def plot(self, frame='Tel', names=None, label=True, units=None, clear=True, annotate=False,
             subarrays=True):
        """ Plot all apertures in this SIAF

        Parameters
        -----------
        names : list of strings
            A subset of aperture names, if you wish to plot only a subset
        subarrays : bool
            Plot all the minor subarrays if True, else just plot the "main" apertures
        label : bool
            Add text labels stating aperture names
        units : str
            one of 'arcsec', 'arcmin', 'deg'
        clear : bool
            Clear plot before plotting (set to false to overplot)
        annotate : bool
            Add annotations for detector (0,0) pixels
        frame : str
            Which coordinate system to plot in: 'Tel', 'Idl', 'Sci', 'Det'
        """
        if clear: plt.clf()
        ax = plt.subplot(111)
        ax.set_aspect('equal')

        # which list of apertures to iterate over?
        if subarrays:
            iterable = self.apertures.values
        else:
            iterable = self._getFullApertures

        for ap in iterable():
            if names is not None:
                if ap.AperName not in names: continue

            ap.plot(frame=frame, label=label, ax=ax, units=None)
            if annotate:
                ap.plotDetectorOrigin(frame=frame)
        ax.set_xlabel('V2 [arcsec]')
        ax.set_ylabel('V3 [arcsec]')

        if frame == 'Tel' or frame == 'Idl':
            # enforce V2 increasing toward the left
            ax.autoscale_view(True, True, True)
            xlim = ax.get_xlim()
            if xlim[1] > xlim[0]: ax.set_xlim(xlim[::-1])
            ax.set_autoscalex_on(True)

        self._last_plot_frame = frame

    def plotDetectorOrigin(self, which='both', frame=None):
        """ Mark on the plot the detector's origin in Det and Sci coordinates

        Parameters
        -----------
        which : str
            Which detector origin to plot: 'both', 'Det', 'Sci'
        frame : str
            Which coordinate system to plot in: 'Tel', 'Idl', 'Sci', 'Det'
            Optional if you have already called plot() to specify a
            coordinate frame.

        """
        if frame is None: frame = self._last_plot_frame
        for ap in self._getFullApertures():
            ap.plotDetectorOrigin(frame=frame, which=which)

    def plotDetectorChannels(self, frame=None):
        """ Mark on the plot the various detector readout channels

        These are depicted as alternating light/dark bars to show the
        regions read out by each of the output amps.

        Parameters
        ----------
        frame : str
            Which coordinate system to plot in: 'Tel', 'Idl', 'Sci', 'Det'
            Optional if you have already called plot() to specify a
            coordinate frame.

        """
        if frame is None: frame = self._last_plot_frame
        for ap in self._getFullApertures():
            ap.plotDetectorChannels(frame=frame)


def plotAllSIAFs(subarrays=True, showorigin=True, showchannels=True, **kwargs):
    """ Plot All instrument """

    for instr in ['NIRCam', 'NIRISS', 'NIRSpec', 'FGS', 'MIRI']:
        aps = SIAF(instr, **kwargs)
        print("{0} has {1} apertures".format(aps.instrument, len(aps)))

        aps.plot(clear=False, subarrays=subarrays, **kwargs)
        if showorigin: aps.plotDetectorOrigin()
        if showchannels: aps.plotDetectorChannels()


def plotMainSIAFs(label=False, darkbg=False, **kwargs):
    if darkbg:
        col_imaging = 'aqua'
        col_coron = 'lime'
        col_msa = 'violet'
    else:
        col_imaging = 'blue'
        col_coron = 'green'
        col_msa = 'magenta'

    nircam = SIAF('NIRCam')
    niriss = SIAF('NIRISS')
    fgs = SIAF('FGS')
    nirspec = SIAF('NIRSpec')
    miri = SIAF('MIRI')

    im_aps = [
        nircam['NRCA5_FULL'],
        nircam['NRCB5_FULL'],
        niriss['NIS-CEN'],
        miri['MIRIM_ILLUM'],
        fgs['FGS1_FULL'],
        fgs['FGS2_FULL']
    ]

    for letter in ['A', 'B']:
        for num in range(5):
            im_aps.append(nircam['NRC{}{}_FULL'.format(letter, num + 1)])

    coron_aps = [
        nircam['NRCA2_MASK210R'],
        nircam['NRCA4_MASKSWB'],
        nircam['NRCA5_MASK335R'],
        nircam['NRCA5_MASK430R'],
        nircam['NRCA5_MASKLWB'],
        nircam['NRCB3_MASKSWB'],
        nircam['NRCB1_MASK210R'],
        nircam['NRCB5_MASK335R'],
        nircam['NRCB5_MASK430R'],
        nircam['NRCB5_MASKLWB'],
        miri['MIRIM_MASK1065'],
        miri['MIRIM_MASK1140'],
        miri['MIRIM_MASK1550'],
        miri['MIRIM_MASKLYOT']
    ]
    msa_aps = [nirspec['NRS_FULL_MSA' + str(n + 1)] for n in range(4)]
    msa_aps.append(nirspec['NRS_S1600A1_SLIT'])  # square aperture

    for aplist, col in zip([im_aps, coron_aps, msa_aps], [col_imaging, col_coron, col_msa]):
        for ap in aplist:
            ap.plot(color=col, frame='Tel', label=label, **kwargs)

    # ensure V2 increases to the left
    ax = plt.gca()
    xlim = ax.get_xlim()
    if xlim[0] < xlim[1]:
        ax.set_xlim(xlim[::-1])

