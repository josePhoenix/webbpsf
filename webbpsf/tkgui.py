#!/usr/bin/env python
"""
webbpsf.tkgui -- Tk/TTk graphical interface wrapping WebbPSF
"""

import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import astropy.io.fits as fits
import Tkinter as tk
import ttk
import tkMessageBox
import tkFileDialog
import logging
_log = logging.getLogger('webbpsf')

try:
    import pysynphot
    _HAS_PYSYNPHOT_INSTALLED = True
    if os.getenv('PYSYN_CDBS') is not None:
        _HAS_PYSYNPHOT_DATA = True
    else:
        _HAS_PYSYNPHOT_DATA = False
except ImportError:
    pysynphot = None
    _HAS_PYSYNPHOT_INSTALLED = False
    _HAS_PYSYNPHOT_DATA = False


import poppy
from webbpsf import webbpsf_core



class ConfigurationError(Exception):
    pass

class Options(object):
    """
    Options for GUI dropdowns

    >>> calculation_options = Options(firstopt='First Option', \
    >>>     secondopt='Second Option')
    >>> calculation_options.firstopt
    'First Option'
    >>> calculation_options.names
    ('First Option', 'Second Option')
    >>> calculation_options.short_names
    ('firstopt', 'secondopt')
    >>> calculation_options.validate('First Option')
    True
    >>> calculation_options.validate('Not an Option')
    False
    """
    def __init__(self, **kwargs):
        """
        Set attributes on the Options instance according to
        kwargname=value for each keyword argument.
        """
        self.ordered_options = []
        for key, value in kwargs.items():
            if key in ['names', 'short_names', 'validate']:
                # prevent accidental name collisions
                raise ConfigurationError("The name {0} is in use by the API")
            setattr(self, key, value)
            self.ordered_options.append((key, value))

    @property
    def names(self):
        """
        Names defined for this set of options
        """
        return [long_name for _, long_name in self.ordered_options]

    @property
    def short_names(self):
        """
        Short (attribute) names defined for this set of options
        """
        return [short_name for short_name, _ in self.ordered_options]

    def validate(self, value):
        """
        Validates that `value` is a valid choice from these options
        """
        return value in self.names

    def __len__(self):
        return len(self.ordered_options)

    def __getitem__(self, item):
        return self.ordered_options[item]

    def __iter__(self):
        return self.names.__iter__()


# Control-building helpers

def _add_labeled_dropdown(controller, name, root, label="Entry:", label_pad=(1, 1), values=None,
                          width=5, position=(0, 0), default=None, **kwargs):
    """Convenient wrapper for adding a labeled Combobox"""

    ttk.Label(root, text=label).grid(
        row=position[0],
        column=position[1],
        sticky='W',
        padx=label_pad
    )

    controller.vars[name] = tk.StringVar()
    controller.widgets[name] = ttk.Combobox(
        root,
        textvariable=controller.vars[name],
        width=width,
        state='readonly'
    )
    controller.widgets[name].grid(row=position[0], column=position[1] + 1, **kwargs)
    controller.widgets[name]['values'] = values

    if default is None:
        default = values[0]
    controller.widgets[name].set(default)


def _add_labeled_entry(controller, name, root, label="Entry:", label_pad=(1, 1), value="",
                       width=5, position=(0, 0), postlabel=None, **kwargs):
    """Convenient wrapper for adding a labeled Entry"""
    ttk.Label(root, text=label).grid(
        row=position[0],
        column=position[1],
        sticky='W',
        padx=label_pad
    )

    controller.vars[name] = tk.StringVar()
    controller.widgets[name] = ttk.Entry(
        root,
        textvariable=controller.vars[name],
        width=width
    )
    controller.widgets[name].insert(0, value)
    controller.widgets[name].grid(row=position[0], column=position[1]+1, **kwargs)

    if postlabel is not None:
        ttk.Label(root, text=postlabel).grid(
            row=position[0],
            column=position[1]+2,
            sticky='W'
        )


class PSFGenerationGUI(object):
    """Base Class for a PSF generation GUI using Tkinter & TTK native widgets"""
    OUTPUT_OPTS = Options(
        oversampled='Oversampled image',
        detector_scale='Detector sampled image',
        both_as_fits='Both as FITS extensions'
    )
    PERFECT_OPD = "Zero OPD (perfect)"

    def __init__(self):
        # set basic options to defaults
        self.field_of_view = 5
        self.nlambda = None
        self.detector_oversampling = 2
        self.fft_oversampling = 2
        self.advanced_options = {
            'parity': 'any',
            'force_coron': False,
            'no_sam': False,
            'psf_vmin': 1e-8,
            'psf_vmax': 1.0,
            'psf_scale': 'log',
            'psf_cmap_str': 'Jet (blue to red)',
            'psf_normalize': 'Total',
            'psf_cmap': matplotlib.cm.jet
        }
        self.sptype = None  # TODO:jlong: reasonable default
        self.iname = None  # ditto
        self.filter = None  # ditto
        self.inst = None  # set during init?
        self.opd_name = None
        self.opd_i = None

        # init the object and subobjects
        self.instrument = {}
        self.widgets = {}
        self.vars = {}
        self.psf_hdulist = None
        # init Tk window
        self.root = tk.Tk()
        self.root.geometry('+50+50')

    def mainloop(self):
        self.root.mainloop()

    def _create_widgets(self):
        """Create a nice GUI using the enhanced widget set provided by
        the ttk extension to Tkinter, available in Python 2.7 or newer
        """
        #---- create the GUIs
        self.root.title("James Webb Space Telescope PSF Calculator")

        frame = ttk.Frame(self.root)

        #-- source properties
        source_labelframe = ttk.LabelFrame(frame, text='Source Properties')
        self._populate_source_properties(source_labelframe)
        source_labelframe.columnconfigure(2, weight=1)
        source_labelframe.grid(row=1, sticky='E,W', padx=10, pady=5)

        #-- instrument config
        lf = ttk.LabelFrame(frame, text='Instrument Config')
        notebook = ttk.Notebook(lf)
        self.widgets['tabset'] = notebook
        notebook.pack(fill='both')

        self._populate_instrument_config(notebook)

        self.ev_update_opd_labels()
        lf.grid(row=2, sticky='E,W', padx=10, pady=5)
        notebook.select(0)

        #-- calculation options

        calc_options_frame = ttk.LabelFrame(frame, text='Calculation Options')
        self._populate_calculation_options(calc_options_frame)
        calc_options_frame.columnconfigure(2, weight=1)
        calc_options_frame.columnconfigure(4, weight=1)
        calc_options_frame.grid(row=5, sticky='E,W', padx=10, pady=15)

        #-- control buttons
        control_buttons_frame = ttk.Frame(frame)

        def addbutton(text, command, pos, disabled=False):
            """Shorthand inner function for adding buttons"""
            self.widgets[text] = ttk.Button(control_buttons_frame, text=text, command=command)
            self.widgets[text].grid(column=pos, row=0, sticky='E')
            if disabled:
                self.widgets[text].state(['disabled'])

        addbutton('Compute PSF', self.ev_calc_psf, 0)
        addbutton('Display PSF', self.ev_display_psf, 1, disabled=True)
        addbutton('Display profiles', self.ev_display_profiles, 2, disabled=True)
        addbutton('Save PSF As...', self.ev_save_as, 3, disabled=True)
        addbutton('More options...', self.ev_options, 4, disabled=False)

        ttk.Button(control_buttons_frame, text='Quit', command=self.quit).grid(column=5, row=0)
        control_buttons_frame.columnconfigure(2, weight=1)
        control_buttons_frame.columnconfigure(4, weight=1)
        control_buttons_frame.grid(row=6, sticky='E,W', padx=10, pady=15)

        frame.grid(row=0, sticky='N,E,S,W')
        frame.columnconfigure(0, weight=1)
        frame.rowconfigure(0, weight=1)
        self.root.columnconfigure(0, weight=1)
        self.root.rowconfigure(0, weight=1)

    def _populate_source_properties(self, source_labelframe):
        if _HAS_PYSYNPHOT_DATA:
            _add_labeled_dropdown(
                self,
                "SpType",
                source_labelframe,
                label='    Spectral Type:',
                values=poppy.specFromSpectralType("", return_list=True),
                default='G0V',
                width=25,
                position=(0, 0),
                sticky='W'
            )
            plot_button = ttk.Button(source_labelframe, text='Plot spectrum',
                                     command=self.ev_plot_spectrum)
            plot_button.grid(row=0, column=2, sticky='E', columnspan=4)

        frame_root = ttk.Frame(source_labelframe)

        _add_labeled_entry(
            self,
            "source_off_r",
            frame_root,
            label='    Source Position: r=',
            value='0.0',
            width=5,
            position=(1, 0),
            sticky='W'
        )
        _add_labeled_entry(
            self,
            "source_off_theta",
            frame_root,
            label='arcsec,  PA=',
            value='0',
            width=3,
            position=(1, 2),
            sticky='W'
        )

        self.vars["source_off_centerpos"] = tk.StringVar()
        self.vars["source_off_centerpos"].set('corner')

        ttk.Label(frame_root, text='deg, centered on ').grid(row=1, column=4)
        pixel = ttk.Radiobutton(frame_root, text='pixel',
                                variable=self.vars["source_off_centerpos"], value='pixel')
        pixel.grid(row=1, column=5)
        corner = ttk.Radiobutton(frame_root, text='corner',
                                 variable=self.vars["source_off_centerpos"], value='corner')
        corner.grid(row=1, column=6)
        frame_root.grid(row=1, column=0, columnspan=5, sticky='W')

    def _populate_calculation_options(self, calc_opts_root):
        row_counter = 0
        _add_labeled_entry(
            self,
            'FOV',
            calc_opts_root,
            label='Field of View:',
            width=3,
            value=str(self.field_of_view),
            postlabel='arcsec/side',
            position=(row_counter, 0),
        )
        row_counter += 1
        _add_labeled_entry(
            self,
            'detector_oversampling',
            calc_opts_root,
            label='Output Oversampling:',
            width=3,
            value=str(self.detector_oversampling),
            postlabel='x finer than instrument pixels       ',
            position=(row_counter, 0)
        )
        row_counter += 1
        _add_labeled_entry(
            self,
            'fft_oversampling',
            calc_opts_root,
            label='Coronagraph FFT Oversampling:',
            width=3,
            value=str(self.fft_oversampling),
            postlabel='x finer than Nyquist',
            position=(row_counter, 0)
        )
        row_counter += 1
        _add_labeled_entry(
            self,
            'nlambda',
            calc_opts_root,
            label='# of wavelengths:',
            width=3,
            value='',
            position=(row_counter, 0),
            postlabel='Leave blank for autoselect'
        )
        row_counter += 1

        _add_labeled_dropdown(
            self,
            "jitter",
            calc_opts_root,
            label='Jitter model:',
            values=['Just use OPDs'],  # 'Just use OPDs', 'Gaussian blur', 'Accurate yet SLOW grid'
            width=20,
            position=(row_counter, 0),
            sticky='W',
            columnspan=2
        )
        row_counter += 1
        _add_labeled_dropdown(
            self,
            "output_format",
            calc_opts_root,
            label='Output Format:',
            values=self.OUTPUT_OPTS.names,
            width=30,
            position=(row_counter, 0),
            sticky='W',
            columnspan=2
        )

        calc_opts_root.grid(row=4, sticky='E,W', padx=10, pady=5)

    def _populate_instrument_config(self, notebook):
        raise NotImplementedError("Subclasses must implement _populate_instrument_config")

    def quit(self):
        """Quit the GUI"""
        confirm_quit = tkMessageBox.askyesno(
            message='Are you sure you want to quit?',
            icon='question',
            title='Confirm quit'
        )
        if confirm_quit:
            self.root.destroy()

    def ev_save_as(self):
        """Event handler for Save As of output PSFs"""
        filename = tkFileDialog.asksaveasfilename(
            initialfile='PSF_%s_%s.fits' % (self.iname, self.filter),
            filetypes=[('FITS', '.fits')],
            defaultextension='.fits',
            parent=self.root
        )
        if len(filename) > 0:
            self.psf_hdulist.writeto(filename)
            print "Saved to %s" % filename

    def ev_options(self):
        """Event handler to open advanced options dialog"""
        dialog = WebbPSFOptionsDialog(
            self.root,
            input_options=self.advanced_options
        )
        self.advanced_options = dialog.results

    def ev_plot_spectrum(self):
        """Event handler for Plot Spectrum button"""
        self._update_from_gui()

        #sptype = self.widgets['SpType'].get()
        #iname = self.widgets[self.widgets['tabset'].select()]
        print "Spectral type is", self.sptype
        print "Selected instrument tab is", self.iname
        #if iname != 'TFI':
            #filter = self.widgets[self.iname+"_filter"].get()
        print "Selected instrument filter is", self.filter

        plt.clf()

        ax1 = plt.subplot(311)
        spectrum = poppy.specFromSpectralType(self.sptype)
        synplot(spectrum)
        ax1.set_ybound(1e-6, 1e8)  # hard coded for now
        ax1.yaxis.set_major_locator(matplotlib.ticker.LogLocator(base=1000))
        legend_font = matplotlib.font_manager.FontProperties(size=10)
        ax1.legend(loc='lower right', prop=legend_font)

        ax2 = plt.subplot(312, sharex=ax1)
        ax2.set_ybound(0, 1.1)
        band = self.inst._getSynphotBandpass(self.inst.filter)
        band.name = "%s %s" % (self.iname, self.inst.filter)
        synplot(band)
        legend_font = matplotlib.font_manager.FontProperties(size=10)
        plt.legend(loc='lower right', prop=legend_font)

        ax3 = plt.subplot(313, sharex=ax1)
        if self.nlambda is None:
            # Automatically determine number of appropriate wavelengths.
            # Make selection based on filter configuration file
            nlambda = self.inst._filter_nlambda_default.get(self.filter, 10)
        else:
            nlambda = self.nlambda
        ax1.set_xbound(0.1, 100)
        plt.draw()
        waves, weights = self.inst._getWeights(spectrum, nlambda=nlambda)

        wave_step = waves[1] - waves[0]
        plot_waves = 1e6 * np.concatenate((
            [waves[0] - wave_step],
            waves,
            [waves[-1] + wave_step]
        ))
        plot_weights = np.concatenate((
            [0],
            weights,
            [0]
        ))

        plt.ylabel("Weight")
        plt.xlabel(r"Wavelength [$\mu$m]")

        ax3.plot(plot_waves, plot_weights, drawstyle='steps-mid')

        ax1.set_xbound(0.1, 100)

        _refresh_window()

    def ev_calc_psf(self):
        """Event handler for PSF Calculations"""
        self._update_from_gui()

        if _HAS_PYSYNPHOT_DATA:
            source = poppy.specFromSpectralType(self.sptype)
        else:
            source = None  # generic flat spectrum

        _refresh_window()

        self.psf_hdulist = self.inst.calcPSF(
            source=source,
            detector_oversample=self.detector_oversampling,
            fft_oversample=self.fft_oversampling,
            fov_arcsec=self.field_of_view,
            nlambda=self.nlambda,
            display=True
        )

        for widgetname in ['Display PSF', 'Display profiles', 'Save PSF As...']:
            self.widgets[widgetname].state(['!disabled'])
        _refresh_window()
        _log.info("PSF calculation complete")

    def ev_display_psf(self):
        """Event handler for displaying the PSF"""
        #self._update_from_gui()
        #if self.PSF_HDUlist is not None:
        plt.clf()
        poppy.display_PSF(
            self.psf_hdulist,
            vmin=self.advanced_options['psf_vmin'],
            vmax=self.advanced_options['psf_vmax'],
            scale=self.advanced_options['psf_scale'],
            cmap=self.advanced_options['psf_cmap'],
            normalize=self.advanced_options['psf_normalize']
        )
        _refresh_window()

    def ev_display_profiles(self):
        """Event handler for displaying PSF encircled energy profiles"""
        #self._update_from_gui()
        poppy.display_profiles(self.psf_hdulist)
        _refresh_window()

    def ev_display_optics(self):
        """Event handler for displaying the optical system"""
        self._update_from_gui()
        plt.clf()
        self.inst.display()
        _refresh_window()

    def ev_display_opd(self):
        """Event handler to display the OPD map"""
        self._update_from_gui()
        if self.inst.pupilopd is None:
            tkMessageBox.showwarning(
                message="You currently have selected no OPD file (i.e. "
                        "perfect telescope) so there's nothing to display.",
                title="Can't Display"
            )
        else:
            opd = fits.getdata(self.inst.pupilopd[0])

            if len(opd.shape) > 2:
                opd = opd[self.opd_i, :, :]  # grab correct slice

            # mask out all pixels which are exactly 0, outside the aperture
            masked_opd = np.ma.masked_equal(opd, 0)
            cmap = matplotlib.cm.jet
            cmap.set_bad('k', 0.8)
            plt.clf()
            plt.imshow(masked_opd, cmap=cmap, interpolation='nearest',
                       vmin=-0.5, vmax=0.5)
            plt.title("OPD from %s, #%d" %
                      (os.path.basename(self.opd_name), self.opd_i))
            cbar = plt.colorbar(orientation='vertical')
            cbar.set_label('microns')

            fig = plt.gcf()
            plt.text(0.4, 0.02,
                     "OPD WFE = %6.2f nm RMS" % (masked_opd.std() * 1000.),
                     transform=fig.transFigure)
            _refresh_window()

    def ev_update_opd_labels(self):
        """Update the descriptive text for all OPD files"""
        for iname in self.instrument.keys():
            self.ev_update_opd_label(
                self.widgets[iname+"_opd"],
                self.widgets[iname+"_opd_label"],
                iname
            )

    def ev_update_opd_label(self, widget_combobox, widget_label, iname):
        """Update the descriptive text displayed about one OPD file"""
        if widget_combobox.get() == self.PERFECT_OPD:
            header_summary = self.PERFECT_OPD
        else:
            filename = os.path.join(self.instrument[iname]._datapath, 'OPD',
                                    widget_combobox.get())
            header_summary = fits.getheader(filename)['SUMMARY']

        self.widgets[iname+"_opd_i"]['state'] = 'readonly'
        widget_label.configure(text=header_summary, width=30)

    def _update_from_gui(self):
        #TODO:jlong: this should update source and calc options, then subclasses only do instruments
        raise NotImplementedError("Subclasses must implement _update_from_gui")


class WebbPSFGUI(PSFGenerationGUI):
    """ A GUI for the Webb PSF Simulator

    Documentation TBD!

    """
    OUTPUT_OPTS = Options(
        oversampled='Oversampled image',
        detector_scale='Detector sampled image',
        both_as_fits='Both as FITS extensions',
        mock_jwst_dms='Mock JWST DMS Output'
    )

    def __init__(self, opdserver=None):
        super(WebbPSFGUI, self).__init__()
        # this list defines the order of instrument tabs in the GUI
        self.instrument_names = ['NIRCam', 'NIRSpec', 'NIRISS', 'MIRI', 'FGS']
        for i in self.instrument_names:
            self.instrument[i] = webbpsf_core.Instrument(i)

        #invoke link to ITM server if provided?
        if opdserver is not None:
            self._enable_opdserver = True
            self._opdserver = opdserver
        else:
            self._enable_opdserver = False

        # create widgets & run
        self._create_widgets()
        self.root.update()

    def _populate_instrument_config(self, notebook):
        for i, iname in enumerate(self.instrument_names):
            page = ttk.Frame(notebook, name="tab_" + iname)
            if i == 0:
                self._deleteme_page = page
            notebook.add(page, text=iname)
            notebook.select(i)  # make it active
            self.widgets[notebook.select()] = iname  # save reverse lookup to string name

            if iname == 'NIRCam':
                lframe = ttk.Frame(page)
                self._deleteme_frame = lframe

                instr_config_label = ttk.Label(
                    lframe,
                    text='Configuration Options for '+iname+', module: '
                )
                instr_config_label.grid(row=0, column=0, sticky='W')

                mname = 'NIRCam module'
                self.vars[mname] = tk.StringVar()
                self.widgets[mname] = ttk.Combobox(lframe, textvariable=self.vars[mname],
                                                   width=2, state='readonly')
                self.widgets[mname].grid(row=0,column=1, sticky='W')
                self.widgets[mname]['values'] = ['A', 'B']
                self.widgets[mname].set('A')

                lframe.grid(row=0, columnspan=2, sticky='W')
            else:
                instr_config_label = ttk.Label(page, text='Configuration Options for '+iname)
                instr_config_label.grid(row=0, columnspan=2, sticky='W')

            ttk.Button(page, text='Display Optics', command=self.ev_display_optics ).grid(column=2, row=0, sticky='E', columnspan=3)

            _add_labeled_dropdown(
                self,
                iname+"_filter",
                page,
                label='    Filter:',
                values=self.instrument[iname].filter_list,
                default=self.instrument[iname].filter,
                width=12,
                position=(1, 0),
                sticky='W'
            )

            if iname == 'NIRSpec' or iname =='MIRI':
                fr2 = ttk.Frame(page)
                #label = 'IFU' if iname !='TFI' else 'TF'
                ttk.Label(fr2, text='   IFU wavelen: ', state='disabled').grid(row=0, column=0)
                self.widgets[iname+"_ifu_wavelen"] = ttk.Entry(fr2, width=5) #, disabledforeground="#A0A0A0")
                self.widgets[iname+"_ifu_wavelen"].insert(0, str(self.instrument[iname].monochromatic))
                self.widgets[iname+"_ifu_wavelen"].grid(row=0, column=1)
                self.widgets[iname+"_ifu_wavelen"].state(['disabled'])
                ttk.Label(fr2, text=' um' , state='disabled').grid(row=0, column=2)
                fr2.grid(row=1,column=2, columnspan=6, sticky='E')

                iname2 = iname+"" # need to make a copy so the following lambda function works right:
                self.widgets[iname+"_filter"].bind('<<ComboboxSelected>>', lambda e: self.ev_update_ifu_label(iname2))

            if len(self.instrument[iname].image_mask_list) >0 :
                masks = self.instrument[iname].image_mask_list
                masks.insert(0, "")

                _add_labeled_dropdown(self, iname+"_coron", page, label='    Coron:', values=masks,  width=12, position=(2,0), sticky='W')

            if len(self.instrument[iname].image_mask_list) >0 :
                masks = self.instrument[iname].pupil_mask_list
                masks.insert(0, "")
                _add_labeled_dropdown(self, iname+"_pupil", page, label='    Pupil:', values=masks,  width=12, position=(3,0), sticky='W')

                fr2 = ttk.Frame(page)
                _add_labeled_entry(self, iname+"_pupilshift_x", fr2, label='  pupil shift in X:', value='0', width=3, position=(3,4), sticky='W')
                _add_labeled_entry(self, iname+"_pupilshift_y", fr2, label=' Y:', value='0', width=3, position=(3,6), sticky='W')

                ttk.Label(fr2, text='% of pupil' ).grid(row=3, column=8)
                fr2.grid(row=3, column=3, sticky='W')

            ttk.Label(page, text='Configuration Options for the OTE').grid(row=4, columnspan=2, sticky='W')
            fr2 = ttk.Frame(page)

            opd_list = self.instrument[iname].opd_list
            opd_list.insert(0, self.PERFECT_OPD)

            if self._enable_opdserver:
                opd_list.append("OPD from ITM Server")
            default_opd = self.instrument[iname].pupilopd if self.instrument[iname].pupilopd is not None else self.PERFECT_OPD
            _add_labeled_dropdown(self, iname+"_opd", fr2, label='    OPD File:', values=opd_list, default=default_opd, width=21, position=(0,0), sticky='W')

            _add_labeled_dropdown(self, iname+"_opd_i", fr2, label=' # ', values= [str(i) for i in range(10)], width=3, position=(0,2), sticky='W')

            self.widgets[iname+"_opd_label"] = ttk.Label(fr2, text=' 0 nm RMS            ', width=35)
            self.widgets[iname+"_opd_label"].grid( column=4,sticky='W', row=0)

            self.widgets[iname+"_opd"].bind(
                '<<ComboboxSelected>>',
                lambda e, inst=iname: self.ev_update_opd_label(
                    self.widgets[inst+"_opd"],
                    self.widgets[inst+"_opd_label"],
                    inst
                )
            )
            ttk.Button(fr2, text='Display', command=self.ev_display_opd).grid(column=5,sticky='E',row=0)

            fr2.grid(row=5, column=0, columnspan=4, sticky='S')

            # ITM interface here - build the widgets now but they will be hidden by default until the ITM option is selected
            fr2 = ttk.Frame(page)
            _add_labeled_entry(self, iname+"_coords", fr2, label='    Source location:', value='0, 0', width=12, position=(1,0), sticky='W')
            units_list = ['V1,V2 coords', 'detector pixels']
            _add_labeled_dropdown(self, iname+"_coord_units", fr2, label='in:', values=units_list, default=units_list[0], width=11, position=(1,2), sticky='W')
            choose_list=['', 'SI center', 'SI upper left corner', 'SI upper right corner', 'SI lower left corner', 'SI lower right corner']
            _add_labeled_dropdown(self, iname+"_coord_choose", fr2, label='or select:', values=choose_list, default=choose_list[0], width=21, position=(1,4), sticky='W')

            ttk.Label(fr2, text='    ITM output:').grid(row=2, column=0, sticky='W')
            self.widgets[iname+"_itm_output"] = ttk.Label(fr2, text='    - no file available yet -')
            self.widgets[iname+"_itm_output"].grid(row=2, column=1, columnspan=4, sticky='W')
            ttk.Button(fr2, text='Access ITM...', command=self.ev_launch_itm_dialog).grid(column=5,sticky='E',row=2)


            fr2.grid(row=6, column=0, columnspan=4, sticky='SW')
            self.widgets[iname+"_itm_coords"] = fr2

    def _update_from_gui(self):
        # get GUI values
        if _HAS_PYSYNPHOT_DATA:
            self.sptype = self.widgets['SpType'].get()
        self.iname = self.widgets[self.widgets['tabset'].select()]

        nlambda_string = self.widgets['nlambda'].get()
        if len(nlambda_string) > 0:
            try:
                self.nlambda = int()
            except ValueError:
                tkMessageBox.showwarning(
                    message="Got a non-integer value for number of wavelengths. Falling back on"
                            "autoselected value.",
                    title="Invalid Number of Wavelengths"
                )
                self.nlambda = None
        else:
            self.nlambda = None

        self.field_of_view = float(self.widgets['FOV'].get())
        self.fft_oversampling = int(self.widgets['fft_oversampling'].get())
        self.detector_oversampling = int(self.widgets['detector_oversampling'].get())

        self.output_format = self.widgets['output_format'].get()

        options = {
            'rebin': not (self.output_format == 'Oversampled PSF only'),
            'mock_dms': (self.output_format == 'Mock full image from JWST DMS'),
            'jitter': self.widgets['jitter'].get(),
        }

        # and get the values that may have previously been set by the 'advanced options' dialog
        if self.advanced_options is not None:
            for a in self.advanced_options.keys():
                options[a] = self.advanced_options[a]

        # configure the relevant instrument object
        self.inst = self.instrument[self.iname]
        self.filter = self.widgets[self.iname+"_filter"].get()
        self.inst.filter = self.filter

        self.opd_name = self.widgets[self.iname+"_opd"].get()
        if self._enable_opdserver and 'ITM' in self.opd_name:
            # get from ITM server
            self.opd_i = 0
            self.inst.pupilopd = self._opdserver.get_OPD(return_as="FITS")
            self.opd_name = "OPD from ITM OPD GUI"

        elif self.opd_name == self.PERFECT_OPD:
            # perfect OPD
            self.opd_name = "Perfect"
            self.inst.pupilopd = None
        else:
            # Regular FITS file version
            self.opd_name = self.widgets[self.iname+"_opd"].get()
            self.opd_i = int(self.widgets[self.iname+"_opd_i"].get())
            # OPD tuple as (path, slice number)
            self.inst.pupilopd = (
                os.path.join(self.inst._datapath, "OPD", self.opd_name),
                self.opd_i
            )

        _log.info("Selected OPD is "+str(self.opd_name))

        if self.iname+"_coron" in self.widgets:
            self.inst.image_mask = self.widgets[self.iname+"_coron"].get()
            self.inst.pupil_mask = self.widgets[self.iname+"_pupil"].get()
            # TODO read in mis-registration options here.

            options['source_offset_r'] = float(self.widgets["source_off_r"].get())
            options['source_offset_theta'] = float(self.widgets["source_off_theta"].get())

            # convert percentages to fractions
            shift_x = self.widgets[self.iname+"_pupilshift_x"].get()
            try:
                options['pupil_shift_x'] = float(shift_x) / 100.0
            except ValueError:
                tkMessageBox.showwarning(
                    message="Got an invalid value for pupil X shift. Using 0% instead.",
                    title="Invalid Pupil X Shift"
                )
                options['pupil_shift_x'] = 0.0

            shift_y = self.widgets[self.iname+"_pupilshift_y"].get()
            try:
                options['pupil_shift_y'] = float(shift_y) / 100.0
            except ValueError:
                tkMessageBox.showwarning(
                    message="Got an invalid value for pupil X shift. Using 0% instead.",
                    title="Invalid Pupil X Shift"
                )
                options['pupil_shift_y'] = 0.0

        self.inst.options = options

    def ev_launch_itm_dialog(self):
        tkMessageBox.showwarning(
            message="ITM dialog box not yet implemented",
            title="Can't Display"
        )

    def ev_update_opd_label(self, widget_combobox, widget_label, iname):
        """Update the descriptive text displayed about one OPD file"""
        if 'ITM' in widget_combobox.get():
            if self._enable_opdserver:
                header_summary = "Get OPD from ITM Server"
                self.widgets[iname+"_itm_coords"].grid() # re-show ITM options
            else:
                header_summary = ("ITM Server is not running "
                                  "or otherwise unavailable.")
                self.widgets[iname+"_itm_coords"].grid_remove()  # hide ITM options

            widget_label.configure(text=header_summary, width=30)
            self.widgets[iname+"_itm_coords"].grid() # re-show ITM options
        else:
            super(WebbPSFGUI, self).ev_update_opd_label(widget_combobox, widget_label, iname)
            self.widgets[iname+"_itm_coords"].grid_remove()  # hide ITM options

#-------------------------------------------------------------------------

class OptionsDialog(object, tk.Toplevel):
    """
    PSF generation GUI options dialog

    Based on example code from
    <http://effbot.org/tkinterbook/tkinter-dialog-windows.htm>
    """

    def __init__(self, parent, title=None, input_options=None):
        self.results = input_options  # in case we cancel this gets returned
        self.vars = {}
        self.widgets = {}

        self.colortables = {
            'Jet (blue to red)': matplotlib.cm.jet,
            'Gray': matplotlib.cm.gray,
            'Heat (black-red-yellow)': matplotlib.cm.gist_heat,
            'Copper (black to tan)': matplotlib.cm.copper,
            'Stern': matplotlib.cm.gist_stern,
            'Prism (repeating rainbow)': matplotlib.cm.prism,
        }

        tk.Toplevel.__init__(self, parent)
        self.transient(parent)
        self.input_options = input_options

        if title:
            self.title(title)

        self.parent = parent
        self.result = None

        body = ttk.Frame(self)
        self.initial_focus = self.body(body)
        body.pack(padx=5, pady=5)

        self.buttonbox()
        self.grab_set()
        if not self.initial_focus:
            self.initial_focus = self

        self.protocol("WM_DELETE_WINDOW", self.cancel)
        self.geometry("+%d+%d" % (parent.winfo_rootx()+50,
                                  parent.winfo_rooty()+50))

        self.initial_focus.focus_set()
        self.wait_window(self)
    #
    # construction hooks

    def body(self, master):
        # create dialog body.  return widget that should have
        # initial focus.  this method should be overridden
        pass

    def buttonbox(self):
        # add standard button box. override if you don't want the
        # standard buttons

        box = ttk.Frame(self)

        w = ttk.Button(box, text="OK", width=10,
                       command=self.ok, default=tk.ACTIVE)
        w.pack(side=tk.LEFT, padx=5, pady=5)
        w = ttk.Button(box, text="Cancel", width=10,
                       command=self.cancel)
        w.pack(side=tk.LEFT, padx=5, pady=5)

        self.bind("<Return>", self.ok)
        self.bind("<Escape>", self.cancel)

        box.pack()

    #
    # standard button semantics

    def ok(self, event=None):
        if not self.validate():
            self.initial_focus.focus_set() # put focus back
            return

        self.withdraw()
        self.update_idletasks()

        self.apply()

        self.cancel()

    def cancel(self, event=None):
        # put focus back to the parent window
        self.parent.focus_set()
        self.destroy()

    #
    # command hooks

    def validate(self):
        return self.apply(test=True)

    def apply(self, test=False):
        raise NotImplementedError("Dialog subclasses must implement apply()")

PROPAGATION_OPTS = Options(
    mft='regular propagation (MFT)',
    full='full coronographic propagation (FFT/SAM)'
)
SEMI_ANALYTIC_OPTS = Options(
    semi_analytic='semi-analytic method if possible',
    basic_fft='basic FFT method always'
)
PARITY_OPTS = Options(
    odd='odd',
    even='even',
    either='either'
)


class WebbPSFOptionsDialog(OptionsDialog):
    def __init__(self, *args, **kwargs):
        self.results = None

        super(WebbPSFOptionsDialog, self).__init__(*args, **kwargs)

    def body(self, master):
        calc_options_labelframe = ttk.LabelFrame(master, text='WebbPSF Advanced Options')

        row_counter = 0
        _add_labeled_dropdown(
            self,
            "force_coron",
            calc_options_labelframe,
            label='Direct imaging calculations use',
            values=PROPAGATION_OPTS,
            default=(PROPAGATION_OPTS.full if self.input_options['force_coron']
                     else PROPAGATION_OPTS.mft),
            width=30,
            position=(row_counter, 0),
            sticky='W'
        )
        row_counter += 1
        _add_labeled_dropdown(
            self,
            "no_sam",
            calc_options_labelframe,
            label='Coronagraphic calculations use',
            values=SEMI_ANALYTIC_OPTS,
            default=(SEMI_ANALYTIC_OPTS.basic_fft if self.input_options['no_sam']
                     else SEMI_ANALYTIC_OPTS.semi_analytic),
            width=30,
            position=(row_counter, 0),
            sticky='W'
        )
        row_counter += 1
        _add_labeled_dropdown(
            self,
            "parity",
            calc_options_labelframe,
            label='Output pixel grid parity is',
            values=PARITY_OPTS,
            default=self.input_options['parity'],
            width=10,
            position=(row_counter, 0),
            sticky='W'
        )

        calc_options_labelframe.grid(row=1, sticky='E,W', padx=10, pady=5)

        return self.widgets['force_coron']  # return widget to have initial focus

    def apply(self, test=False):
        results = {
            'force_coron': self.vars['force_coron'].get() == PROPAGATION_OPTS.full,
            'no_sam': self.vars['no_sam'].get() == SEMI_ANALYTIC_OPTS.basic,
            'parity': self.vars['parity'].get(),
            'psf_scale': self.vars['psf_scale'].get(),
            'psf_cmap_str': self.vars['psf_cmap'].get(),
            'psf_normalize': self.vars['psf_normalize'].get(),
            'psf_cmap': self.colortables[str(self.vars['psf_cmap'].get())],
        }

        try:
            results['psf_vmax'] = float(self.vars['psf_vmax'].get())
            results['psf_vmin'] = float(self.vars['psf_vmin'].get())
        # TODO:jlong: validate interactively
        # http://stackoverflow.com/questions/4140437/
        except ValueError:
            tkMessageBox.showwarning(
                message="Ensure the PSF scale values are valid floating-point literals (e.g. "
                        "'0.005', '1.2e-5')",
                title="Invalid PSF Plot Scale Values"
            )
            return False

        if not test:
            self.results = results
        return True

#-------------------------------------------------------------------------


def synplot(thing, waveunit='micron', label=None, **kwargs):
    """
    Plot a single PySynPhot object (either SpectralElement or SourceSpectrum)
    versus wavelength, with nice axes labels.

    Really just a simple convenience function.
    """

    # convert to requested display unit.
    wave = thing.waveunits.Convert(thing.wave, waveunit)

    if label is None:
        label = thing.name

    if isinstance(thing, pysynphot.spectrum.SourceSpectrum):
        artist = plt.loglog(wave, thing.flux, label=label, **kwargs)
        plt.xlabel("Wavelength [%s]" % waveunit)
        if str(thing.fluxunits) == 'flam':
            plt.ylabel("Flux [%s]" % ' erg cm$^{-2}$ s$^{-1}$ Ang$^{-1}$')
        else:
            plt.ylabel("Flux [%s]" % thing.fluxunits)
    elif isinstance(thing, pysynphot.spectrum.SpectralElement):
        artist = plt.plot(wave, thing.throughput, label=label, **kwargs)
        plt.xlabel("Wavelength [%s]" % waveunit)
        plt.ylabel("Throughput")
        plt.gca().set_ylim(0, 1)
    else:
        _log.error("Don't know how to plot that object...")
        artist = None
    return artist


def _refresh_window():
    """
    Force the window to refresh, and optionally to show itself
    if hidden (for recent matplotlibs)
    """
    plt.draw()
    plt.show(block=False)


def tkgui():
    """Launch the Tk GUI to WebbPSF"""
    # enable log message printout
    logging.basicConfig(
        level=logging.INFO,
        format='%(name)-10s: %(levelname)-8s %(message)s'
    )

    if _HAS_PYSYNPHOT_INSTALLED and not _HAS_PYSYNPHOT_DATA:
        _log.warn(
            "Environment variable PYSYN_CDBS not set in environment! "
            "Set PYSYN_CDBS to the directory containing the CDBS data. "
            "For more information, see `Installing or updating pysynphot' in "
            "installation.rst."
        )

    # start the GUI
    gui = WebbPSFGUI()
    gui.mainloop()
