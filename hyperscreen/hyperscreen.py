#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""hyperscreen is an improved background rejection algorithm
for """

from astropy.io import fits
from astropy.table import Table

import pandas as pd
import numpy as np
np.seterr(divide='ignore')

import matplotlib as mpl


import warnings
warnings.filterwarnings("ignore",category =RuntimeWarning)

class HRCevt1:
    '''
    A more robust HRC EVT1 file. Includes explicit 
    columns for every status bit, as well as calculated 
    columns for the f_p, f_b plane for your boomerangs. 
    Check out that cool new filtering algorithm!
    '''

    def __init__(self, evtfile, as_dataframe = True):
        self.filename = evtfile
        self.hdulist = fits.open(evtfile)
        self.data = Table(self.hdulist[1].data)
        self.header = self.hdulist[1].header
        self.gti = self.hdulist[2].data
        self.hdulist.close() # Don't forget to close your fits file! 

        fp_u, fb_u, fp_v, fb_v = self.calculate_fp_fb()
        
        self.gti.starts = self.gti['START']
        self.gti.stops = self.gti['STOP']

        self.data["fp_u"] = fp_u
        self.data["fb_u"] = fb_u
        self.data["fp_v"] = fp_v
        self.data["fb_v"] = fb_v

        # Make individual status bit columns with legible names
        self.data["AV3 corrected for ringing"] = self.data["status"][:,0]
        self.data["AU3 corrected for ringing"] = self.data["status"][:,1]
        self.data["Event impacted by prior event (piled up)"] = self.data["status"][:,2]
        # Bit 4 (Python 3) is spare
        self.data["Shifted event time"] = self.data["status"][:,4]
        self.data["Event telemetered in NIL mode"] = self.data["status"][:,5]
        self.data["V axis not triggered"] = self.data["status"][:,6]
        self.data["U axis not triggered"] = self.data["status"][:,7]
        self.data["V axis center blank event"] = self.data["status"][:,8]
        self.data["U axis center blank event"] = self.data["status"][:,9]
        self.data["V axis width exceeded"] = self.data["status"][:,10]
        self.data["U axis width exceeded"] = self.data["status"][:,11]
        self.data["Shield PMT active"] = self.data["status"][:,12]
        # Bit 14 (Python 13) is hardware spare
        self.data["Upper level discriminator not exceeded"] = self.data["status"][:,14]
        self.data["Lower level discriminator not exceeded"] = self.data["status"][:,15]
        self.data["Event in bad region"] = self.data["status"][:,16]
        self.data["Amp total on V or U = 0"] = self.data["status"][:,17]
        self.data["Incorrect V center"] = self.data["status"][:,18]
        self.data["Incorrect U center"] = self.data["status"][:,19]
        self.data["PHA ratio test failed"] = self.data["status"][:,20]
        self.data["Sum of 6 taps = 0"] = self.data["status"][:,21]
        self.data["Grid ratio test failed"] = self.data["status"][:,22]
        self.data["ADC sum on V or U = 0"] = self.data["status"][:,23]
        self.data["PI exceeding 255"] = self.data["status"][:,24]
        self.data["Event time tag is out of sequence"] = self.data["status"][:,25]
        self.data["V amp flatness test failed"] = self.data["status"][:,26]
        self.data["U amp flatness test failed"] = self.data["status"][:,27]
        self.data["V amp saturation test failed"] = self.data["status"][:,28]
        self.data["U amp saturation test failed"] = self.data["status"][:,29]
        self.data["V hyperbolic test failed"] = self.data["status"][:,30]
        self.data["U hyperbolic test failed"] = self.data["status"][:,31]
        self.data["Hyperbola test passed"] = np.logical_not(np.logical_or(self.data['U hyperbolic test failed'], self.data['V hyperbolic test failed']))
        self.data["Hyperbola test failed"] = np.logical_or(self.data['U hyperbolic test failed'], self.data['V hyperbolic test failed'])


        self.obsid = self.header["OBS_ID"]
        self.obs_date = self.header["DATE"]
        self.target = self.header["OBJECT"]
        self.detector = self.header["DETNAM"]
        self.grating = self.header["GRATING"]
        self.exptime = self.header["EXPOSURE"]


        self.numevents = len(self.data["time"])
        self.hyperbola_passes = np.sum(np.logical_or(self.data['U hyperbolic test failed'], self.data['V hyperbolic test failed']))
        self.hyperbola_failures = np.sum(np.logical_not(np.logical_or(self.data['U hyperbolic test failed'], self.data['V hyperbolic test failed'])))

        if self.hyperbola_passes + self.hyperbola_failures != self.numevents:
            print("Warning: Number of Hyperbola Test Failures and Passes ({}) does not equal total number of events ({}).".format(self.hyperbola_passes + self.hyperbola_failures, self.numevents))

        if as_dataframe is True:
            # Multidimensional columns don't grok with Pandas
            self.data.remove_column('status')
            self.data = self.data.to_pandas()


    def __str__(self):
        if isinstance(self.data, pd.core.frame.DataFrame):
            return "HRC EVT1 object with {} events. Data is packaged as a Pandas Dataframe".format(self.numevents)
        else: 
            return "HRC EVT1 object with {} events. Data is packaged as an Astropy Table".format(self.numevents)


    def calculate_fp_fb(self):
        '''
        Calculate the Fine Position (fp) and normalized central tap
        amplitude (fb) for the HRC U- and V- axes.

        Parameters
        ----------
        data : Astropy Table
            Table object made from an HRC evt1 event list. Must include the
            au1, au2, au3 and av1, av2, av3 columns.

        Returns
        -------
        fp_u, fb_u, fp_v, fb_v: float
            Calculated fine positions and normalized central tap amplitudes
            for the HRC U- and V- axes
        '''
        a_u = self.data["au1"]  # otherwise known as "a1"
        b_u = self.data["au2"]  # "a2"
        c_u = self.data["au3"]  # "a3"

        a_v = self.data["av1"]
        b_v = self.data["av2"]
        c_v = self.data["av3"]

        with np.errstate(invalid='ignore'):
            # Do the U axis
            fp_u = ((c_u - a_u) / (a_u + b_u + c_u))
            fb_u = b_u / (a_u + b_u + c_u)

            # Do the V axis
            fp_v = ((c_v - a_v) / (a_v + b_v + c_v))
            fb_v = b_v / (a_v + b_v + c_v)

        return fp_u, fb_u, fp_v, fb_v


    def threshold(self, img, bins):
        nozero_img = img.copy()
        nozero_img[img == 0] = np.nan
        
        # This is a really stupid way to threshold
        median = np.nanmedian(nozero_img)
        thresh = median*5
        
        thresh_img = nozero_img
        thresh_img[thresh_img < thresh] = np.nan
        thresh_img[:int(bins[1]/2),:] = np.nan
    #     thresh_img[:,int(bins[1]-5):] = np.nan
        return thresh_img, thresh

    def tapscreen(self):
        '''
        Grant Tremblay's new algorithm. Screens events on a tap-by-tap basis.
        '''
        if not isinstance(self.data, pd.core.frame.DataFrame):
            raise ValueError("This HRCevt1 object is not a Pandas Dataframe. You must initialize it with 'as_dataframe=True'")

        data = self.data

        taprange = range(2,58)

        bins = [200,200] # number of bins
        
        # Instantiate these empty dictionaries to hold our results
        u_axis_survivals = {}
        v_axis_survivals = {}
        
        for tap in taprange:
            
            tapmask_u = data[data['crsu'] == tap].index.values
            tapmask_v = data[data['crsv'] == tap].index.values

            # If either of these are zero, skip this iteration of the for loop
            # if len(tapmask_u) or len(tapmask_v) == 0:
            #     continue
            
            # DO U AXIS
            keep_u = np.isfinite(data['fb_u'][tapmask_u])
            
            hist_u, xbounds_u, ybounds_u = np.histogram2d(data['fb_u'][tapmask_u][keep_u], data['fp_u'][tapmask_u][keep_u], bins=bins)
            thresh_hist_u, thresh_u = self.threshold(hist_u, bins=bins)
            
            posx_u = np.digitize(data['fb_u'][tapmask_u], xbounds_u)
            posy_u = np.digitize(data['fp_u'][tapmask_u], ybounds_u)
            hist_mask_u = (posx_u > 0) & (posx_u <= bins[0]) & (posy_u > -1) & (posy_u <= bins[1])
            
            hhsub_u = thresh_hist_u[posx_u[hist_mask_u] - 1, posy_u[hist_mask_u] - 1] # Values of the histogram where the points are
            pass_fb_u = data['fb_u'][tapmask_u][hist_mask_u][np.isfinite(hhsub_u)]
            
            u_axis_survivals["U Axis Tap {:02d}".format(tap)] = pass_fb_u.index.values
            
            
            # DO V AXIS
            
            keep_v = np.isfinite(data['fb_v'][tapmask_v])
            
            hist_v, xbounds_v, ybounds_v = np.histogram2d(data['fb_v'][tapmask_v][keep_v], data['fp_v'][tapmask_v][keep_v], bins=bins)
            thresh_hist_v, thresh_v = self.threshold(hist_v, bins=bins)
            
            posx_v = np.digitize(data['fb_v'][tapmask_v], xbounds_v)
            posy_v = np.digitize(data['fp_v'][tapmask_v], ybounds_v)
            hist_mask_v = (posx_v > 0) & (posx_v <= bins[0]) & (posy_v > -1) & (posy_v <= bins[1])
            
            hhsub_v = thresh_hist_v[posx_v[hist_mask_v] - 1, posy_v[hist_mask_v] - 1] # Values of the histogram where the points are
            pass_fb_v = data['fb_v'][tapmask_v][hist_mask_v][np.isfinite(hhsub_v)]
            
            v_axis_survivals["V Axis Tap {:02d}".format(tap)] = pass_fb_v.index.values
            if len(v_axis_survivals) == 0:
                print("HMMM V")           
        
        # print(v_axis_survivals)

        u_all_survivals = np.concatenate([x for x in u_axis_survivals.values()])
        v_all_survivals = np.concatenate([x for x in v_axis_survivals.values()])

        all_survivals = np.intersect1d(u_all_survivals, v_all_survivals)
        survival_mask = np.isin(self.data.index.values, all_survivals)
        failure_mask = np.logical_not(survival_mask)

        num_survivals = sum(survival_mask)
        num_failures = sum(failure_mask)

        percent_tapscreen_rejected = round(((num_failures / self.numevents) * 100),2)

        if num_survivals + num_failures != self.numevents:
            print("WARNING!!! TOTAL NUMBER OF FAILURES AND SURVIVALS DOES NOT EQUAL TOTAL EVENTS. SOMETHING IS WRONG!")

        legacy_hyperbola_test_survivals = sum(self.data['Hyperbola test passed'])
        legacy_hyperbola_test_failures = sum(self.data['Hyperbola test failed'])
        percent_legacy_hyperbola_test_rejected = round(((legacy_hyperbola_test_failures / self.numevents) * 100),2)

        percent_improvement_over_legacy_test = round((percent_tapscreen_rejected - percent_legacy_hyperbola_test_rejected),2)

        tapscreen_results_dict = {"U Axis Survivals by Tap": u_axis_survivals,
                          "V Axis Survivals by Tap": v_axis_survivals,
                          "U Axis All Survivals": u_all_survivals,
                          "V Axis All Survivals": v_all_survivals,
                          "All Survivals (event indices)": all_survivals,
                          "All Survivals (boolean mask)": survival_mask,
                          "All Failures (boolean mask)": failure_mask,
                          "Percent rejected by Tapscreen": percent_tapscreen_rejected,
                          "Percent rejected by Hyperbola": percent_legacy_hyperbola_test_rejected,
                          "Percent improvement": percent_improvement_over_legacy_test
                          }

        return tapscreen_results_dict



    def hyperbola(self,fb, a, b, h):
        '''Create a simple hyperbola'''
        hyperbola = b * np.sqrt(((fb - h)**2 / a**2) - 1)

        return hyperbola


    def legacy_hyperbola_test(self, tolerance=0.035):
        '''
        Apply the hyperbolic test.
        '''

        # Remind the user what tolerance they're using
        # print("{0: <25}| Using tolerance = {1}".format(" ", tolerance))

        # Set hyperbolic coefficients, depending on whether this is HRC-I or -S
        if self.detector == "HRC-I":
            a_u = 0.3110
            b_u = 0.3030
            h_u = 1.0580

            a_v = 0.3050
            b_v = 0.2730
            h_v = 1.1
            # print("{0: <25}| Using HRC-I hyperbolic coefficients: ".format(" "))
            # print("{0: <25}|    Au={1}, Bu={2}, Hu={3}".format(" ", a_u, b_u, h_u))
            # print("{0: <25}|    Av={1}, Bv={2}, Hv={3}".format(" ", a_v, b_v, h_v))

        if self.detector == "HRC-S":
            a_u = 0.2706
            b_u = 0.2620
            h_u = 1.0180

            a_v = 0.2706
            b_v = 0.2480
            h_v = 1.0710
            # print("{0: <25}| Using HRC-S hyperbolic coefficients: ".format(" "))
            # print("{0: <25}|    Au={1}, Bu={2}, Hu={3}".format(" ", a_u, b_u, h_u))
            # print("{0: <25}|    Av={1}, Bv={2}, Hv={3}".format(" ", a_v, b_v, h_v))

        # Set the tolerance boundary ("width" of the hyperbolic region)

        h_u_lowerbound = h_u * (1 + tolerance)
        h_u_upperbound = h_u * (1 - tolerance)

        h_v_lowerbound = h_v * (1 + tolerance)
        h_v_upperbound = h_v * (1 - tolerance)

        # Compute the Hyperbolae
        with np.errstate(invalid='ignore'):
            zone_u_fit = self.hyperbola(self.data["fb_u"], a_u, b_u, h_u)
            zone_u_lowerbound = self.hyperbola(self.data["fb_u"], a_u, b_u, h_u_lowerbound)
            zone_u_upperbound = self.hyperbola(self.data["fb_u"], a_u, b_u, h_u_upperbound)

            zone_v_fit = self.hyperbola(self.data["fb_v"], a_v, b_v, h_v)
            zone_v_lowerbound = self.hyperbola(self.data["fb_v"], a_v, b_v, h_v_lowerbound)
            zone_v_upperbound = self.hyperbola(self.data["fb_v"], a_v, b_v, h_v_upperbound)

        zone_u = [zone_u_lowerbound, zone_u_upperbound]
        zone_v = [zone_v_lowerbound, zone_v_upperbound]

        # Apply the masks
        # print("{0: <25}| Hyperbolic masks for U and V axes computed".format(""))

        with np.errstate(invalid='ignore'):
            # print("{0: <25}| Creating U-axis mask".format(""), end=" |")
            between_u = np.logical_not(np.logical_and(
                self.data["fp_u"] < zone_u[1], self.data["fp_u"] > -1 * zone_u[1]))
            not_beyond_u = np.logical_and(
                self.data["fp_u"] < zone_u[0], self.data["fp_u"] > -1 * zone_u[0])
            condition_u_final = np.logical_and(between_u, not_beyond_u)

            # print(" Creating V-axis mask")
            between_v = np.logical_not(np.logical_and(
                self.data["fp_v"] < zone_v[1], self.data["fp_v"] > -1 * zone_v[1]))
            not_beyond_v = np.logical_and(
                self.data["fp_v"] < zone_v[0], self.data["fp_v"] > -1 * zone_v[0])
            condition_v_final = np.logical_and(between_v, not_beyond_v)


        mask_u = condition_u_final
        mask_v = condition_v_final

        hyperzones = {"zone_u_fit": zone_u_fit,
                        "zone_u_lowerbound": zone_u_lowerbound,
                        "zone_u_upperbound": zone_u_upperbound,
                        "zone_v_fit": zone_v_fit,
                        "zone_v_lowerbound": zone_v_lowerbound,
                        "zone_v_upperbound": zone_v_upperbound}

        hypermasks = {"mask_u": mask_u, "mask_v": mask_v}

        # print("{0: <25}| Hyperbolic masks created".format(""))
        # print("{0: <25}| ".format(""))
        return hyperzones, hypermasks


    def get_gtimask(self):

        all_gti_indices = np.array([])

        for i in range(len(self.gti.starts)):
            goodtimes = np.where((self.data['time'] > self.gti.starts[i]) & (self.data['time'] < self.gti.stops[i]))
            all_gti_indices = np.concatenate([all_gti_indices, goodtimes[0]]).astype(int)

        gtimask = [True if i in all_gti_indices else False for i in range(len(self.data['time']))]

        self.gtimask = gtimask
        return self.gtimask



    def boomerang(self, mask=None, show=True, save=False, savedir='./'):

        self.fig, self.ax = plt.subplots(figsize=(12,8))

        if mask is not None:
            frame = self.ax.scatter(self.data['fb_u'][mask], self.data['fp_u'][mask], c=self.data['sumamps'][mask], cmap='plasma', s=0.5, rasterized=True)
        else:
            frame = self.ax.scatter(self.data['fb_u'], self.data['fp_u'], c=self.data['sumamps'], cmap='plasma', s=0.5, rasterized=True)

        self.ax.set_title('{} | {} | ObsID {} | {} ksec | {} counts'.format(self.target, self.detector, self.obsid, round(self.exptime/1000, 1), self.numevents))
        self.ax.set_ylabel(r'Fine Position $f_p$ $(C-A)/(A + B + C)$')
        self.ax.set_xlabel(r'Normalized Central Tap Amplitude $f_b$ $B / (A+B+C)$')

        self.cbar = plt.colorbar(frame, pad=-0.005)
        self.cbar.set_label("SUMAMPS")

        if show is True:
            plt.show()

        if save is True:
            savepath = '{}{}_{}_boomerang.pdf'.format(savedir, self.obsid, self.detector)
            self.fig.savefig(savepath, dpi=150, bbox_inches='tight')
            print('Saved boomerang figure to: {}'.format(savepath))


    def quicklook(self, masked_x=None, masked_y=None, title=None, show=True, return_img_data=False):

        nbins = (300,300)

        if masked_x is not None and masked_y is not None:
            x = masked_x[400:]
            y = masked_y[400:]
            img_data, yedges, xedges = np.histogram2d(y,x,nbins)
        else:
            x = self.data['x'][400:]
            y = self.data['y'][400:]
            img_data, yedges, xedges = np.histogram2d(y,x, nbins)

        extent = [xedges[0], xedges[-1], yedges[0], yedges[-1]]


        if show is True:
            self.fig, self.ax = plt.subplots(figsize=(8,8))
            self.ax.grid(False)
            plt.imshow(img_data, extent=extent, norm=LogNorm(), interpolation=None, cmap='viridis', origin='lower')
            #plt.imshow(img_data,  interpolation=None, cmap='magma', origin='lower')
            if title is None:
                self.ax.set_title("ObsID {} | {} | {} | {:,} events".format(self.obsid, self.target, self.detector, self.numevents))
            else:
                self.ax.set_title("{}".format(title))
            self.ax.set_xlabel("Sky X")
            self.ax.set_ylabel("Sky Y")
            plt.show(block=True)


        if return_img_data is True:
            return img_data, extent





def styleplots():

    mpl.rcParams['agg.path.chunksize'] = 10000

    # Make things pretty
    plt.style.use('ggplot')

    labelsizes = 15

    plt.rcParams['font.size'] = labelsizes
    plt.rcParams['axes.titlesize'] = 12
    plt.rcParams['axes.labelsize'] = labelsizes
    plt.rcParams['xtick.labelsize'] = labelsizes
    plt.rcParams['ytick.labelsize'] = labelsizes


def image(x, y, nbins=(300,300), truncate_start=False, cmap=mpl.cm.magma, title='Set Title', show=True, savepath=None):

        if truncate_start is True:
            x = x[400:]
            y = y[400:]

        img_data, yedges, xedges = np.histogram2d(y,x,nbins)

        extent = [xedges[0], xedges[-1], yedges[0], yedges[-1]]

        fig, ax = plt.subplots(figsize=(8,8))
        # cmap.set_bad('white')
        ax.imshow(img_data, norm=LogNorm(), interpolation=None, cmap=cmap, origin='lower')
        ax.set_title(title)
        ax.grid(False)
        ax.set_ylabel('y')
        ax.set_xlabel('x')

        if show is True:
            plt.show()

        if savepath is not None:
            plt.savefig(savepath, dpi=300)
            print("Saved plot to {}".format(savepath))

        plt.close()
