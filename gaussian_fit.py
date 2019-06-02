#!/usr/bin/python

from bitstring import *

from struct import *
import sys
from optparse import OptionParser
import matplotlib.pyplot as plt
#import progressbar as pb
import numpy as np
import time
import math
# local modules
#import gen_utils as util
import scipy.optimize as op
from scipy.stats import chisquare, kstest, ks_2samp, ttest_ind
from scipy.interpolate import splrep, splev
from scipy.signal import argrelextrema
import os.path
import glob
from numpy.polynomial import polynomial as P
import pylab
import psrchive
from astropy.modeling import models, fitting
import emcee
import corner
import utils


def CommandLineParse():
    usage = "usage: %prog -i <input file>"
    parser = OptionParser(usage)
    parser.add_option("-i", "--infile", action="store", type="string",
                      help="EPOS file to be read.",dest="infile")
    parser.add_option("-o", "--outfile", action="store", type="string",
                      help="Output file name.",dest="ofile")
    parser.add_option("-d", "--directory", action="store", type="string",
                      help="EPOS file to be read.",dest="directory")
    (options, args) = parser.parse_args()
    if len(sys.argv) == 1:
        parser.error('Specify an input file')
    return(options);

def smooth(data):
    smoothed = []
    normed = []
    for i in range(0, len(data)):
        if i < len(data)-2 :
            smoothed.append(np.mean(data[i:i+3]))
        elif i >= len(data)-2 :
            smoothed.append( np.mean(np.concatenate((data[i::], data[0:i - 1024 + 3]), axis=0)))

    for i in range(0, len(data)):
        normed.append((smoothed[i]) / (max(smoothed)) )    
    return normed

class InteractiveSelection:
    def __init__(self, figurename):
        self.press = None
        self.fig = figurename
        self.primitives = []

    def connect(self):
        self.cidpress = self.fig.canvas.mpl_connect('button_press_event', self.onclick)
        self.cidrelease = self.fig.canvas.mpl_connect('button_release_event', self.onrelease)

    def onclick(self, event):
        print 'button=%d, x=%d, y=%d, xdata=%f, ydata=%f'%(event.button, event.x, event.y, event.xdata, event.ydata)
        self.x0 = event.xdata
        self.y0 = event.ydata
        self.press = self.x0, self.y0, event.xdata, event.ydata

    def onrelease(self, event):
        print 'button_off=%d, x=%d, y=%d, xdata=%f, ydata=%f'%(event.button, event.x, event.y, event.xdata, event.ydata)
        self.x1 = event.xdata
        self.y1 = event.ydata

        self.center = self.x0 + (self.x1 - self.x0)/2
        self.amplitude = 1 * self.y0
        self.width =  abs(self.x1 - self.x0)
        self.primitives += [self.amplitude, self.center, self.width]
        print "self.primitives", self.primitives
        self.fig.canvas.draw()
    
    def gaussian_fitting(self, archive, off_pulse_left, off_pulse_right):
#        for i in range(0, len(self.primitives)/3):
#            if i ==0:
#                gg_init = models.Gaussian1D(self.primitives[i*3], self.primitives[i*3+1], self.primitives[i*3+2])
#            else:
#                gg_init = gg_init + models.Gaussian1D(self.primitives[i*3], self.primitives[i*3+1], self.primitives[i*3+2])
#        fitter = fitting.LevMarLSQFitter()
#        gg_fit = fitter(gg_init, np.linspace(0, len(archive), len(archive)), archive)
#        gaussian = gg_fit(np.linspace(0, len(archive), len(archive)))

        gaussian = utils.GaussianModel(archive, self.primitives)
        residuals = archive - gaussian
        utils.plot(archive, gaussian)
        
##########################
        border = utils.OnPulseBorders()
        border.ONPulseRegion(archive, off_pulse_left, off_pulse_right)
###########################



        
        self.KS_test = ks_2samp(archive[off_pulse_left:off_pulse_right], residuals[border.main_component[0]:border.main_component[1]])
        print self.KS_test[0], self.KS_test[1]

        while (self.KS_test[1] <= 0.05):
            if max(np.abs(residuals[border.main_component[0]:border.main_component[1]])) in residuals:
                self.primitives.append(archive[np.where(residuals == max(np.abs(residuals[border.main_component[0]:border.main_component[1]])) )].tolist()[0] )
                self.primitives.append(np.where(residuals == max(np.abs(residuals[border.main_component[0]:border.main_component[1]])) )[0] )
                self.primitives.append(10)
            else:
                self.primitives.append(archive[np.where(residuals == -max(np.abs(residuals[border.main_component[0]:border.main_component[1]])) )].tolist()[0])
                self.primitives.append(np.where(residuals == -max(np.abs(residuals[border.main_component[0]:border.main_component[1]])) )[0] )
                self.primitives.append(10)

            for i in range(0, len(self.primitives)/3):
                if i ==0:
                    gg = models.Gaussian1D(self.primitives[i*3], self.primitives[i*3+1], self.primitives[i*3+2])
                else:
                    gg = gg + models.Gaussian1D(self.primitives[i*3], self.primitives[i*3+1], self.primitives[i*3+2])
            fitter = fitting.LevMarLSQFitter()
            gg_fit_KS_test = fitter(gg, np.linspace(0, len(archive), len(archive)), archive)
            gaussian = gg_fit_KS_test(np.linspace(0, len(archive), len(archive)))
            residuals = archive - gaussian
            self.KS_test = ks_2samp(archive[off_pulse_left:off_pulse_right], residuals[border.main_component[0]:border.main_component[1]])
            print "KS_test", self.KS_test[0], self.KS_test[1]

        utils.plot(archive, gaussian)

class errorbars:
    def __init__(self, initial_values, archive, off_pulse_left, off_pulse_right):
        border = utils.OnPulseBorders()
        border.ONPulseRegion(archive, off_pulse_left, off_pulse_right)

        test = []
        test_theta = []
        for i in range(0, len(initial_values)/3):
            test = test + ["amplitude_ml_"+str(i)+",", "mean_ml_"+str(i)+",", "stddev_ml_"+str(i)+","]
            test_theta = test_theta + ["amplitude_"+str(i)+",", "mean_"+str(i)+",", "stddev_"+str(i)+","]
        test = test + ["sigma"]
        test_theta = test_theta + ["sigma"]
        print " ".join(str(x) for x in test_theta)
        theta_text = " ".join(str(x) for x in test_theta)
        test_ml = " ".join(str(x) for x in test)

        def lnlike(theta, archive):
            theta_text = theta
            for i in range(0, len(initial_values)/3):
                if i ==0:
                    gg_mcmc = models.Gaussian1D(theta[i*3], theta[i*3+1], theta[i*3+2])
                else:
                    gg_mcmc = gg_mcmc + models.Gaussian1D(theta[i*3], theta[i*3+1], theta[i*3+2])

            fitter = fitting.LevMarLSQFitter()
            gg_fit_mcmc = fitter(gg_mcmc, np.linspace(0, len(archive), len(archive)), archive)
            gaussian = gg_fit_mcmc(np.linspace(0, len(archive), len(archive)))
#            return -0.5*((np.sum(archive-gaussian)**2)/theta[-1]**2 )
            return np.exp(-0.5*((np.sum(archive-gaussian)**2)/theta[-1]**2  ))
#            return np.exp(-0.5*np.sum((archive-gaussian)**2*(1./theta[-1]**2)))

#######################
        nll = lambda *args: -lnlike(*args)
        result = op.minimize(nll, initial_values+[1], args=(archive))
        test_ml = result["x"]
        print "results", result["x"]

####comment#############        
#        for i in range(0, len(initial_values)/3):
#            if i ==0:
#                g = models.Gaussian1D(result["x"][i*3], result["x"][i*3+1], result["x"][i*3+2])
#            else:
#                g = g + models.Gaussian1D(result["x"][i*3], result["x"][i*3+1], result["x"][i*3+2])
#        fitter = fitting.LevMarLSQFitter()
#        g_fit = fitter(g, np.linspace(0, len(archive), len(archive)), archive)
#        gaussian = g_fit(np.linspace(0, len(archive), len(archive)) )

        gaussian = utils.GaussianModel(archive, result["x"])
#######################
#######################
        def lnprior(theta):
            levels = 10/100 #%
            theta_text = theta
            for i in range(0, len(initial_values)/3):
                if (0 < theta[-1] < 1) and (theta[i*3] - theta[i*3]*levels < theta[i*3] < theta[i*3] + theta[i*3]*levels) and (theta[i*3+1] - theta[i*3+1]*levels < theta[i*3+1] < theta[i*3+1] + theta[i*3+1]*levels) and (theta[i*3+2] - theta[i*3+2]*levels < theta[i*3+1] < theta[i*3+2] + theta[i*3+2]*levels):
                    return 0.0
            return -np.inf
        
        def lnprob(theta, archive):
            lp = lnprior(theta)
            if not np.isfinite(lp):
                return -np.inf
            return lp + lnlike(theta, archive)
    
        ndim, nwalkers = len(initial_values)+1, 500
        pos = [result["x"] + 1e-4*np.random.randn(ndim) for i in range(nwalkers)]
        sampler = emcee.EnsembleSampler(nwalkers, ndim, lnprob, args=(archive,))
        sampler.run_mcmc(pos, 5000)

        samples = sampler.chain[:, 500:, :].reshape((-1, ndim))
        print "SAMPLE SHAPE", samples.shape, len(samples)

        fig = corner.corner(samples, truths=initial_values+[1])
        fig.savefig('/homes/egraikou/pulsars/B1913+16/B1913+16_akaike.eps', dpi=100, bbox_inches='tight')

############################################
############################################
        gaussian_mcmc = []
        for theta_text in samples[np.random.randint(len(samples), size=100)]:
            for i in range(0, len(initial_values)/3):
                if i ==0:
                    gg_error = models.Gaussian1D(theta_text[i*3], theta_text[i*3+1], theta_text[i*3+2])
                else:
                    gg_error = gg_error + models.Gaussian1D(theta_text[i*3], theta_text[i*3+1], theta_text[i*3+2])
            fitter = fitting.LevMarLSQFitter()
            gg_error_fit = fitter(gg_error, np.linspace(0, len(archive), len(archive)), archive)
            gaussian_mcmc.append(gg_error_fit(np.linspace(0, len(archive), len(archive)) ))
        
        gaussian_1sigma_lower = []
        gaussian_1sigma_upper = []
        gaussian_2sigma_lower = []
        gaussian_2sigma_upper = []

        for k in range(0, len(archive)):
            gaussian_1sigma_lower.append(np.percentile([i[k] for i in gaussian_mcmc], 16))
            gaussian_1sigma_upper.append(np.percentile([i[k] for i in gaussian_mcmc], 84))
            gaussian_2sigma_lower.append(np.percentile([i[k] for i in gaussian_mcmc], 2.5))
            gaussian_2sigma_upper.append(np.percentile([i[k] for i in gaussian_mcmc], 97.5))


        print "local maximum 1 sigma lower limit", utils.local_max(gaussian_1sigma_lower, border.main_component)
        print "local maximum 1 sigma upper limit", utils.local_max(gaussian_1sigma_upper, border.main_component)
        print "local maximum 2 sigma lower limit", utils.local_max(gaussian_2sigma_lower, border.main_component)
        print "local maximum 2 sigma upper limit", utils.local_max(gaussian_2sigma_upper, border.main_component)
        print "local maximum", utils.local_max(gaussian, border.main_component)


        utils.width(gaussian, gaussian_1sigma_lower, gaussian_1sigma_upper, gaussian_2sigma_lower, gaussian_2sigma_upper, border.main_component)
        utils.separation(gaussian, gaussian_1sigma_lower, gaussian_1sigma_upper, gaussian_2sigma_lower, gaussian_2sigma_upper, border.main_component)
        utils.relative_height(gaussian, gaussian_1sigma_lower, gaussian_1sigma_upper, gaussian_2sigma_lower, gaussian_2sigma_upper, border.main_component)

        fig3 = plt.figure()
        plt.plot(np.linspace(0, len(archive), len(archive)), archive)
        plt.plot(np.linspace(0, len(archive), len(archive)), gaussian, 'r-')
        plt.plot(np.linspace(0, len(archive), len(archive)), gaussian_1sigma_upper, 'r--')
        plt.plot(np.linspace(0, len(archive), len(archive)), gaussian_1sigma_lower, 'r--')
        plt.plot(np.linspace(0, len(archive), len(archive)), gaussian_2sigma_upper, 'b:')
        plt.plot(np.linspace(0, len(archive), len(archive)), gaussian_2sigma_lower, 'b:')
#        plt.axhline(max(gaussian)*95/100, color = 'k', linestyle='-', linewidth = 3 )
#        plt.axhline(max(gaussian_2sigma_lower)*95/100, color = 'k', linestyle='-', linewidth = 3 )
        plt.show()

################################
################################
def main():
    options = CommandLineParse()

    arch = psrchive.Archive_load(str(options.directory)+str(options.infile))
    arch.fscrunch()
    arch.tscrunch()
    arch.pscrunch()
    arch.remove_baseline()
    arch.centre()
    data = arch.get_data()
    #########
#    initial_archive = np.array(smooth(data[0, 0, 0, :]))
    initial_archive = data[0, 0, 0, :]

    border = utils.OnPulseBorders()
    border.OFFPulseRegion(initial_archive)
#    border.ONPulseRegion(initial_archive)

    fig1 = plt.figure(1)
    plt.plot(np.linspace(0, len(data[0, 0, 0, :]), len(data[0, 0, 0, :])), initial_archive)
    plt.xlim(0, 1024)

    dr = InteractiveSelection(fig1)
    dr.connect()
    plt.show()

    dr.gaussian_fitting(initial_archive, border.off_pulse_left, border.off_pulse_right)
    error = errorbars(dr.primitives, initial_archive, border.off_pulse_left, border.off_pulse_right)

if __name__ == "__main__":
    main()
