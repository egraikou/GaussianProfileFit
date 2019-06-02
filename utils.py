"""
This code contains utility code that is used in gaussian fitting code.

"""

import numpy as np
import time
import math
import psrchive
import matplotlib.pyplot as plt
from astropy.modeling import models, fitting

def plot(archive, gaussian):
    fig1 = plt.figure()
    plt.plot(np.linspace(0, len(archive), len(archive)), archive)
    plt.plot(np.linspace(0, len(archive), len(archive)), gaussian)
    plt.plot(np.linspace(0, len(archive), len(archive)), archive - gaussian )
#    plt.axvline(self.component_border[0], color = 'r', linestyle='-', linewidth = 3)
#    plt.axvline(self.component_border[1], color = 'r', linestyle='-', linewidth = 3)
    plt.xlim(0, len(archive))
    plt.show()


class OnPulseBorders:
    def OFFPulseRegion(self, archive):
        fig2 = plt.figure()
        plt.title('Pick the OFF pulse region, check the command line')
        plt.plot(np.linspace(0, len(archive), len(archive)), archive)
        plt.xlim(0, len(archive))
        plt.pause(0.01)
        self.off_pulse_left = input("Enter the beggining of the off pulse region: ")
        self.off_pulse_right = input("Enter the end of the off pulse region: ")
        plt.close(fig2)

    def ONPulseRegion(self, archive, off_pulse_left, off_pulse_right):
        self.component_border = []
        self.main_component = []
        self.interpulse = []
        rms_off = np.std(archive[off_pulse_left:off_pulse_right])    
        for i in range(np.array(archive).tolist().index(max(archive)), 3, -3):
            if sum(archive[i-3:i])*3 < 3 * rms_off * 3:
                left_border = i
                self.main_component += [left_border ]
                break

        for i in range(np.array(archive).tolist().index(max(archive)), len(archive)-3, 3):
            if sum(archive[i:i+3])*3 < 3 * rms_off * 3:
                right_border = i
                self.main_component += [right_border ]
                break

        print self.main_component
#        print int(math.log10(np.std(archive[0:self.main_component[0]]))), int(math.log10(rms_off))

        if int(math.log10(np.std(archive[0:self.main_component[0]]))) > int(math.log10(rms_off)):
            print "the interpulse is between 0, main component", int(math.log10(np.std(archive[0:self.main_component[0]]))), int(math.log10(rms_off))
            print "the interpulse is between 0, main component", np.std(archive[0:self.main_component[0]]), rms_off
            for i in range(np.array(archive).tolist().index(max(archive[0:self.main_component[0]])), 3, -3  ):
                if sum(archive[i-3:i])*3 < 3 * rms_off *3:
                    left_border = i
                    self.interpulse += [left_border]
                    break
            for i in range(np.array(archive).tolist().index(max(archive[0:self.main_component[0]])), len(archive)-3, 3  ):
                if sum(archive[i:i+3])*3 < 3 * rms_off * 3:
                    right_border = i
                    self.interpulse += [right_border ]
                    break
   
        elif int(math.log10(np.std(archive[0:self.main_component[0]]))) < int(math.log10(rms_off)):
            print "the interpulse is between main component and the end of the archive", int(math.log10(np.std(archive[0:self.main_component[0]]))), int(math.log10(rms_off))
            print "the interpulse is between main component and the end of the archive", np.std(archive[0:self.main_component[0]]), rms_off
            for i in range(np.array(archive).tolist().index(max(archive[self.main_component[-1]:-1] )), 3, -3  ):
                if sum(archive[i-3:i])*3 < 3 * rms_off *3:
                    left_border = i
                    self.interpulse += [left_border]
                    break

            for i in range(np.array(archive).tolist().index(max(archive[self.main_component[-1]:-1])), len(archive)-3, 3  ):
                if sum(archive[i:i+3])*3 < 3 * rms_off * 3:
                    right_border = i 
                    self.interpulse += [right_border ]
                    break
        print self.interpulse

        return self.main_component


def GaussianModel(data, GaussianParameters):
    """Gaussian parameters is an array that contains the parameters od each Gaussian [amplitude1, mean1, wdth1, amplitude2, mean2, width2, ...] """
    gg = []
    for i in range(0, len(GaussianParameters)/3):
        if i==0:
            gg = models.Gaussian1D(GaussianParameters[i*3], GaussianParameters[i*3+1], GaussianParameters[i*3+2])
        else:
            gg = gg + models.Gaussian1D(GaussianParameters[i*3], GaussianParameters[i*3+1], GaussianParameters[i*3+2])

    fitter = fitting.LevMarLSQFitter()
    gg_fit = fitter(gg, np.linspace(0, len(data), len(data)), data)
    gaussian = gg_fit(np.linspace(0, len(data), len(data)))
    return gaussian

def local_max(data, OnPulseBorders):
#    print "inside local_max", np.array(OnPulseBorders[0]), np.array(OnPulseBorders[1])
    local_max = []
    local_max_index = []
    for i in range(1, len(data)-1):
        if (data[i-1] < data[i]) and (data[i+1] < data[i]) and (OnPulseBorders[0] < i) and (i < OnPulseBorders[-1]):
            local_max.append(data[i])
            local_max_index.append(i)
    return local_max_index


def find_nearest(archive, value):
    idx = (np.abs(np.array(archive)-value)).argmin()
    return archive[idx], idx


def relative_height(gaussian, gaussian_1sigma_lower, gaussian_1sigma_upper, gaussian_2sigma_lower, gaussian_2sigma_upper, MainPulseEdges):
    leading_component = gaussian[local_max(gaussian, MainPulseEdges)[0]]
    leading_component_1sigma = (gaussian_1sigma_lower[local_max(gaussian, MainPulseEdges)[0]] - gaussian_1sigma_upper[local_max(gaussian, MainPulseEdges)[0]])/2
    leading_component_2sigma = (gaussian_2sigma_lower[local_max(gaussian, MainPulseEdges)[0]] - gaussian_2sigma_upper[local_max(gaussian, MainPulseEdges)[0]])/2
    trailing_component = gaussian[local_max(gaussian, MainPulseEdges)[-1]]
    trailing_component_1sigma = (gaussian_1sigma_lower[local_max(gaussian, MainPulseEdges)[-1]] - gaussian_1sigma_upper[local_max(gaussian, MainPulseEdges)[-1]])/2
    trailing_component_2sigma = (gaussian_2sigma_lower[local_max(gaussian, MainPulseEdges)[-1]] - gaussian_2sigma_upper[local_max(gaussian, MainPulseEdges)[-1]])/2
    print "relative height", trailing_component - leading_component, "+/-", math.sqrt(leading_component_1sigma**2+ trailing_component_1sigma**2), "(1sigma)", "+/-", math.sqrt(leading_component_2sigma**2+ trailing_component_2sigma**2), "(2sigma)"

def separation(gaussian, gaussian_1sigma_lower, gaussian_1sigma_upper, gaussian_2sigma_lower, gaussian_2sigma_upper, MainPulseEdges):
    leading_component = local_max(gaussian, MainPulseEdges)[0]
    trailing_component = local_max(gaussian, MainPulseEdges)[-1]
    leading_component_1sigma = (local_max(gaussian_1sigma_upper, MainPulseEdges)[0] - local_max(gaussian_1sigma_lower, MainPulseEdges)[0])/2
    leading_component_2sigma = (local_max(gaussian_2sigma_upper, MainPulseEdges)[0] - local_max(gaussian_2sigma_lower, MainPulseEdges)[0])/2

    trailing_component_1sigma = (local_max(gaussian_1sigma_upper, MainPulseEdges)[-1] - local_max(gaussian_1sigma_lower, MainPulseEdges)[-1])/2
    trailing_component_2sigma = (local_max(gaussian_2sigma_upper, MainPulseEdges)[-1] - local_max(gaussian_2sigma_lower, MainPulseEdges)[-1])/2

    print "separation", trailing_component - leading_component, "degrees", "+/-", math.sqrt(leading_component_1sigma**2+ trailing_component_1sigma**2), "degrees", "(1sigma)", "+/-", math.sqrt(leading_component_2sigma**2+ trailing_component_2sigma**2), "degrees", "(2sigma)"


def width(gaussian, gaussian_1sigma_lower, gaussian_1sigma_upper, gaussian_2sigma_lower, gaussian_2sigma_upper, MainPulseEdges):
    levels = [5, 10, 20, 30, 40, 50, 60, 70, 80, 90, 95]
    for level in levels:
        position = []
        for local_pick_position in local_max(gaussian, MainPulseEdges):
            if max(gaussian)*level/100 < gaussian[local_pick_position]:
                position.append(local_pick_position)

        width_right = position[-1] + find_nearest(gaussian[position[-1]:-1], max(gaussian)*level/100)[1]
        width_left = find_nearest(gaussian[0:position[0]], max(gaussian)*level/100)[1]

        width_right_1sigma = (find_nearest(gaussian_1sigma_lower[position[-1]:-1], max(gaussian_1sigma_lower)*level/100)[1] - find_nearest(gaussian_1sigma_upper[position[-1]:-1], max(gaussian_1sigma_upper)*level/100)[1])/2
        width_left_1sigma = (find_nearest(gaussian_1sigma_lower[0:position[0]], max(gaussian_1sigma_lower)*level/100)[1] - find_nearest(gaussian_1sigma_upper[0:position[0]], max(gaussian_1sigma_upper)*level/100)[1])/2
        width_right_2sigma = (find_nearest(gaussian_2sigma_lower[position[-1]:-1], max(gaussian_2sigma_lower)*level/100)[1] - find_nearest(gaussian_2sigma_upper[position[-1]:-1], max(gaussian_2sigma_upper)*level/100)[1])/2
        width_left_2sigma = (find_nearest(gaussian_2sigma_lower[0:position[0]], max(gaussian_2sigma_lower)*level/100)[1] - find_nearest(gaussian_2sigma_upper[0:position[0]], max(gaussian_2sigma_upper)*level/100)[1])/2
        
        print "width", level,'%:', width_right - width_left, "degrees", "+/-", math.sqrt(width_right_1sigma**2+width_left_1sigma**2), "degrees", "(1sigma)", "+/-", math.sqrt(width_right_2sigma**2+width_left_2sigma**2), "degrees", "(2 sigma)"

