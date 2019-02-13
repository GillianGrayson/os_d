from config.run import *
from copy import deepcopy
from plot_data.generator import *
import numpy as np
from infrastructure.path import *
import math
from itertools import combinations
import statsmodels.api as sm


class FitPlane(Generator):

    def generate_plot_data(self, space):

        configs = space.configs
        xs = space.xs
        ys = space.ys

        bin_begin = 1e-10
        num_decades = 15
        bin_end = bin_begin * pow(10, num_decades)
        num_bin_per_decade = 10
        num_bins = num_bin_per_decade * num_decades
        bin_borders = np.zeros(num_bins + 1)
        bin_centers = np.zeros(num_bins)
        for bin_id in range(0, num_bins + 1):
            bin_borders[bin_id] = bin_begin * pow(10, ((bin_id - 1) / num_bin_per_decade))
            if bin_id < num_bins:
                bin_centers[bin_id] = bin_begin * pow(10, ((bin_id - 1 + 0.5) / num_bin_per_decade))
        bin_diff = np.diff(bin_borders)

        non_inc_count = 0

        data = []
        config_id = 0
        for x in xs:
            x_data = []
            for y in ys:
                config = configs[config_id]
                path = get_data_path(config)

                pdf = np.zeros(num_bins)

                num_trajectories = config.auxiliary['num_trajectories']
                for traj_id in range(0, num_trajectories):
                    fn = path + '/' + 'jump_times_' + str(traj_id) + '_' + get_file_suffix(config) + '.txt'
                    jump_times = np.loadtxt(fn)
                    diffs = np.diff(jump_times)
                    for d in diffs:
                        if d >= bin_begin and d <= bin_end:
                            bin_id = np.floor((np.log10(d) - np.log10(bin_begin)) * num_bins / (np.log10(bin_end) - np.log10(bin_begin) + 1e-8))
                            pdf[bin_id] += 1
                        else:
                            non_inc_count += 1

                norm = np.sum(pdf)
                for bin_id in range(0, num_bins):
                    pdf[bin_id] /= (norm * bin_diff[bin_id])

                norm_check = 0.0
                for bin_id in range(0, num_bins):
                    norm_check += pdf[bin_id] * bin_diff[bin_id]


def power_law_regression(x, y, x_min, x_max):
    ids = [id for id in range(0, len(x)) if y[id] > 0 and x_min < x[id] and x[id] < x_max]
    x = list(np.array(x)[ids])
    y = list(np.array(y)[ids])

    count = len(x)

    log_x = np.log(x)
    log_y = np.log(y)

    exog = sm.add_constant(log_x)
    results = sm.OLS(log_y, exog).fit()

    R2 = results.rsquared
    intercept = results.params[0]
    slope = results.params[1]

    log_y_fit = [(curr_x * slope + intercept) for curr_x in log_x]
    y_fit = [pow(10.0, curr_y) for curr_y in log_y_fit]

    return [slope, R2, count, y_fit]


def find_power_law_range(x, y):
    i_max = 0
    j_max = 0
    ids = [id for id in range(0, len(y)) if y[id] > 0]

    x = list(np.array(x)[ids])
    y = list(np.array(y)[ids])

    target_len = 1.0
    R2_max = 0.0

    for i in range(0, len(x)):
        for j in range(i, len(x)):
            [slope, R2, count, y_fit] = power_law_regression(x, y, x[i], x[j])

            if count < 10:
                continue

            curr_len = np.log10(x[j]) - np.log10(x[i])

            if R2 > 0.98 and curr_len > 1.0:
                curr_len += 120.0 * (R2 - 0.98)
                if abs(target_len - curr_len) <= 0.01:
                    if R2 > R2_max:
                        target_len = curr_len
                        R2_max = R2
                        i_max = x[i]
                        j_max = x[j]

                elif target_len + 1 < curr_len:
                    target_len = curr_len
                    R2_max = R2
                    i_max = x[i]
                    j_max = x[j]

    return [i_max, j_max]