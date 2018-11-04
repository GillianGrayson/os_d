from config.run import *
from copy import deepcopy
from plot_data.generator import *
import numpy as np
from infrastructure.path import *
import math
from itertools import combinations


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

    xx = x
    y0 = y
    ids = [id for id in range(0, len(x)) if y[id] > 0 and x_min < x[id] and x[id] < x_max]
    x = list(np.array(x)[ids])
    y = list(np.array(y)[ids])

    cnt = len(x)

    logx = np.log(x)
    logy = np.log(y)