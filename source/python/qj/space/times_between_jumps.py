from config.run import *
from copy import deepcopy
from space.space import *
import numpy as np


class TBJ(Space):

    def xy_space(self, config):

        xs = []
        ys = []
        configs = []

        if config.run.quantum_system is QuantumSystem.photonic:
            x_start = 0.05
            x_shift = 0.05
            x_num = 100
            xs = np.linspace(x_start, x_start + (x_num - 1) * x_shift, x_num)

            y_start = 0.05
            y_shift = 0.05
            y_num = 100
            ys = np.linspace(y_start, y_start + (y_num - 1) * y_shift, y_num)

            for x in xs:
                for y in ys:
                    curr_config = deepcopy(config)
                    curr_config.params.driving['jcs_drv_ampl'] = x
                    curr_config.params.driving['jcs_drv_part_1'] = 0.98 * y
                    curr_config.params.driving['jcs_drv_part_2'] = 1.00 * y

                    configs.append(curr_config)

        elif config.run.quantum_system is QuantumSystem.dimer:
            x_start = 0.01
            x_shift = 0.01
            x_num = 1
            xs = np.linspace(x_start, x_start + (x_num - 1) * x_shift, x_num)

            y_start = 0.01
            y_shift = 0.01
            y_num = 1
            ys = np.linspace(y_start, y_start + (y_num - 1) * y_shift, y_num)

            for x in xs:
                for y in ys:
                    curr_config = deepcopy(config)
                    curr_config.params.driving['dimer_drv_ampl'] = x
                    curr_config.params.driving['dimer_prm_U'] = y

                    configs.append(curr_config)

        self.xs = xs
        self.ys = ys
        self.configs = configs


    def x_space(self, config):
        pass


    def y_space(self, config):
        pass