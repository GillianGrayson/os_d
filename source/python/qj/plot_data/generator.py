from abc import ABCMeta, abstractmethod


class Generator(metaclass=ABCMeta):

    @abstractmethod
    def generate_plot_data(self, space):
        pass