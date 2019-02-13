from abc import ABCMeta, abstractmethod


class Space(metaclass=ABCMeta):

    def __init__(self):
        self.xs = []
        self.ys = []
        self.configs = []

    @abstractmethod
    def x_space(self, config):
        pass

    @abstractmethod
    def y_space(self, config):
        pass

    @abstractmethod
    def xy_space(self, config):
        pass
