from enum import Enum


class QuantumSystem(Enum):
    dimer = 0
    photonic = 1


class Process(Enum):
    lyapunov_std = 0
    evolution_std = 1
    corr_dim_std = 2
    sigma_std = 3
    evolution_deep = 4
    lyapunov_deep = 5
    lyapunov_all = 6


class Propagation(Enum):
    quantum_jumps = 0
    runge_kutta_4 = 1


class Run:
    def __init__(self, quantum_system, process, propagation):
        self.quantum_system = quantum_system
        self.process = process
        self.propagation = propagation