"""Cooling temerature scheduler for multiplicative simulated annealing"""

import logging
from typing import Dict

import numpy as np

from importlib import reload
from runtime.tcr_pmhc_experiments.peptide_generation_with_sim_anneal.temperature import scheduler
reload(scheduler)
TemperatureScheduler = scheduler.TemperatureScheduler

class MultiplicativeScheduler(TemperatureScheduler):
    
    def __init__(
            self, T_0: float, T_f: float, a: float,
            metric_history: Dict, region: str, **kwargs):
        #super().__init__()
        super().__init__(T_0, T_f, metric_history, region)
        self.a = a
        self.schedule = self.create_schedule(self.T_0, self.T_f, self.a, self.iterations)

    def create_schedule(
            self, T_0: float, T_f: float, a: float,
            iterations: int):
        """Multiplicative temperature schedule"""
        return (T_0 - T_f) * np.power(a, np.arange(iterations)) + T_f
    
    def get_next_temperature(self):
        return self.schedule[self.k]





