"""Module for temperature scheduling"""
from abc import ABCMeta, abstractmethod
from argparse import Namespace
import logging
from typing import Dict, List

import numpy as np

class TemperatureScheduler(metaclass=ABCMeta):
    def __init__(
            self, T_0: float, T_f: float,
            metric_history: Dict, region: str):
        self.metric_history = metric_history
        self.T_0 = T_0
        self.T_f = T_f
        self.iterations = np.shape(metric_history["acceptance"])[0]
        self.k = 1
        self.region = region
        self.metric_history[f"{self.region} temperature"][0] = T_0
        self.metric_history[f"{self.region} temperature"][1] = T_0
        
    @abstractmethod
    def get_next_temperature(self):
        logging.error("TemperatureScheduler class not intended to be run. "\
                      "This class just serves as a super template for "\
                      "other temperature schedulers")
        pass

    def get_temperature(self):
        return self.metric_history[f"{self.region} temperature"][self.k]

    def step(self):
        logging.debug(f"Temperature_scheduler step {self.k}")
        if self.k >= self.iterations-1:
            logging.info("Scheduler step ignored at end of annealing")
        else:
            self.metric_history[f"{self.region} temperature"][self.k+1] = self.get_next_temperature()
        self.k += 1
        
    
def get_schedule(
        args: Namespace):
    """Get temperature schedule from possiblities"""
    if args.schedule == "multiplicative":
        schedule = multiplicative(
            args.T0, args.Tf, args.a, args.iters)
    try:
        return schedule
    except:
        logging.error('No valid schedule type specified')
        return

def get_TSA_temperature(
        energies: List,
        temps: List,
        acceptance: List,
        k_A: float) -> float:
    """Get the thermodynamic-SA-prescribed temperature"""
    delta_E = np.diff(energies)
    proposed_delta_S =  delta_E/temps
    proposed_inc_S = np.argwhere(proposed_delta_S < 0)
    tot_inc_S = np.sum(proposed_delta_S[proposed_inc_S])
    acc_delta_E = np.sum(delta_E*acceptance[1:])
    logging.info(f"Proposed delta S: {tot_inc_S}")
    logging.info(f"Accepted delta E: {acc_delta_E}")
    if acc_delta_E <= 0:
        return temps[0]
    if tot_inc_S == 0:
        return temps[0]
    return -k_A * (acc_delta_E/tot_inc_S)


def multiplicative(
    T_0: float,
    T_f: float,
    a: float,
    iterations: int
):
    """Multiplicative temperature schedule"""
    return (T_0 - T_f) * np.power(a,np.arange(iterations)) + T_f