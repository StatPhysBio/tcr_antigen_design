"""Cooling temerature scheduler for thermodynamic simulated annealing"""

import logging
from typing import Dict

import numpy as np

#from protein_holography.anneal.temperature.scheduler import TemperatureScheduler

from importlib import reload
from runtime.tcr_pmhc_experiments.peptide_generation_with_sim_anneal.temperature import scheduler
reload(scheduler)
TemperatureScheduler = scheduler.TemperatureScheduler

class ThermodynamicScheduler(TemperatureScheduler):
    def __init__(
            self, T_0: float, T_f: float, k_A: float,
            region: str, metric_history: Dict, **kwargs):
        super().__init__(T_0, T_f, metric_history, region)
        logging.info("Updated thermodynamic scheduler")
        self.region = region
        self.k_A = k_A
        self.tot_inc_S = 0.
        self.tot_acc_delta_E = 0.

    def get_next_temperature(self):
        logging.info("New scheme")
        """Get the TSA-prescribed temperature"""
        opt_dir = -1

        logging.debug(f"Curr step {self.k}")
        curr_delta_E = (
            self.metric_history[f"{self.region} pnE"][self.k] -
            self.metric_history[f"{self.region} pnE"][self.k-1])
        logging.debug(f"Curr delta E = {curr_delta_E}")
        # add accepted delta E
        if self.metric_history["acceptance"][self.k]:
            self.tot_acc_delta_E += curr_delta_E
        logging.debug(f"Curr accepted delta E = {self.tot_acc_delta_E}")
        if opt_dir * curr_delta_E > 0:
            self.tot_inc_S += curr_delta_E/self.metric_history[f"{self.region} temperature"][self.k-1]

        # return statements

        if opt_dir * self.tot_acc_delta_E >= 0.:
            logging.info("Setting T = T_0 because accepted delta E < 0")
            return self.T_0
        if self.tot_inc_S == 0.:
            logging.info("Setting T = T_0 because deltaS+ = 0")
            return self.T_0
        new_T = -self.k_A * (opt_dir * self.tot_acc_delta_E / # minus sign for opposite energy
                            (opt_dir * self.tot_inc_S))
        
        logging.info(f"New T = {new_T}")
        logging.debug(f"Temps = {self.metric_history[f'{self.region} temperature']}")
        return new_T
    
    def old_get_next_temperature(self):
        """Get the TSA-prescribed temperature"""
        # get associated histories
        energies = self.metric_history[f"{self.region} pnE"][:self.k+2]
        temps = self.metric_history[f"{self.region} temperature"][:self.k+1]
        acceps = self.metric_history["acceptance"][1:self.k+2]

        #calculate all conditional or terms in fraction
        delta_E = np.diff(energies)
        logging.info(f"Current step = {self.k}")
        logging.info(f"temps = {temps}")
        acc_delta_E = np.sum(delta_E * acceps)
        if acc_delta_E <= 0:
            return temps[0]
        logging.debug(f"Accepted delta E = {acc_delta_E}")
        proposed_delta_S = -delta_E/temps
        logging.debug(f"propsoed delta S= {proposed_delta_S}")
        inc_S_steps = proposed_delta_S > 0
        print(inc_S_steps)
        tot_inc_S = -np.sum(proposed_delta_S[inc_S_steps])
        if tot_inc_S == 0:
            return temps[0]
        logging.debug(f"Total increase in S = {tot_inc_S}")
        new_T =- self.k_A * (acc_delta_E / tot_inc_S)
        logging.debug(f"new T = {new_T}")
        if new_T == 0.:
            new_T = temps[0]
        return new_T
        
        

    