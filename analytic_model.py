# class to encapsulate analytic functions which were derived in 
# Blok, NK et al., Faraday Disc. (2015) XXXCITATION

import simulation; reload(simulation)
import numpy as np

class Faraday(simulation.Simulation):

    def decay_formula(self): #XXX NAME
        """
        Calculates the predicted coherence decay.

        in: void (takes input parameters from supplied parameters)
        rtype: float/np.array; nuclear spin coherence
        """
        c = self.get_param('coupling')
        r = self.get_param('average_repump_time')
        N = self.get_param('entangling_attempts')
        p_nv_in_1 = self.get_param('pflip')
        return (1 - p_nv_in_1 + p_nv_in_1 * np.exp(-(2 * np.pi * c * r)**2 / 2))**N

    def faraday_decay_constant(self):
        """
        Returns analytic result when solving self.decay_formula() = 1/e 
        for the number of entangling attempts N. 
        
        in: void (takes input parameters from supplied parameters)
        out: N_{1/e}. rtype: float/np.array; the number of attempts until the nuclear spin dephases
        """

        c = self.get_param('coupling')
        r = self.get_param('average_repump_time')
        p_nv_in_1 = self.get_param('pflip')

        return -1/(np.log(1 - p_nv_in_1 + p_nv_in_1 * np.exp(-0.5 * (2 * np.pi * c * r)**2)))
    
    def get_carbon_state_fidelity(self):
        """
        in: void (takes input from supplied parameters)
        out: void (creates self.cspin_fidelity attribute for later use) 
        """
        self.cspin_fidelity = self.decay_formula() / 2. + 0.5