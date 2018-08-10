# Monte-Carlo simulation class.

import numpy as np
from simulation import Simulation


class QMC(Simulation):

    def __init__(self, params = None):

        Simulation.__init__(self, params)

    def compute_fidelity_from_scratch(self):
        """
        wrapper function to run a full cycle of MC simulation with the current parameter setting
        creates the attribute self.cspin_fidelity as function of entangling attempts
        
        in: void
        out: void
        """
        self.get_nv_state_during_sequence()
        self.calc_nuclear_phase_distribution()
        self.get_carbon_state_fidelity()
        return

    def get_nv_state_during_sequence(self):
        """
        Creates boolean flags which determine the NV spin state during each entangling attempt
        
        in: void
        out: void
        """
        self._generate_random_numbers()

        init_infidelity = self.get_param('init_infidelity') 
        mw_infidelity   = self.get_param('mw_infidelity') 
        sin2_alpha      = self.get_param('pflip') 

        # initial states for the NV ground state
        self._nv_in_0 = self._rns_nv_init < (1. - init_infidelity)
        self._nv_in_p1 = self._rns_nv_init >= (1. - init_infidelity / 2.)
        self._nv_in_m1 = ~(self._nv_in_0 | self._nv_in_p1)

        # mw pi-pulse flip errors
        mw_failed = self._rns_mw > 1. - mw_infidelity

        # state projection after first mw pulse
        nv_alpha_flip = self._rns_nv_alpha > sin2_alpha

        # state was in ms = -1 after initial microwave pulse and remains there after the second pulse
        mw_can_fail = mw_failed & (~self._nv_in_p1)
        self._mw_m1 = mw_can_fail & ((self._nv_in_m1 & (~nv_alpha_flip)) | (nv_alpha_flip & self._nv_in_0))
        
        # state was in ms = 0 after initial pulse and remains there
        self._mw_0 = mw_can_fail & ((self._nv_in_0 & (~nv_alpha_flip)) | (nv_alpha_flip & self._nv_in_m1))

        # if _nv_mw_repump: the nv is pumped to ms = 0. Need to draw a random repump time
        self._nv_repump = self._mw_m1 | (self._nv_in_0 & (~ nv_alpha_flip)) | self._nv_in_p1 | (self._nv_in_m1 & nv_alpha_flip)

    def calc_nuclear_phase_distribution(self):
        """
        estimates the acquired nuclear spin phase per entangling attempt.

        in: void
        out: void
        """

        decoupling_duration = self.get_param('larmor_period') * self.get_param('larmor_order')
        T = self.get_param('T') 
        coupling_strength = self.get_param('coupling')


        # calculate phases applied to carbon spin state for all electron states
        init_phase_e_inp1 = coupling_strength * T
        init_phase_e_inm1 = -init_phase_e_inp1

        # choose a rotating frame with the average frequency: (omega_0 + omega_m1)/2
        mw_phase_e_in0 = -decoupling_duration*coupling_strength 
        mw_phase_e_inm1 = -mw_phase_e_in0
        mw_phase_e_inp1 = 1.5 * 2 * decoupling_duration*coupling_strength    

        # repumping phases and timing mismatches:
        repump_jitter = -1 * coupling_strength * self._get_static_repump_jitter() 
        repump_phases = -1 * self._repump_phase()

        ## add up all nuclear spin phases for each entangling attempt
        carbon_phases  = np.zeros(np.shape(self._nv_in_m1))
        carbon_phases += self._nv_in_m1 * init_phase_e_inm1 
        carbon_phases += self._nv_in_p1 * (init_phase_e_inp1 + mw_phase_e_inp1)
        carbon_phases += self._mw_0 * mw_phase_e_in0 
        carbon_phases += self._mw_m1 * mw_phase_e_inm1 
        carbon_phases += self._nv_repump * (repump_phases + repump_jitter)

        self.carbon_phases = carbon_phases * 2 * np.pi 

    def get_carbon_state_fidelity(self):
        """
        calculates the nuclear spin state fidelities according to:
        F = (1 + cos(phase))/2

        in: void
        out: void (creates instance attribute cspin_fidelity for further use instead)
        """

        # access a certain nuclear phase with self.carbon_phases[entangling_attempt][MC_repetition]
        # cumsum calculates the phase that is built up from all entangling attempts
        # average over all MC repetitions via np.average(matrix, axis = 1)
        self.cspin_fidelity = np.cumsum(self.carbon_phases, axis = 0)
        self.cspin_fidelity = np.cos(self.cspin_fidelity)
        self.cspin_fidelity = (np.average(self.cspin_fidelity, axis = 1) + 1.)/2.

    def find_attempts_from_fidelity(self, target):
        """
        finds the number of attempts that returns a fidelity closest to the target fidelity
        commonly used to find N_{1/e}.

        in: 0.0 < target < 1.0: desired fidelity value
        out: int, number of entangling attempts closest to target.

        Note: No need for binary search as data is typically small. 
        """
        return (np.abs(self.cspin_fidelity - target)).argmin()

    def _generate_random_numbers(self):
        """
        Generates and stores matrices of random numbers for later use in MC sim.
        
        in: void
        out: void
        """
        nrows = self.get_param('entangling_attempts')
        ncols = self.get_param('repetitions')
        self._rns_mw            = np.random.rand(nrows, ncols)
        self._rns_nv_alpha      = np.random.rand(nrows, ncols)
        self._rns_nv_init       = np.random.rand(nrows, ncols)
        self._rns_nv_repump     = np.random.rand(nrows, ncols)

    def _repump_phase(self):
        """
        converts uniformly distributed random numbers for NV repumping to the target exponential distribution:
        1-e^(- t / tau) = RN

        from this, obtain the repump time t: t = -ln(1 - RN) * tau
        tau is the avg. repump duration of the experiment

        in: void
        out: 2D np.array containing exponentially distributed repumping times
        """

        coupling = self.get_param('coupling')
        repump_duration = self.get_param('average_repump_time')
        return (-np.log(1 - self._rns_nv_repump) - 1) * coupling * (repump_duration) 

    def _get_static_repump_jitter(self):
        """
        Creates offsets in the average repumping time that are static for each experimental repetition

        in: void
        out: 1D np.array
        """
        repump_offset   = self.get_param('repump_time_offset')
        repump_sigma    = self.get_param('repump_time_jitter')
        repetitions     = self.get_param('repetitions')

        if repump_sigma <= 0.0:
            return np.ones(repetitions)*repump_offset
        else:
            # assume that we only observe an increase in repump offsets as laser intensity decreases
            return np.abs(np.random.normal(loc=repump_offset,
                                           scale=repump_sigma,
                                           size=length)) 
