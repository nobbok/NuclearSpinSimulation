# abstract base class for analytic/MC simulation
# provides a parameter interface for inheriting classes

import abc

class Simulation(object):
    __metaclass__ = abc.ABCMeta

    def __init__(self, params = None):
        """
        loads standard settings from parameters.py

        kwargs:
        type params: dict
        """

        #load standard settings
        self._params = self._load_params() #load default parameters from parameters.py
        self._create_parameter_map()

        self.cspin_fidelity = None

        ### loop through custom parameter dictionary to overwrite defaults from parameters.py
        if params:
            for k in params:
                self.set_param(k , params[k])

    def _load_params(self):
        import parameters; reload(parameters)
        return parameters.params

    def get_param(self, parameter_name):
        route = self._parameter_map[parameter_name]
        return self._get_param_from_map(*route)

    def get_params(self):
        return self._params

    def set_param(self, parameter_name, val):
        self._set_param_from_map(self._parameter_map[parameter_name], val)

    def _get_param_from_map(self, *route):
        if len(route) == 1:
            return self._params[route[0]]
        return self._get_param_from_map(*route[:-1])[route[-1]]

    def _set_param_from_map(self, route, val):
        d = self._params
        for k in route[:-1]:
            d = d[k]
        d[route[-1]] = val

    def _create_parameter_map(self):
        """
        updates the parameter dict with a map to all keys in parameters
        This mpa is used to set/get parameters
        Note that naming collisions in the parameters are not allowed currently. 
        TODO: add JSON support to fix naming collisions (e.g. coherence time for NV and nuclear spin)
        """
        def helper(k, curKey, curParams):
            if k and type(curParams) != dict: # found an actual simulation parameter
                self._parameter_map[k] = [key for key in curKey]
                return

            else:
                for key in curParams:
                    curKey.append(key)
                    helper(key, curKey, curParams[key])
                    curKey.pop()

        self._parameter_map = {}
        helper('', [], self._params)

    @abc.abstractmethod
    def get_carbon_state_fidelity(self):
        """
        required function for all inheriting classes
        """
        return