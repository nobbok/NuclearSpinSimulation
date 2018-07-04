params = {}

params['carbon_params'] = {
    'coupling'  : 80e3, # in Hz
    'T2'        : 60e-3, # in seconds // not used yet
}

params['nv_params'] = {
    'average_repump_time'   : 220e-9, # in seconds 
    'repump_time_jitter'    : 0.0,    # in seconds
    'repump_time_offset'    : 0.0,    # in seconds 
    'pflip'                 : 0.5,    # probability
    'mw_infidelity'         : 0.008,  # probability
    'init_infidelity'       : 0.001,  # probability
}

params['simulation_params'] = {
    'entangling_attempts'   : 2000,
    'repetitions'           : 1000,
    'larmor_period'         : 2.256e-6, # in s; @ magnetic Field B = 414 G
    'larmor_order'          : 1, 
    'T'                     : 2.5e-6,   #in s. See main paper Fig. 5  
}