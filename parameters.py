params = {}

params['carbon_params'] = {
    'coupling'  : 80e3, # in Hz
    'T2'        : 60e-3, # in seconds
}

params['nv_params'] = {
    'average_repump_time'   : 220e-9, # in seconds 
    'repump_time_jitter'    : 0.0,    # in seconds
    'repump_time_offset'    : 0.0,    # in seconds 
    'pflip'                 : 0.5,    #
    'mw_infidelity'         : 0.008,
    'init_infidelity'       : 0.001,
}

params['simulation_params'] = {
    'entangling_attempts'   : 2000,
    'repetitions'           : 1000,
    'larmor_period'         : 2.256e-6, # in s; @ magnetic Field B = 414 G
    'larmor_order'          : 1, 
    'T'                     : 2.5e-6,   #in s. See main paper Fig. 5
    'do_carbon_pi'          : False,    
}