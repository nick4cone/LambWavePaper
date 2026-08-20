import numpy as np
import scipy
from scipy import optimize

def lamb_speed_by_wavefront(ps, time, lon, lon_idxs):
    # Script tuned for test case 1, with the pressure perturbation.
    # Compute at the equator, so lat_ind should reflect this
    # for the given resolution
    lamb_times = np.zeros(len(lon_idxs))

    for i in np.arange(len(lon_idxs)):

        lon_idx = lon_idxs[i]

        # Time series of ps at this point
        ps_t = ps[:, lon_idx]

        # Find the global maximum here.
        max_idx = np.argmax(ps_t)

        # Find the first negative beyond the maximum
        neg_idx = max_idx
        while ps_t[neg_idx] > 0:
            neg_idx += 1

        # Find the zero crossing with linear regression
        slope = (ps_t[neg_idx]-ps_t[neg_idx-1])/(time[neg_idx]-time[neg_idx-1])
        lamb_times[i] = time[neg_idx-1] - ps_t[neg_idx-1]/slope

    # Now, fit to all the lamb times
    lon_vals = lon[lon_idxs]

    lamb_times_s = lamb_times*24*60*60
    dist_m = 6371220*(270-lon_vals)*np.pi/180.0
    slope, _ = np.polyfit(lamb_times_s, dist_m, 1)

    return lamb_times, slope

def lamb_wave_by_zero_crossing_ps(ps, time, lon, lon_idxs):
    # Script tuned for test case 1, with the pressure perturbation.
    # Compute at the equator, so lat_ind should reflect this
    # for the given resolution
    lamb_times = np.zeros(len(lon_idxs))

    for i in np.arange(len(lon_idxs)):

        lon_idx = lon_idxs[i]

        # Time series of ps at this point
        ps_t = ps[:, lon_idx]

        # Find the global maximum here.
        max_idx = np.argmax(ps_t)

        # Find the first negative beyond the maximum
        neg_idx = max_idx
        while ps_t[neg_idx] > 0:
            neg_idx += 1

        # Find the zero crossing with linear regression
        slope = (ps_t[neg_idx]-ps_t[neg_idx-1])/(time[neg_idx]-time[neg_idx-1])
        lamb_times[i] = time[neg_idx-1] - ps_t[neg_idx-1]/slope

    # Now, fit to all the lamb times
    lon_vals = lon[lon_idxs]

    lamb_times_s = lamb_times*24*60*60
    dist_m = 6371220*(270-lon_vals)*np.pi/180.0
    slope, _ = np.polyfit(lamb_times_s, dist_m, 1)

    return lamb_times, slope

        

def find_lamb_wave(field, time, lon, lon_idxs):
    # Identify the slope the Lamb wave on a Hovmuller diagram!
    # This determines the speed (given some factors)
    # This is for a certain lon value
    # Define a tolerance above which we assume we have 

    # Inputs:
    # field: field that is being investigated. Assume this is of the form (time, lon)
    # so, we have already selected the equator (lat=0) to examine
    # lon_vals: Indices of longitude value(s) to examine
    # extreme_type: max or min
    # tol: A positive value above/below which to stop the values looked at.

    # Output:
    # times: the times the Lamb wave reaches the lon values 

    #if extreme_type == 'max':
    #    tol = -tol

    # Determine the range of values to look at
    #time = field.time
    
    lamb_times = np.zeros(len(lon_idxs))

    for i in np.arange(len(lon_idxs)):

        lon_idx = lon_idxs[i]
        
        # Convert data to a 1D array.
        data = field[:, lon_idx]

        # Determine the global max at this point.
        max_idx = np.argmax(data)
        lamb_times[i] = time[max_idx]

    lon_vals = lon[lon_idxs]

    lamb_times_s = lamb_times*24*60*60
    dist_m = 6371220*(270-lon_vals)*np.pi/180.0
    slope, _ = np.polyfit(lamb_times_s, dist_m, 1)

    return lamb_times, slope

# Compute speeds for models with a pressure vertical coordinate
def p_coord_speeds(zTs, H, g, kappa):
    # zTs: Range of model top heights
    # H: Isothermal scale height

    roots = np.zeros_like(zTs)
    heights = np.zeros_like(zTs)
    speeds = np.zeros_like(zTs)

    # Positive value as a guess for root solving
    x0 = 0.5

    # Tolearance on interval to avoid singularities!
    eps = 1e-6
    
    for i in np.arange(len(zTs)):
        zT = zTs[i]

        # Define the residual
        if zT < (8*H/3):
            def R_func(x):
                return np.tan(x) - 4*kappa*H*x/(zT*(1 - 2*kappa + (2*H*x/zT)**2))

            #x_root = scipy.optimize.fsolve(R_func, x0)
            
            x_root = scipy.optimize.root_scalar(R_func, bracket = [eps,np.pi/2-eps], method='brentq')
            x_root = x_root.root
            
            #x_root = scipy.optimize.brentq(R_func, 0.01, np.pi/2)
            
            # Back out the height
            h_star = 4*kappa*H/((2*H*x_root/zT)**2 + 1)
            
        else:
            def R_func(x):
                return np.tanh(x) - 4*kappa*H*x/(zT*(1 - 2*kappa - (2*H*x/zT)**2))
        
            # Compute the singularity point 
            x_up = zT * np.sqrt(1-2*kappa)/(2*H)
                
            #x_root = scipy.optimize.fsolve(R_func, x0)
            x_root = scipy.optimize.root_scalar(R_func, bracket = [eps, x_up-eps], method='brentq')
            x_root = x_root.root
            
            #x_root = scipy.optimize.brentq(R_func, 0.01, np.pi/2)
    
            # Back out the height
            h_star = 4*kappa*H/(1 - (2*H*x_root/zT)**2)

    

        # Comptue the speed
        roots[i] = x_root
        heights[i] = h_star
        speeds[i] = np.sqrt(g*h_star)


    return roots, heights, speeds

    