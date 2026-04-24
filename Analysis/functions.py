import numpy as np

def find_lamb_wave(field, lat_val, lon_vals, extreme_type, tol):
    # Identify the point of the Lamb wave on a Hovmuller diagram!
    # This is for a certain lon value
    # Define a tolerance above which we assume we have 

    # Inputs:
    # field: field that is being investigated. Assume this is of the form (time, lat, lon)
    # lon: Indices of longitude value(s) to examine
    # extreme_type: max or min
    # tol: A positive value above/below which to stop the values looked at.

    # Output:
    # times: the times the Lamb wave reaches the lon values 

    if extreme_type == 'max':
        tol = -tol

    # Determine the range of values to look at
    time = field.time
    
    lamb_times = np.zeros(len(lon_vals))

    #print(len(lon_vals))

    
    for i in np.arange(len(lon_vals)):

        lon = lon_vals[i]
        #print('\n', lon)
        
        # Convert data to a 1D array.
        data = field[:, lat_val, lon]

        # Determine the time range to investigate
        idx = 0
        while (data[idx] - tol) < 0:
            idx += 1

        #print(idx)

        # Find the min/max over this time
        if extreme_type == 'min':
            arg_min = np.argmin(np.asarray(data[0:idx+1]))
            lamb_times[i] = time[arg_min]
            #lamb_times[i] = time[idx]
            #print(arg_min)
        elif extreme_type == 'max':
            lamb_times[i] = np.argmax(np.asarray(data[0:idx+1]))

    lamb_times_s = lamb_times*24*60*60
    dist_m = 6371220*(270-lon_vals)*np.pi/180.0
    slope, _ = np.polyfit(lamb_times_s, dist_m, 1)

    return lamb_times, slope
        

    