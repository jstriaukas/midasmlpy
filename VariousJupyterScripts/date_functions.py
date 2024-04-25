# 3-11 work content
from datetime import datetime
from unicodedata import numeric
from sklearn.preprocessing import StandardScaler
import numpy as np

# Finds time difference between two given times
#time input style: [year,month,day,hour,min,sec]
def diff_time_mf(time1, time2, origin, units = ["auto", "secs", "mins", "hours", "days", "weeks"]):
    ini_time1 = time1
    ini_time2 = time2
    if(len(ini_time1) != len(ini_time2)):
        print("time list length doesn't equal!")
        quit()
    result = []
    for i in range(len(ini_time1)):
        time1 = ini_time1[i]
        time2 = ini_time2[i]
        try:
            time1 = datetime.strptime(time1, '%m-%d-%y')
            time2 = datetime.strptime(time2, '%m-%d-%y')
        except:
            time1 = datetime.strptime(time1, '%m/%d/%y')
            time2 = datetime.strptime(time2, '%m/%d/%y')            
        delta_time = time1 - time2
        delta_days = delta_time.days
        delta_sec = delta_time.seconds
        if units == "secs":
            result.append(86400 * delta_days + delta_sec)
        elif units == "mins":
            result.append(1440 * delta_days + delta_sec / 60)
        elif units == "hours":
            result.append(24 * delta_days + delta_sec / 3600)
        elif units == "days":
            result.append(delta_days + delta_sec / 86400)
        elif units == "weeks":
            result.append(delta_days / 7 + delta_sec / 604800)
        else:
            result.append(86400 * delta_days + delta_sec)
    return result

# Lag a number by a given period and unit
def lag_num(x_lag, period, unit):
    if type(x_lag) == int or float:
        return x_lag
    try:
        multiplier = float(x_lag[:-1])
    except:
        print('The description of lags cannot be recognized. The format should be 3m, 1q, etc')
        quit()
    #Convert multiplier to daily frequency (business days)
    ndaysPerYear = 264
    ndaysPerQuarter = 66
    ndaysPerMonth = 22
    nhoursPerDay = 8
    if x_lag[-1] == 'y':
        multiplier = multiplier * ndaysPerYear
    if x_lag[-1] == 'q':
        multiplier = multiplier * ndaysPerQuarter
    if x_lag[-1] == 'm':
        multiplier = multiplier * 1
    if x_lag[-1] == 'h':
        multiplier = multiplier / nhoursPerDay
    if x_lag[-1] == 's':
        multiplier = multiplier /  (nhoursPerDay*60*60)
    if unit == 1:
        x_lag = round(multiplier / (ndaysPerYear * period))
    if unit == 2:
        x_lag = round(multiplier / (ndaysPerMonth * period))
    if unit == 3:
        x_lag = round(multiplier / period)
    if unit == 4:
        x_lag = round(multiplier / (period / nhoursPerDay))
    if unit == 5:
        x_lag = round(multiplier / (period / nhoursPerDay / 60))
    if unit == 6:
        x_lag = round(multiplier / (period / nhoursPerDay / 60 / 60))
    return x_lag

# Calculate the mode of a given dataset
def mode_midasml(data):
    nobs = len(data)
    data = data.tolist()
    data.sort()
    count = 1
    countMax = 1
    dataMode = data[0]
    for t in range(1, nobs):
        if data[t] == data[t - 1]:
            count += 1
        else:
            if count > countMax:
                countMax = count
                dataMode = data[t - 1]
            count = 1
    if count > countMax:
        countMax = count
        dataMode = data[nobs - 1]
    return {
        "dataMode": dataMode,
        "countMax": countMax
    }
    
import statistics
import numpy as np

# Calculate the frequency of data points in a given date vector
def data_freq(DateVec):
    ini_DV = DateVec
    DateVec = np.array(DateVec)
    shape_vec = DateVec.shape
    for i in range(len(DateVec) - 1):
        DateVec[i]=DateVec[i + 1] - DateVec[i]
        total_days = DateVec[i][0] * 360 + DateVec[i][1] * 30 + DateVec[i][2]
        if total_days >= 0:
            if DateVec[i][2] < 0:
                DateVec[i][2] += 30
                DateVec[i][1] -= 1
                if DateVec[i][1] < 0:
                    DateVec[i][1] += 12
                    DateVec[i][0] -=1
            elif DateVec[i][1] < 0:
                DateVec[i][1] += 12
                DateVec[i][0] -=1
            
    DateVec = np.array(DateVec[:-1])
    y = DateVec.shape
    x= DateVec[:,0]
    # Check annual or lower frequency
    modeUse = mode_midasml(DateVec[:,0])
    modeUse = modeUse["dataMode"]
    if modeUse >= 1:
        period = modeUse
        unit = 1
        return {
            "period": period,
            "unit": unit
        }
        
    # Check monthly frequency, quarter = 3 months, semiannual = 6 months
    modeUse = mode_midasml(DateVec[:,1])
    modeUse = modeUse["dataMode"]
    if modeUse < 0:
        modeUse = modeUse + 12
    if modeUse >= 1:
        period = modeUse
        unit =2
        return {
            "period": period,
            "unit": unit
        }
        
    # Check daily frequency, week = 7 days, biweekly = 14 days
    modeUse = mode_midasml(DateVec[:,2])
    modeUse = modeUse["dataMode"]
    if modeUse < 0:
        modeUse += 30
    if modeUse >= 1:
        period = modeUse
        unit = 3
        return {
            "period": period,
            "unit": unit
        }
    
    # Check hourly frequency
    modeUse = mode_midasml(DateVec[:,3])
    modeUse = modeUse["dataMode"]
    if modeUse < 0:
        modeUse += 24
    if modeUse >= 1:
        period = modeUse
        unit = 4
        return {
            "period": period,
            "unit": unit
        }
        
  # Check minutely frequency
    modeUse = mode_midasml(DateVec[:,4])
    modeUse = modeUse["dataMode"]
    if modeUse < 0:
        modeUse += 60
    if modeUse >= 1:
        period = modeUse
        unit = 5
        return {
            "period": period,
            "unit": unit
        }
  
  # Check secondly frequency  
    elapse = diff_time_mf(ini_DV[:,5], ini_DV[:,5], origin = "1970-01-01", units = "secs")
    period = statistics.mean(elapse)
    unit = 6
    return {
        "period": period,
        "unit": unit
    }

# Check if two dates match   
def dateMatch(x, y):
    n = len(x)
    x_out = n * [0]
    for i in range(n):
        i_x = x[i]
        if i_x in y:
            x_out[i] = i_x
        else:
            i_match = False
            i_x_k = i_x
            k = 1
            while i_match == False and k <= 5:
                i_x_k = datetime.strptime(i_x_k, '%m-%d-%y')
                i_x_k.day -= 1
                if i_x_k in y:
                    i_match = True
                k = k + 1
            if k < 5:
                x_out[i] = i_x_k
            else:
                x_out[i] = i_x
    return x_out

from datetime import datetime

# Check if a given date is in the correct format
# ISO Date Format: 2021-07-07 06:36:55
# Formatted Time: 23/06/21 10:33:02
def corr_datetype(date):
    for i in range(len(date)):
        x = date[i].split(" ")
        x_len = len(x)
        if x_len == 1:
            date[i] += " 00:00:00"
        x = date[i].split(" ")
        if "/" in x[0]:
            x_conv = x[0].split("/")
            year = x_conv[2]
            month = x_conv[0]
            day = x_conv[1]
            date[i] = year + "-" + month + "-" + day + " " + x[1]
    return date

# Calculate the beginning of a month for a given date
def monthBegin(x):
    x = corr_datetype(x)
    df = datetime.strptime(x, '%Y-%m-%d %H:%M:%S')
    df = df.replace(day = 1, hour = 0, minute = 0, second = 0)
    return df

import calendar
from datetime import datetime
from unicodedata import numeric
import statistics
import numpy as np
from datetime import datetime
import calendar

# Calculate the end of a month for a given date
def monthEnd(x):
    x = corr_datetype(x)
    df = datetime.strptime(x, '%Y-%m-%d %H:%M:%S')
    year = df.year
    month = df.month + 1
    last = calendar.monthrange(year, month)
    last_day = last[1]
    df = df.replace(day = last_day, hour = 0, minute = 0, second = 0)
    return df

# Check if a given dataset contains any missing values
def is_na(data):
    data_temp = []
    for i, y in enumerate(data):
        if y is not None:
            data_temp.append(y)
        else: continue
    return data_temp

# Convert a string into a date vector    
def date_vec(s):
    w = 6
    h = len(s)
    mat = [[0 for x in range(w)] for y in range(h)]
    for i in range(h): 
        mat[i][0] = s[i].year
        mat[i][1] = s[i].month
        mat[i][2] = s[i].day
        mat[i][3] = s[i].hour
        mat[i][4] = s[i].minute
        mat[i][5] = s[i].second
    return mat

# Handle mixed frequency data

def mixed_freq_data(data_y, data_ydate, data_x, data_xdate, x_lag, y_lag, horizon):
    """
    Function to process mixed frequency data.

    Args:
        data_y (list): List of values for the dependent variable.
        data_ydate (list): List of dates for the dependent variable.
        data_x (list): List of values for the independent variable.
        data_xdate (list): List of dates for the independent variable.
        x_lag (int): Number of lags for the independent variable.
        y_lag (int): Number of lags for the dependent variable.
        horizon (int): Forecast horizon.

    Returns:
        dict: Dictionary containing the processed data and other information.
    """

    # Preprocess data_y and data_ydate
    data_y = is_na(data_y)
    data_ydate = is_na(data_ydate)
    data_ydate_vec = corr_datetype(data_ydate)
    for i in range(len(data_ydate_vec)):
        data_ydate_vec[i] = datetime.strptime(data_ydate_vec[i], '%Y-%m-%d %H:%M:%S')
    data_ydate_vec = date_vec(data_ydate_vec)

    # Preprocess data_x and data_xdate
    data_x = is_na(data_x)
    data_xdate = is_na(data_xdate)
    data_xdate_vec = corr_datetype(data_xdate)
    for i in range(len(data_xdate_vec)):
        data_xdate_vec[i] = datetime.strptime(data_xdate_vec[i], '%Y-%m-%d %H:%M:%S')
    data_xdate_vec = date_vec(data_xdate_vec)

    # Get the numerical representation of data_ydate and data_xdate
    data_ydate_num = data_ydate
    data_xdate_num = data_xdate

    # Define date format units
    date_format = ["year(s)", "month(s)", "day(s)", "hour(s)", "minute(s)", "second(s)"]

    # Get the period and unit of data_ydate_vec and data_xdate_vec
    period_y = data_freq(data_ydate_vec)["period"]
    unit_y = data_freq(data_ydate_vec)["unit"]
    period_x = data_freq(data_xdate_vec)["period"]
    unit_x = data_freq(data_xdate_vec)["unit"]

    # Convert lag values to appropriate units
    y_lag = lag_num(y_lag, period_y, unit_y)
    x_lag = lag_num(x_lag, period_x, unit_x)
    horizon = lag_num(horizon, period_x, unit_x)

    if y_lag < 0:
        print("y_lag cann't be negative.")
        quit()
    if x_lag < 0:
        print("x_lag cann't be negative")
        quit()

    # Get the number of observations
    nobs = len(data_ydate_num)

    # Initialize the values for data_y and data_ydate
    est_y = data_y
    est_ydate = data_ydate_num

    # Initialize the lagged values for data_y and data_ydate
    est_lag_y = [[None for x in range(y_lag)] for y in range(nobs)]
    est_lag_ydate = [[None for x in range(y_lag)] for y in range(nobs)]

    # Initialize the values for data_x and data_xdate
    est_x = [[None for x in range(x_lag)] for y in range(nobs)]
    est_xdate = [[None for x in range(x_lag)] for y in range(nobs)]
    for t in range(nobs):
        # Create lagged values of y
        est_lag_y[t][:] = data_y[t - 1 : t - y_lag -1: -1]
        est_lag_ydate[t][:] = data_ydate_num[t - 1  : t - y_lag -1: -1]

        # Create lagged values of x
        for i in range(len(data_xdate_num)):
            if data_xdate_num[i] >= data_ydate_num[t]:
                loc = i
                break
        if loc == None:
            loc = len(data_xdate_num)
        
        else:
            est_x[t][:] = data_x[loc - horizon  : loc -horizon - x_lag : -1]
            est_xdate[t][:] = data_xdate_num[loc - horizon  : loc -horizon - x_lag : -1]

    # All values of y dont always have values for x, so the values returning an empty list of x are removed in the lists
    indices_to_keep = [i for i, row in enumerate(est_x) if row]

    # Filter all lists based on these indices
    est_x = [est_x[i] for i in indices_to_keep]
    est_xdate = [est_xdate[i] for i in indices_to_keep]
    est_y = [est_y[i] for i in indices_to_keep]
    est_ydate = [est_ydate[i] for i in indices_to_keep]
    est_lag_y = [est_lag_y[i] for i in indices_to_keep]
    est_lag_ydate = [est_lag_ydate[i] for i in indices_to_keep]

    # Return the processed data and other information as a dictionary
    return {
        "est_y": est_y,
        "est_ydate": est_ydate,
        "est_x": est_x,
        "est_xdate": est_xdate,
        "est_lag_y": est_lag_y,
        "est_lag_ydate": est_lag_ydate
    }


def legendre_polynomial(x, d):
    """
    This function calculates the Legendre polynomial of order n for the input x.
    For degrees 2 and above, the function uses the recursive relation to calculate the polynomial.

    Args:
        x (float): The input value.
        d (int): The degree of the Legendre polynomial.
    
    Returns:
        float: The value of the Legendre polynomial at x for degree d.
    """
    if d == 0:
        return 1
    elif d == 1:
        return x
    else:
        return ((2 * d - 1) * x * legendre_polynomial(x, d - 1) - (d - 1) * legendre_polynomial(x, d - 2)) / d

# Create matrix of legendre polynomials for lags and degrees
def legendre_matrix_create(x_lags, legendre_degree = 3):
    """
    This function creates a matrix of Legendre polynomials for the input x_lags and legendre_degree.

    Args:
        x_lags (int): The number of lags.
        legendre_degree (int): The degree of the Legendre polynomial. Default is 3.
    
    Returns:
        numpy.ndarray: A matrix of Legendre polynomials for the given lags and degrees.
    """
    # Create list of equally spaced values between -1 and 1 for use in legendre polynomials
    x_values = np.linspace(-1, 1, num=x_lags)
    # Loop through the lags and degrees to create a lag x degree matrix of legendre polynomials
    legendre_matrix = [[legendre_polynomial(x, n) for n in range(legendre_degree)] for x in x_values]
    return np.array(legendre_matrix)

def data_transform(Y, Y_date, X, X_date, x_lags, y_lags, horizon, legendre_degree=3, standardize = True):
    """
    This function takes the Y and X data and creates the X_tilde matrix by multiplying the X matrix with the legendre matrix.
    
    Parameters:
    - Y (ndarray): The target variable
    - Y_date (ndarray): The dates of the target variable
    - X (ndarray): The predictors
    - X_date (ndarray): The dates of the predictors
    - x_lags (int): The number of lags to be created for the predictors
    - y_lags (int): The number of lags to be created for the target variable
    - horizon (int): The number of periods to forecast ahead
    - legendre (int): The degree of the legendre polynomial
    - standardize (bool): Whether to standardize the predictors or not

    Returns:
    - dict: A dictionary containing the Y, X_tilde, and Y_lagged matrices

    """
    # Ensure data structure
    Y_date = Y_date.astype(str)
    Y = Y
    X_date = X_date.astype(str)
    X = X.astype(np.float32)
    # initialize empty X_tilde
    X_tilde = []
    legendre = legendre_matrix_create(x_lags, legendre_degree)
    # Standardize inputs if needed
    if standardize:
        X = (X - np.mean(X, axis=0)) / np.std(X, axis=0)
    for i in range(X.shape[1]):
        # First get the correct data structure
        result = mixed_freq_data(Y, Y_date, X[:,i], X_date, x_lags, y_lags, horizon)
        X_matrix = np.array(result['est_x'], dtype=float)
        # Add the dot-product of X_matrix and the legendre matrix to the X_tilde matrix
        X_tilde.append(np.dot(X_matrix, legendre))
    # Create arrays
    X_tilde = np.array(X_tilde)
    Y = np.array(result['est_y'])
    Y_lagged = np.array(result['est_lag_y'])
    # convert 3d X_tilde array into 2d with shape (predictors*legendre) x observations 
    X_tilde = X_tilde.reshape(-1, X_tilde.shape[1]).T
    return {
        "Y": Y,
        "X_tilde": X_tilde,
        "Y_lagged": Y_lagged
    }

#if __name__ == "__main__":
#    with open("data_files/input_data_refdate.txt") as file_data_refdate:
#        data_refdate = [line.rstrip() for line in file_data_refdate]
#    with open("data_files/input_data_x.txt") as file_data_x:
#        data_x = [line.rstrip() for line in file_data_x]
#    with open("data_files/input_data_xdate.txt") as file_data_xdate:
#        data_xdate = [line.rstrip() for line in file_data_xdate]
#    x_lag = 12
#    horizon = 1
#    est_start = ["1990-01-01"]
#    est_end = ["2002-03-01"]
#    result = mixed_freq_data_single(data_refdate, data_x, data_xdate, x_lag, horizon, est_start, #est_end, disp_flag = True)
#    print(result)

#if __name__ == "__main__":
#    with open("input_data_x.txt") as file_data_x:
#        data_x = [line.rstrip() for line in file_data_x]
#    with open("input_data_xdate.txt") as file_data_xdate:
#        data_xdate = [line.rstrip() for line in file_data_xdate]
#    with open("input_data_y.txt") as file_data_x:
#        data_y = [line.rstrip() for line in file_data_x]
#    with open("input_data_ydate.txt") as file_data_xdate:
#        data_ydate = [line.rstrip() for line in file_data_xdate]
#    x_lag = 12
#    y_lag = 1
#    horizon = 1
#    est_start = ["1990-01-01"]
#    est_end = ["2002-03-01"]
#    result = mixed_freq_data(data_y, data_ydate, data_x, data_xdate, x_lag, y_lag, horizon, #est_start, est_end, disp_flag = True)
    
 #   print(result)

