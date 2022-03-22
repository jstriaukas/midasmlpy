# 3-11 work content
from datetime import datetime
from unicodedata import numeric

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

def monthBegin(x):
    x = corr_datetype(x)
    df = datetime.strptime(x, '%Y-%m-%d %H:%M:%S')
    df = df.replace(day = 1, hour = 0, minute = 0, second = 0)
    return df

import calendar
def monthEnd(x):
    x = corr_datetype(x)
    df = datetime.strptime(x, '%Y-%m-%d %H:%M:%S')
    year = df.year
    month = df.month + 1
    last = calendar.monthrange(year, month)
    last_day = last[1]
    df = df.replace(day = last_day, hour = 0, minute = 0, second = 0)
    return df

def is_na(data):
    data_temp = []
    for i, y in enumerate(data):
        if y is not None:
            data_temp.append(y)
        else: continue
    return data_temp
    
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

def mixed_freq_data(data_y, data_ydate, data_x, data_xdate, x_lag, y_lag, horizon, est_start, est_end, disp_flag = True):
    data_y = is_na(data_y)
    data_ydate = is_na(data_ydate)
    data_x = is_na(data_x)
    data_xdate = is_na(data_xdate)
    data_ydate_vec = data_ydate
    data_ydate_vec = corr_datetype(data_ydate_vec)
    for i in range(len(data_ydate_vec)):
        data_ydate_vec[i] = datetime.strptime(data_ydate_vec[i], '%Y-%m-%d %H:%M:%S')
    data_xdate_vec = data_xdate
    data_xdate_vec = corr_datetype(data_xdate_vec)
    for i in range(len(data_xdate_vec)):
        data_xdate_vec[i] = datetime.strptime(data_xdate_vec[i], '%Y-%m-%d %H:%M:%S')
    est_start = corr_datetype(est_start)
    est_start = datetime.strptime(est_start[0], '%Y-%m-%d %H:%M:%S')
    est_end = corr_datetype(est_end)
    est_end = datetime.strptime(est_end[0], '%Y-%m-%d %H:%M:%S')
    data_xdate_vec = date_vec(data_xdate_vec)
    data_ydate_vec = date_vec(data_ydate_vec)
    data_ydate_num = data_ydate
    data_xdate_num = data_xdate
    date_format = ["year(s)", "month(s)", "day(s)", "hour(s)", "minute(s)", "second(s)"]
    period_y = data_freq(data_ydate_vec)["period"]
    unit_y = data_freq(data_ydate_vec)["unit"]
    period_x = data_freq(data_xdate_vec)["period"]
    unit_x = data_freq(data_xdate_vec)["unit"]
    y_lag = lag_num(y_lag, period_y, unit_y)
    x_lag = lag_num(x_lag, period_x, unit_x)
    horizon = lag_num(horizon, period_x, unit_x)
    if y_lag < 0:
        print("y_lag cann't be negative.")
        quit()
    if x_lag < 0:
        print("x_lag cann't be negative")
        quit()
    min_date_y = data_ydate_num[y_lag + 1]
    min_date_x = data_xdate_num[max(0, x_lag + horizon - 1)]
    if min_date_y > min_date_x:
        min_date = min_date_y
    else:
        min_date = min_date_x
    max_date_y = data_ydate_num[len(data_ydate_num) - 1]
    max_date_x = data_xdate_num[len(data_xdate_num) - 1]
    if horizon < 0:
        max_date_x = data_xdate_vec[len(data_xdate_vec)[0]][:]
        max_date_x[unit_x] = max_date_x[unit_x] + period_x * horizon
        max_date_x = str(max_date_x[0]) + "-" +str(max_date_x[1]) + "-" +str(max_date_x[2]) + " " +str(max_date_x[3]) + ":" +str(max_date_x[4]) + ":" +str(max_date_x[5])
        max_date_x = datetime.strptime(max_date_x, '%Y-%m-%d %H:%M:%S')
    if max_date_y > max_date_x:
        max_date = max_date_x
    else:
        max_date = max_date_x
    if est_start is None:
        est_start = min_date
    else:
        if est_start < min_date:
            print("Start date cannot be earlier than possible due to lagged regressors. Reset start date to most recent possible.")
            est_start = min_date
    if est_end is None:
        est_end = max_date
    else:
        if est_end == max_date:
            print("Terminal date cannot be later than largest date account for lags. Reset to largest date.")
            est_end = max_date
    tol = 1e-10
    for i in range(len(data_ydate_num)):
        if data_ydate_num[i] >= est_start:
            loc_start = i
            break
    for i in range(len(data_ydate_num)):
        if data_ydate_num[i] >= est_end:
            loc_end = i
            break
    est_y = data_y[loc_start : loc_end]
    est_ydate = data_ydate_num[loc_start : loc_end]
    for i in range(len(data_ydate_num)):
        if data_ydate_num[i] >= max_date:
            loc_forecast_end = i
            break
    if loc_end + 1 <= loc_forecast_end:
        out_y = data_y[loc_end + 1: loc_forecast_end + 1]
        out_ydate = data_ydate_num[loc_end + 1: loc_forecast_end + 1]
        n_forecast = len(out_y)
    else:
        out_y, out_ydate = None
        n_forecast = 0
    nobs = loc_end - loc_start
    est_lag_y = [[None for x in range(y_lag)] for y in range(nobs)]
    est_lag_ydate = [[None for x in range(y_lag)] for y in range(nobs)]
    for m in range(y_lag):
        est_lag_y[:][m] = data_y[loc_start - m - 1:loc_end - m]
        est_lag_ydate[:][m] = data_ydate_num[loc_start - m - 1:loc_end - m]
    if loc_end + 1 <= loc_forecast_end:
        out_lag_y = [[None for x in range(y_lag)] for y in range(n_forecast)]
        out_lag_ydate = [[None for x in range(y_lag)] for y in range(n_forecast)]
        for m in range(y_lag):
            out_lag_y[:][m] = data_y[loc_end - m:loc_forecast_end - m - 1]
            out_lag_ydate[:][m] = data_ydate_num[loc_end - m:loc_forecast_end - m - 1]
    else:
        out_lag_y, out_lag_ydate = None
    est_x = [[None for x in range(y_lag)] for y in range(nobs)]
    est_xdate = [[None for x in range(y_lag)] for y in range(nobs)]
    for t in range(nobs):
        for i in range(len(data_xdate_num)):
            if data_xdate_num[i] >= est_ydate[t]:
                loc = i
                break
        if loc == None:
            loc = len(data_xdate_num)
        if loc - horizon > len(data_x):
            nobs = t - 1
            est_y = est_y[:nobs]
            est_ydate = est_ydate[:nobs]
            est_lag_y = est_lag_y[:nobs]
            est_lag_ydate = est_lag_ydate[:nobs]
            est_x = est_x[:nobs]
            est_xdate = est_xdate[:nobs]
            max_date = est_ydate[len(est_ydate)]
            print("Horizon is a large negative number. Observations are further truncated to max date possible")
            break
        else:
            est_x[t][:] = data_x[loc - horizon - 1 : loc -horizon - x_lag : -1]
            est_xdate[t][:] = data_xdate_num[loc - horizon - 1 : loc -horizon - x_lag : -1]
    if loc_end + 1 <= loc_forecast_end:
        out_x = [[None for x in range(x_lag)] for y in range(n_forecast)]
        out_xdate = [[None for x in range(x_lag)] for y in range(n_forecast)]
        for t in range(n_forecast):
            for i in range(len(data_xdate_num)):
                if data_xdate_num[i] >= out_ydate[t]:
                    loc = i
                    break
            if loc is None:
                loc = len(data_xdate_num)
            if loc - horizon > len(data_x):
                n_forecast = t- 1
                out_y = out_y[:n_forecast]
                out_ydate = out_ydate[:n_forecast]
                out_lag_y = out_lag_y[:n_forecast]
                out_lag_ydate = out_lag_ydate[:n_forecast]
                out_x = out_x[:n_forecast]
                out_xdate = out_xdate[:n_forecast]
                break
            else:
                out_x[t][:] = data_x[loc - horizon - 1: loc - horizon - x_lag: -1]
                out_xdate[t][:] = data_xdate_num[loc - horizon - 1: loc - horizon - x_lag: -1]
    else:
        out_x, out_xdate = None
    if disp_flag:
        print("Frequency of Data Y:", period_y, date_format[unit_y])
        print("Frequency of Data X:", period_x, date_format[unit_x])
        print("Start Date: ", est_start)
        print("Start Date: ", est_end)
        print("Mixed frequency regression time frame:")
        for m in [0, 1, nobs - 1]:
            print("Reg Y(", est_ydate[m], ")`s on: ")
            if y_lag == 1:
                print("Y(", est_lag_ydate[m][0], ")`s")
            if y_lag == 2:
                print("Y(", est_lag_ydate[m][0], ")`s Y(", est_lag_ydate[m][len(est_lag_ydate)[1]], ")`s")
            if y_lag == 3:
                print("Y(", est_lag_ydate[m][0],")`s Y(", est_lag_ydate[m][1],")`s ... Y(", est_lag_ydate[m][len(est_lag_ydate)[1]], ")`s")
            if x_lag == 1:
                print("X(", est_xdate[m][:], ")`s")
            if x_lag == 2:
                print("X(", est_xdate[m][0], ")`s X(", est_xdate[m][len(est_xdate)[1]], ")`s")
            if x_lag == 3:
                print("X(", est_xdate[m][0], ")`s X(", est_xdate[m][1], ")`s ... X(", est_xdate[m][len(est_xdate)[1]], ")`s")
    return {
        "est_y": est_y,
        "est_ydate": est_ydate,
        "est_x": est_x,
        "est_xdate": est_xdate,
        "est_lag_y": est_lag_y,
        "est_lag_ydate": est_lag_ydate,
        "out_y": out_y,
        "out_ydate": out_ydate,
        "out_x": out_x,
        "out_xdate": out_xdate,
        "out_lag_y": out_lag_y,
        "out_lag_ydate": out_lag_ydate,
        "x_lag": x_lag,
        "y_lag": y_lag,
        "min_date": min_date,
        "max_date": max_date
    }
    
def mixed_freq_data_single(data_refdate, data_x, data_xdate, x_lag, horizon, est_start, est_end, disp_flag = True):
    data_refdate = is_na(data_refdate)
    data_x = is_na(data_x)
    data_xdate = is_na(data_xdate)
    data_refdate = corr_datetype(data_refdate)
    data_xdate = corr_datetype(data_xdate)
    est_start = corr_datetype(est_start)
    est_end = corr_datetype(est_end)
    for i in range(len(data_refdate)):
        data_refdate[i] = datetime.strptime(data_refdate[i], '%Y-%m-%d %H:%M:%S')
    for i in range(len(data_xdate)):
        data_xdate[i] = datetime.strptime(data_xdate[i], '%Y-%m-%d %H:%M:%S')
    data_refdate_vec = data_refdate
    data_xdate_vec = data_xdate
    est_start = datetime.strptime(est_start[0], '%Y-%m-%d %H:%M:%S')
    est_end = datetime.strptime(est_end[0], '%Y-%m-%d %H:%M:%S')
    data_refdate_vec = date_vec(data_refdate_vec)
    data_xdate_vec = date_vec(data_xdate_vec)
    data_refdate_num = data_refdate
    data_xdate_num = data_xdate
    
    date_format = ['year(s)','month(s)','day(s)','hour(s)','minute(s)','second(s)']
    period_ref = data_freq(data_refdate_vec)["period"]
    unit_ref = data_freq(data_refdate_vec)["unit"]
    period_x = data_freq(data_xdate_vec)["period"]
    unit_x = data_freq(data_xdate_vec)["unit"]
    
    ref_lag = lag_num(1, period_ref, unit_ref)
    x_lag = lag_num(x_lag, period_x, unit_x)
    horizon = lag_num(horizon, period_x, unit_x)
    if ref_lag < 0:
        print('ref.lag cannot be negative.')
        quit()
    if x_lag < 0:
        print('x.lag cannot be negative')
        quit()
    # Minimum and maximum dates that data support
    min_date_ref = data_refdate_num[ref_lag]
    min_date_x = data_xdate_num[max(0, x_lag + horizon - 1)]
    if min_date_ref > min_date_x:
        min_date = min_date_ref
    else:
        min_date = min_date_x
    max_date_ref = data_refdate_num[len(data_refdate_num) - 1]
    max_date_x = data_xdate_num[len(data_xdate_num) - 1]
    if horizon < 0:
        max_date_x = data_xdate_vec[len(data_xdate_vec)[0] - 1][:]
        max_date_x[unit_x] = max_date_x[unit_x] + period_x * horizon
        max_date_x = str(max_date_x[0]) + "-" +str(max_date_x[1]) + "-" +str(max_date_x[2]) + " " +str(max_date_x[3]) + ":" +str(max_date_x[4]) + ":" +str(max_date_x[5])
        max_date_x = datetime.strptime(max_date_x, '%Y-%m-%d %H:%M:%S')
    if max_date_ref > max_date_x:
        max_date = max_date_x
    else:
        max_date = max_date_ref
    # Check and set default sample period
    if est_start is None:
        est_start = min_date
    else:
        if est_start < min_date:
            est_start = min_date
    if est_end is None:
        est_end = max_date
    else:
        if est_end > max_date:
            print('Terminal date cannot be later than largest date accounting for lags. Reset to largest date possible: ', max_date)
            est_end = max_date
    # Construct reference date data
    tol = 1e-10
    loc_start = 0
    for i in range(len(data_refdate_num)):
        
        if data_refdate_num[i] >= est_start:
            loc_start = i
            break
    for i in range(len(data_refdate_num)):
        if data_refdate_num[i] >= est_end:
            loc_end = i
            break
    est_refdate = data_refdate_num[loc_start:loc_end]
    
    for i in range(len(data_refdate_num)):
        if data_refdate_num[i] >= max_date:
            loc_forecast_end = i
            break
    if loc_end + 1 <= loc_forecast_end:
        out_refdate = data_refdate_num[loc_end: loc_forecast_end]
        n_forecast = len(out_refdate)
    else:
        out_refdate = None
        n_forecast = len(out_refdate)
    nobs = loc_end - loc_start
    est_x = [[None for x in range(x_lag)] for y in range(nobs)]
    est_xdate = [[None for x in range(x_lag)] for y in range(nobs)]
    for t in range(nobs):
        for i in range(len(data_xdate_num)):
            x = data_xdate_num[i]
            y = est_refdate[t]
            if x >= y:
                loc = i
                break
        if loc is None:
            loc = len(data_xdate_num)
        if loc - horizon > len(data_x):
            nobs = t - 1
            est_refdate = est_refdate[:nobs]
            est_x = est_x[:nobs]
            est_xdate = est_xdate[:nobs]
            max_date = est_refdate[len(est_refdate)]
            print("Horizon is a large negative number. Observations are further truncated to max date possible", max_date)
            break
        else:
            est_x[t][:] = data_x[loc - horizon - 1 : loc -horizon - x_lag : -1]
            est_xdate[t][:] = data_xdate_num[loc - horizon - 1 : loc -horizon - x_lag : -1]
    if loc_end + 1 <= loc_forecast_end:
        out_x = [[None for x in range(x_lag)] for y in range(n_forecast)]
        out_xdate = [[None for x in range(x_lag)] for y in range(n_forecast)]
        for t in range(n_forecast):            
            for i in range(len(data_xdate_num)):
                if data_xdate_num[i] >= out_refdate[t]:
                    loc = i
                    break
            if loc is None:
                loc = len(data_xdate_num)
            if loc - horizon > len(data_x):
                n_forecast = t - 1
                out_refdate = out_refdate[:n_forecast]
                out_x = out_x[:n_forecast]
                out_xdate = out_xdate[n_forecast]
                break
            else:
                out_x[t][:] = data_x[loc - horizon - 1 : loc - horizon - x_lag : -1]
                out_xdate[t][:] = data_xdate_num[loc - horizon - 1 : loc - horizon - x_lag : -1]
    else:
        out_x, out_xdate = None
    if disp_flag:
        # Display mixed frequency data
        print('Frequency of Reference Date:', period_ref, date_format[unit_ref])
        print('Frequency of Data X:',period_x,date_format[unit_x]) 
        print('Start Date: ', est_start) 
        print('Terminal Date: ', est_end) 
        
        # Display timeframe of mixed frequency regression
        for m in [0, 1, nobs - 1]:
            print("Reg date(", est_refdate[m], ")`s on: ")
            if x_lag == 1:
                print("X(", est_xdate[m][0], ")`s")
            if x_lag == 2:
                print("X(", est_xdate[m][0], ")`s X(", est_xdate[m][10], ")`s ...")
            if x_lag == 3:
                print("X(", est_xdate[m][0],")`s Y(", est_xdate[m][1],")`s X(", est_xdate[m][10], ")`s ...")
            if x_lag > 3:
                print("X(", est_xdate[m][0],")`s Y(", est_xdate[m][1],")`s ... X(", est_xdate[m][10], ")`s ...")
           
    return {
        "est_refdate": est_refdate,
        "est_x": est_x,
        "est_xdate": est_xdate,
        "out_refdate": out_refdate,
        "out_x": out_x,
        "out_xdate": out_xdate,
        "x_lag": x_lag,
        "min_date": min_date,
        "max_date": max_date,
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