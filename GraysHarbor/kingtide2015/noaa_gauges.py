
from pylab import *
from scipy.interpolate import interp1d
from clawpack.geoclaw import util
import datetime

def fetch_gauges(make_plots=False):
    station_westport = '9441102'
    MLW_westport = -1.068 # m relative to MTL
    MHW_westport = 1.068  # m relative to MTL
    
    station_aberdeen = '9441187'
    MLW_aberdeen = -1.21  # m relative to MTL
    MHW_aberdeen = 1.21  # m relative to MTL
    
    begin_date = datetime.datetime(2015,12,21)
    end_date = datetime.datetime(2015,12,26)
    time_zone = 'GMT'
    datum = 'MTL'
    units = 'metric'
    cache_dir = '.'
    
    time_westport, eta_westport, etap_westport = \
        util.fetch_noaa_tide_data(station_westport, begin_date, end_date,
                                  time_zone, datum, units, cache_dir)
    
    time_aberdeen, eta_aberdeen, etap_aberdeen = \
        util.fetch_noaa_tide_data(station_aberdeen, begin_date, end_date,
                                  time_zone, datum, units, cache_dir)

    # Save eta and eta prime for ocean forcing:
    # convert to elapsed seconds:
    dt_westport = time_westport - time_westport[0]
    t_westport = array([dt.item().total_seconds() for dt in dt_westport])

    # lag time so forcing in ocean reaches Westport at desired time:
    dt_ocean_to_westport = 2700.   # from tests
    t_ocean = t_westport - dt_ocean_to_westport
    t = t_ocean

    # decrease amplitude in ocean based on observed growth:
    etap_ocean = etap_westport / 1.1
    eta = etap_ocean

    dt_aberdeen = time_aberdeen - time_aberdeen[0]
    t_aberdeen = array([dt.item().total_seconds() for dt in dt_aberdeen])

    etap_ocean_fcn = interp1d(t_ocean,etap_ocean,kind='linear')
    print('Ocean level at t=0 is %.3f m relative to %s' \
            % (etap_ocean_fcn(0), datum))


    if make_plots:
        figure(250, figsize=(13,6))
        clf()
        plot(time_westport, etap_westport, 'b', label='Westport %s' % station_westport)
        plot(time_westport, MHW_westport*ones(time_westport.shape),'b--')
        plot(time_westport, MLW_westport*ones(time_westport.shape),'b--')

        plot(time_aberdeen, etap_aberdeen, 'r', label='Aberdeen %s' % station_aberdeen)
        plot(time_aberdeen, MHW_aberdeen*ones(time_aberdeen.shape),'r--')
        plot(time_aberdeen, MLW_aberdeen*ones(time_aberdeen.shape),'r--')

        grid(True)
        xticks(rotation=20)
        legend()
        title('NOAA predicted tides in Grays Harbor')
        ylabel('meters relative to %s' % datum)

        figure(251, figsize=(13,6))
        clf()
        plot(t_westport, etap_westport, 'b', label='Westport %s' % station_westport)
        plot(t_westport, MHW_westport*ones(t_westport.shape),'b--')
        plot(t_westport, MLW_westport*ones(t_westport.shape),'b--')

        plot(t_aberdeen, etap_aberdeen, 'r', label='Aberdeen %s' % station_aberdeen)
        plot(t_aberdeen, MHW_aberdeen*ones(t_aberdeen.shape),'r--')
        plot(t_aberdeen, MLW_aberdeen*ones(t_aberdeen.shape),'r--')

        plot(t_ocean, etap_ocean, 'k--', label='Ocean forcing')

        grid(True)
        xticks(rotation=20)
        legend()
        title('NOAA predicted tides in Grays Harbor')
        ylabel('meters relative to %s' % datum)


    return t_ocean,etap_ocean,t_westport,etap_westport,t_aberdeen,etap_aberdeen


def make_tidedata(t, eta, fname):
    # difference:
    etaprime = (eta[2:]-eta[:-2])/(t[2:]-t[:-2])  # central diff at interior t
    etaprime = hstack(((eta[1]-eta[0])/(t[1]-t[0]), etaprime, \
                       (eta[-1]-eta[-2])/(t[-1]-t[-2])))

    d = vstack((t, eta, etaprime)).T

    mt = d.shape[0]
    header = '%i  # mt' % mt
    savetxt(fname, d, header=header, comments='')
    print('Created ',fname)


if __name__=='__main__':

    t_ocean,etap_ocean,t_westport,etap_westport,t_aberdeen,etap_aberdeen = \
        fetch_gauges(make_plots=True)

    fname = 'tidedata_dec2015.txt'
    make_tidedata(t_ocean, etap_ocean, fname)
    


