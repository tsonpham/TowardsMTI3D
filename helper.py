from operator import itemgetter
from obspy.core.util import AttribDict
from obspy import read_events
from obspy.geodetics import gps2dist_azimuth
from mpl_toolkits.basemap import Basemap
from pyproj import Geod
import numpy as np

def make_stats(event=None, station=None):
    """
    Map event and station object to stats with attributes.

    :param event: ObsPy `~obspy.core.event.event.Event` object
    :param station: station object with attributes latitude, longitude and
        elevation
    :return: ``stats`` object with station and event attributes
    """
    ## support function 1
    def __get_event_origin_prop(h):
        def wrapper(event):
            try:
                r = (event.preferred_origin() or event.origins[0])[h]
            except IndexError:
                raise ValueError('No origin')
            if r is None:
                raise ValueError('No origin ' + h)
            if h == 'depth':
                r = r / 1000
            return r
        return wrapper
    ## support function 2
    def __get_event_magnitude(event):
        try:
            return (event.preferred_magnitude() or event.magnitudes[0])['mag']
        except IndexError:
            raise ValueError('No magnitude')
    ## support function 3
    def __get_event_id(event):
        evid = event.get('resource_id')
        if evid is not None:
            evid = str(evid)
        return evid
    ## cataloge 1
    _STATION_GETTER = (('station', itemgetter('name')),
                       ('station_latitude', itemgetter('latitude')),
                       ('station_longitude', itemgetter('longitude')))
    _SACHEADER_GETTER = (('station', itemgetter('kstnm')),
                         ('station_latitude', itemgetter('stla')),
                         ('station_longitude', itemgetter('stlo')))
    ## cataloge 2
    _EVENT_GETTER = (
        ('event_latitude', __get_event_origin_prop('latitude')),
        ('event_longitude', __get_event_origin_prop('longitude')),
        ('event_depth', __get_event_origin_prop('depth')),
        ('event_magnitude', __get_event_magnitude),
        ('event_time', __get_event_origin_prop('time')),
        ('event_id', __get_event_id))
    ## process
    stats = AttribDict({})
    if event is not None:
        for key, getter in _EVENT_GETTER:
            stats[key] = getter(event)
    if station is not None:
        try:
            for key, getter in _STATION_GETTER:
                stats[key] = getter(station)
        except KeyError:
            for key, getter in _SACHEADER_GETTER:
                stats[key] = getter(station)
    ## calculate 
    dist, baz, az = gps2dist_azimuth(stats.station_latitude,
                                    stats.station_longitude,
                                    stats.event_latitude,
                                    stats.event_longitude)
    stats.update({'distance': dist/1e3, 'back_azimuth': baz, 'azimuth': az})
    return stats

def prepare_stats(quakeml,stationtxt,t0,vred,window):
    ## event time
    event = read_events('ga2021sqogij.xml')[0]
    event_time = (event.preferred_origin() or event.origins[0])['time']

    ## time window of interest
    t0 = 20
    vred = 7.5
    window = 120

    ## form a list of STATS objects from event and input configuration file
    objstats = []
    with open('stations.txt') as fp:
        for line in fp:
            if line.strip().startswith('#'): continue
            ## create a station object having a name, latitude and longitude
            tokens = line.split()
            tmp = {'name': tokens[0], 'latitude': float(tokens[1]),
                   'longitude': float(tokens[2])}
            s = make_stats(event, tmp)
            s.update({'t0':t0,'vred':vred,'window':window})
            objstats.append(s)
    
    return objstats

def plot_lune_frame(ax, frame_color='k', grid_color='lightgray', 
                    clvd_left=True, clvd_right=True, lon_0=0):
    g = Geod(ellps='sphere')
    bm=Basemap(projection='hammer',lon_0=lon_0)
    ## Make sure that the axis has equal aspect ratio
    ax.set_aspect('equal')
    ## Plot meridian grid lines
    lats = np.arange(-90, 91)
    for lo in range(-30, 31, 10):
        lons = np.ones(len(lats)) * lo
        x, y = bm(lons, lats)
        ax.plot(x, y, lw=0.5, c=grid_color)
    x0, y = bm(-30*np.ones(len(lats)), lats)
    x1, y = bm(30*np.ones(len(lats)), lats)
    # ax.fill_betweenx(y, x0, x1, lw=0)
    ## Plot the left most meridian boundary
    lons = np.ones(len(lats)) * -30
    x, y = bm(lons, lats)
    ax.plot(x, y, lw=1, color=frame_color)
    ## Plot the right most meridian boundary
    lons = np.ones(len(lats)) * 30
    x, y = bm(lons, lats)
    ax.plot(x, y, lw=1, c=frame_color)
    ## Plot parallel grid lines
    lons = np.arange(-30, 31)
    for la in range(-90, 91, 10):
        lats = np.ones(len(lons)) * la
        x, y = bm(lons, lats)
        ax.plot(x, y, lw=0.5, c=grid_color)
    ## Put markers on special mechanism
    #-- isotropic points
    x, y = bm(0, 90)
    ax.plot(x, y, 'o', c=frame_color, ms=2)
    ax.annotate('+ISO', xy=(x, y*1.03), ha='center')
    ax.plot(x, 0, 'o', c=frame_color, ms=2)
    ax.annotate('-ISO', xy=(x, -y*.03), ha='center', va='top')
    #-- CLVD points
    x, y = bm(30, 0)
    ax.plot(x, y, 'o', c=frame_color, ms=2)
    if clvd_right: ax.annotate('-CLVD', xy=(1., 0.5), xycoords='axes fraction', \
                rotation='vertical', va='center')
    x, y = bm(-30, 0)
    ax.plot(x, y, 'o', c=frame_color, ms=2)
    if clvd_left: ax.annotate('+CLVD', xy=(0, 0.5), xycoords='axes fraction', \
                rotation='vertical', ha='right', va='center')
    # -- Double couple point
    x, y = bm(0, 0)
    ax.plot(x, y, 'o', c=frame_color, ms=2)
    ax.annotate('DC', xy=(x*.99, y*.99), ha='right', va='top')
    # -- LVD
    lvd_lon = 30
    lvd_lat = np.degrees(np.arcsin(1/np.sqrt(3)))
    x, y = bm(-lvd_lon, lvd_lat)
    ax.plot(x, y, 'o', c=frame_color, ms=2)
    x, y = bm(lvd_lon, 90-lvd_lat)
    ax.plot(x, y, 'o', c=frame_color, ms=2)
    arc = g.npts(-lvd_lon, lvd_lat, lvd_lon, 90-lvd_lat, 50)
    x, y = bm([p[0] for p in arc], [p[1] for p in arc])
    ax.plot(x, y, lw=1, c=frame_color)

    x, y = bm(-lvd_lon, lvd_lat-90)
    ax.plot(x, y, 'o', c=frame_color, ms=2)
    x, y = bm(lvd_lon, -lvd_lat)
    ax.plot(x, y, 'o', c=frame_color, ms=2)
    arc = g.npts(-lvd_lon, lvd_lat-90, lvd_lon, -lvd_lat, 50)
    x, y = bm([p[0] for p in arc], [p[1] for p in arc])
    ax.plot(x, y, lw=1, c=frame_color)
               
    return bm
    
def _mt2lune(mxx, myy, mzz, mxy, mxz, myz):
    m33 = np.array([[mxx, mxy, mxz], [mxy, myy, myz], [mxz, myz, mzz]])
    eivals = np.linalg.eigvals(m33)
    eivals.sort()
    
    ## lune longitude calculated from the eigen value triple
    nom = -eivals[0] - eivals[2] + 2 * eivals[1]
    den = np.sqrt(3) * (eivals[2]- eivals[0])
    gamma = np.arctan2(nom, den) / np.pi * 180

    ## lune latitude calculated from the eigen value triple
    nom = np.sum(eivals)
    den = np.sqrt(3) * np.sqrt(np.sum(eivals**2))
    beta = np.arccos(nom / den) / np.pi * 180

    ## orientation angles determined from the eigen vector triple
    return gamma, 90 - beta
mt2lune = np.vectorize(_mt2lune)