import numpy as np

from astropy.io import ascii
from astropy.coordinates import SkyCoord
from astropy.table import Table
import astropy.units as u

import sunpy.map
import sunpy.coordinates 
from sunpy.physics.differential_rotation import _calc_P_B0_SD

import matplotlib.pyplot as plt
import matplotlib.animation as an
import pandas as pd
import pickle
from glob import glob

from sunpy.database import Database
from sunpy.net import vso
from datetime import timedelta
tenmins = timedelta(minutes=10)

#set up the database
db = Database('sqlite:///sunpydb.sqlite', default_waveunit=u.AA)

# read in the stereo files
files = glob('/storage2/EUVI/preped/*.fts')

# read in the eis pointing data
eis_pdata = ascii.read("eis_pointing_table.txt")
eis_pdata['XCEN'].unit = u.arcsec
eis_pdata['YCEN'].unit = u.arcsec
eis_pdata['FOVX'].unit = u.arcsec
eis_pdata['FOVY'].unit = u.arcsec
eis_pdata.sort('DATE_OBS')

time = [sunpy.time.parse_time(tt.strip()) for tt in eis_pdata['DATE_OBS']]

#things for plotting
maps = []
eis_times = []
stereo_coords = []
eis_boxes = []

nn = 0
# begin a loop over the EIS pointing data
for i, atime in enumerate(time):
    B0 = _calc_P_B0_SD(atime)['b0']

    
    # center of EIS field of view
    cen = SkyCoord(eis_pdata['XCEN'].quantity[i], eis_pdata['YCEN'].quantity[i],
               frame='helioprojective', dateobs=atime, B0=B0)
               
    ## transforms cen into hgs
    hgs = cen.transform_to('heliographicstonyhurst')

    print "Longitude is %r" % hgs.lon.value
    if hgs.lon > 10*u.deg or  hgs.lon < -30*u.deg:
        print "not in field of view"
        #'raise ValueError("Out of Chesse error, get a better satellite")

    else:
        res = db.query(vso.attrs.Time(atime, atime+tenmins))
        if not res:
            continue
        
        stereo_map = sunpy.map.Map(res[0])
        # get coords for the field of view box form EIS
        x_box = [eis_pdata['XCEN'].quantity[0] - eis_pdata['FOVX'].quantity[0]/2.0,
                 eis_pdata['XCEN'].quantity[0] + eis_pdata['FOVX'].quantity[0]/2.0]
        y_box = [eis_pdata['YCEN'].quantity[0] - eis_pdata['FOVY'].quantity[0]/2.0,
                 eis_pdata['YCEN'].quantity[0] + eis_pdata['FOVY'].quantity[0]/2.0]
    
        coords = zip(x_box, y_box)
    
        b_coord = SkyCoord(coords, frame = 'helioprojective', B0=B0, dateobs = atime)
    
        # EIS coords transform
        bhgs = b_coord.transform_to('heliographicstonyhurst')
        bhgs.B0 = stereo_map.heliographic_latitude
        bhgs.L0 = stereo_map.heliographic_longitude
        bhgs.D0 = stereo_map.dsun
        EIS_hpc = bhgs.transform_to('helioprojective')


        # stereo coords transform
        print nn
        nn += 1
        print atime
        hgs.B0 = stereo_map.heliographic_latitude
        hgs.L0 = stereo_map.heliographic_longitude
        hgs.D0 = stereo_map.dsun
        stereo_hpc = hgs.transform_to('helioprojective')
        maps.append(stereo_map)
        stereo_coords.append(stereo_hpc)
        eis_boxes.append(EIS_hpc)
        eis_times.append(atime)

        if nn > 200:
            break

# create an ascii table
my_table = Table([stereo_coords, eis_boxes, eis_times], 
                 names=['stereo_coords', 'eis_boxes', 'eis_times'])

f = open("/storage2/EUVI/data_table.pik", 'wb')
pickle.dump(my_table,f)
f.close()


# begin the figure
fig = plt.figure()
ax = plt.subplot(projection=stereo_map.wcs)

im = maps[0].plot()
removes = []

## routine for defining the coords
#def run(i):
#    global removes
#    while removes:
#        removes.pop(0).remove()
#    # update the data
#    amap = maps[i]
#    stereo_coords[i]
#    eis_boxes[i]
#    eis_times[i]
#    im.set_array(amap.data)
#    # make the rectangle
#    w = (eis_boxes[i][1].Tx - eis_boxes[i][0].Tx).to(u.deg).value
#    h = (eis_boxes[i][1].Ty - eis_boxes[i][0].Ty).to(u.deg).value
#    rect = plt.Rectangle((eis_boxes[i][0].Tx.to(u.deg).value, eis_boxes[i][0].Ty.to(u.deg).value), 
#                         w, h, color='white', fill=False, transform=ax.get_transform('world'))
#    ax.add_artist(rect)
#    ax.set_title('Stereo EUVI 30.4 nm %s' % eis_times[i])
#    removes.append(rect)
#
##len(eis_boxes)
#ani = an.FuncAnimation(fig, run, np.arange(len(eis_times)), interval = 100)
#plt.show()



## plotting routines
## plot the stereo images
#fig = plt.figure()
#ax = plt.subplot(projection=stereo_map.wcs)
#im = stereo_map.plot()
#plt.plot(stereo_hpc.Tx.to(u.deg), stereo_hpc.Ty.to(u.deg), 'o',
#         transform=ax.get_transform('world'))
## put the box on the images
#w = (EIS_hpc[1].Tx - EIS_hpc[0].Tx).to(u.deg).value
#h = (EIS_hpc[1].Ty - EIS_hpc[0].Ty).to(u.deg).value
#rect = plt.Rectangle((EIS_hpc[0].Tx.to(u.deg).value, EIS_hpc[0].Ty.to(u.deg).value), #
#        w, h, color='white', fill=False, transform=ax.get_transform('world'))
#ax.add_artist(rect)
#plt.show()






