import numpy as np

from astropy.io import ascii
from astropy.coordinates import SkyCoord
from astropy.table import Table
import astropy.units as u

import sunpy.map
import sunpy.coordinates 
from sunpy.physics.differential_rotation import _calc_P_B0_SD

import matplotlib.pyplot as plt
plt.ioff()
import matplotlib.animation as an
import pandas as pd
import pickle
from glob import glob

from sunpy.database import Database
from sunpy.net import vso
from datetime import timedelta
tenmins = timedelta(minutes=10)

#set up the database of stereo files
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
stereo_boxes = []
eis_coords = []

nn = 0
# begin a loop over the EIS pointing data
for i, atime in enumerate(time):
    B0 = _calc_P_B0_SD(atime)['b0']


    # center of EIS field of view
    cen = SkyCoord(eis_pdata['XCEN'].quantity[i], eis_pdata['YCEN'].quantity[i],
               frame='helioprojective', dateobs=atime, B0=B0)
    eis_coords.append(cen)

    ## transforms cen into hgs
    cen_hgs = cen.transform_to('heliographic_stonyhurst')

    if any(np.isnan((cen_hgs.lon.value, cen_hgs.lat.value))):
        continue
    
    if cen_hgs.lon > 10*u.deg or cen_hgs.lon < -30*u.deg:
        pass
        #raise ValueError("Out of Chesse error, get a better satellite")
    else:
        res = db.query(vso.attrs.Time(atime, atime+tenmins))
        if not res:
            continue

        stereo_map = sunpy.map.Map(res[0])
        # get coords for the field of view box form EIS
        x_box = [eis_pdata['XCEN'].quantity[i] - eis_pdata['FOVX'].quantity[i]/2.0,
                 eis_pdata['XCEN'].quantity[i] + eis_pdata['FOVX'].quantity[i]/2.0]
        y_box = [eis_pdata['YCEN'].quantity[i] - eis_pdata['FOVY'].quantity[i]/2.0,
                 eis_pdata['YCEN'].quantity[i] + eis_pdata['FOVY'].quantity[i]/2.0]

        coords = zip(x_box, y_box)

        b_coord = SkyCoord(coords, frame='helioprojective', B0=B0, dateobs=atime)

        # EIS coords transform
        bhgs = b_coord.transform_to('heliographic_stonyhurst')
        stereo_box_hpc = bhgs.transform_to(stereo_map.coordinate_frame)

        check = [stereo_box_hpc.Tx.value, stereo_box_hpc.Ty.value]
        if np.any(np.isnan(check)):
            continue

        # stereo coords transform
        stereo_hpc = cen_hgs.transform_to(stereo_map.coordinate_frame)
        # horrific hack
        stereo_map.stereo_box = stereo_box_hpc
        stereo_map.eis_coords = cen
        stereo_map.units = stereo_map.spatial_units
        maps.append(stereo_map)
        stereo_coords.append(stereo_hpc)
        stereo_boxes.append(stereo_box_hpc)
        eis_times.append(atime)
        
        nn += 1
        if nn > 200:
            break

# create an ascii table
my_table = Table([stereo_coords, stereo_boxes, eis_times], 
                 names=['stereo_coords', 'stereo_boxes', 'eis_times'])

f = open("/storage2/EUVI/data_table.pik", 'wb')
pickle.dump(my_table,f)
f.close()


# begin the figure
fig = plt.figure()
ax = plt.subplot(projection=stereo_map)

im = maps[0].plot(axes=ax)
removes = []

# routine for defining the coords
def run(i):
    global removes
    while removes:
        removes.pop(0).remove()
    
    # update the data
    amap = maps[i]
    im.set_array(amap.data)
    
    # make the rectangle
    w = (stereo_boxes[i][1].Tx - stereo_boxes[i][0].Tx)
    h = (stereo_boxes[i][1].Ty - stereo_boxes[i][0].Ty)
    
    rect = amap.draw_rectangle(u.Quantity([stereo_boxes[i][0].Tx, stereo_boxes[i][0].Ty]),
                               w, h, transform=ax.get_transform('world'))
        
    ax.set_title('STEREO EUVI 30.4 nm {}'.format(eis_times[i]))
    text = ax.text(10, amap.data.shape[1] - 100,
            "EIS Pointing ({}, {})".format(eis_coords[i].Tx, eis_coords[i].Ty),
            color='white')
    
    removes += rect
    removes.append(text)
    

ani = an.FuncAnimation(fig, run, np.arange(len(eis_times)), interval=50)
ani.save('eis_trans.mp4')
plt.show()

fig = plt.figure()
def overplot(fig, ax, sunpy_map):
    box = sunpy_map.stereo_box
    # make the rectangle
    w = (box[1].Tx - box[0].Tx)
    h = (box[1].Ty - box[0].Ty)
    rect = sunpy_map.draw_rectangle(u.Quantity([box[0].Tx, box[0].Ty]),
                               w, h, transform=ax.get_transform('world'))
    cen = sunpy_map.eis_coords
    text = ax.text(10, sunpy_map.data.shape[1] - 100,
                   "EIS Pointing ({}, {})".format(cen.Tx, cen.Ty),
                   color='white')
    rect.append(text)
    return rect

mc = sunpy.map.Map(maps, cube=True)
mc.peek(plot_function=overplot)
plt.show()



# #plotting routines
# #plot the stereo images
# fig = plt.figure()
# ax = plt.subplot(projection=stereo_map.wcs):
# im = stereo_map.plot()
# plt.plot(stereo_hpc.Tx.to(u.deg), stereo_hpc.Ty.to(u.deg), 'o',
#         transform=ax.get_transform('world'))
# put the box on the images
# w = (EIS_hpc[1].Tx - EIS_hpc[0].Tx).to(u.deg).value
# h = (EIS_hpc[1].Ty - EIS_hpc[0].Ty).to(u.deg).value
# rect = plt.Rectangle((EIS_hpc[0].Tx.to(u.deg).value, EIS_hpc[0].Ty.to(u.deg).value), #
#        w, h, color='white', fill=False, transform=ax.get_transform('world'))
# ax.add_artist(rect)
# plt.show()






