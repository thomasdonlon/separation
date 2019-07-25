#=======================================
#SEPARATION APPLICATION PLOTTING PROGRAM
#=======================================

#plot out a result similar to the last figure of Weiss et al. 2018
#functionality within python3 (and hopefully less obfuscated by omitting large sections of NewbyTools)

#wrapra turns the program from 360->0 ra plotting into 180->-180 ra plotting, for easier visualization of data that crosses the ra=0 axis

import numpy as np
import matplotlib.pyplot as plt
from astropy.coordinates import SkyCoord
from astropy import units as u

#=======================================
#CONSTANT DEFINITIONS
#=======================================

#enter the stripe numbers here in the same order as they appear in star_files & stream_files
stripes = [80, 81, 82, 83, 84, 85, 86]

figwidth = 20
figheight = 10

#=======================================
#HELPER FUNCTION DEFINITIONS (ADAPTED FROM NEWBYTOOLS)
#=======================================
deg = 180.0 / np.pi
rad = np.pi / 180.0
surveyCenterRa = 185.0
surveyCenterDec = 32.5
lCP = 122.932 * rad
raGP = 192.8594813 * rad
decGP = 27.1282511 * rad

def angle_bounds(ra, dec, wrapra=False):
    #inputs should be in degrees: easily fixable, but not necessary in this case
    to_array = False
    if type(ra) != type(np.array([])): #WARNING: assumes identical length for ra & dec, i.e. both are arrays or neither is
        #I know it's lazy coding but whatever, if it breaks contact me and yell at me
        ra = np.array([ra])
        dec = np.array([dec])
        to_array = True

    for i in range(len(ra)):
        if not(wrapra):
            if ra[i] < 0:
                ra[i] += 360
            elif ra[i] > 360:
                ra[i] -= 360
        else:
            if ra[i] > 180:
                ra[i] -= 360

    for i in range(len(dec)):
        if dec[i] < -90:
            dec[i] += 180
        elif dec[i] > 90:
            dec[i] -= 180

    if to_array:
        ra = ra[0]
        dec = dec[0]

    return ra, dec

def get_eta (wedge):
    """ Get the eta value that corresponds to the given stripe value """
    ss = 2.5
    if wedge <= 46:  eta = wedge * ss - 57.5
    else:  eta = wedge * ss - 57.5 - 180.0
    return eta

def GCToEq (mu_deg, nu_deg, wedge, wrapra=False):
    """ Converts Stripe mu, nu into equatorial ra, dec.  Called 'atGCToEq' in at SurveyGeometry.c"""

    node = (surveyCenterRa - 90.0)*rad
    eta = get_eta(wedge)
    inc = (surveyCenterDec + eta)*rad
    mu, nu = (mu_deg*rad), (nu_deg*rad)
    # Rotation
    x2 = np.cos(mu - node)*np.cos(nu)
    y2 = np.sin(mu - node)*np.cos(nu)
    z2 = np.sin(nu)
    x1 = x2
    y1 = y2*np.cos(inc) - z2*np.sin(inc)
    z1 = y2*np.sin(inc) + z2*np.cos(inc)
    ra = np.arctan2(y1,x1) + node
    dec = np.arcsin(z1)
    ra, dec = angle_bounds(ra*deg, dec*deg, wrapra=wrapra)
    return ra, dec

def EqToGC (ra_deg, dec_deg, wedge, wrapra=False):  #produces lists...  anglebounds2!!!  ROTATE
    """ Converts equatorial ra,dec into Great Circle mu, nu; 'atSurveyGeometry.c' in
    m31.phys.rpi.edu:/p/prd/astrotools/v5_18/Linux-2-4-2-3-2/src"""
    node = (surveyCenterRa - 90.0)*rad
    eta = get_eta(wedge)
    inc = (surveyCenterDec + eta)*rad
    ra, dec = (ra_deg*rad), (dec_deg*rad)
    # Rotation
    x1 = np.cos(ra-node)*np.cos(dec)
    y1 = np.sin(ra-node)*np.cos(dec)
    z1 = np.sin(dec)
    x2 = x1
    y2 = y1*np.cos(inc) + z1*np.sin(inc)
    z2 = -y1*np.sin(inc) + z1*np.cos(inc)
    mu = np.arctan2(y2,x2) + node
    nu = np.arcsin(z2)
    nu, mu = angle_bounds(nu*deg, mu*deg, wrapra=wrapra)
    return mu,nu

def numtocolor(c, stripe):
    out = []

    if stripe == 80:
        for n in c:
            if n == 0:
                out.append('k')
            elif n == 1:
                out.append('magenta')
            elif n == 2:
                out.append('g')
            elif n == 3:
                out.append('r')
            else:
                out.append('black')
    elif stripe == 81:
        for n in c:
            if n == 0:
                out.append('k')
            elif n == 1:
                out.append('b')
            elif n == 2:
                out.append('orange')
            elif n == 3:
                out.append('g')
            else:
                out.append('r')
    elif stripe == 82:
        for n in c:
            if n == 0:
                out.append('k')
            elif n == 1:
                out.append('orange')
            elif n == 2:
                out.append('r')
            elif n == 3:
                out.append('b')
            else:
                out.append('g')
    elif stripe == 83:
        for n in c:
            if n == 0:
                out.append('k')
            elif n == 1:
                out.append('g')
            elif n == 2:
                out.append('black')
            elif n == 3:
                out.append('teal')
            else:
                out.append('b')
    elif stripe == 84:
        for n in c:
            if n == 0:
                out.append('k')
            elif n == 1:
                out.append('black')
            elif n == 2:
                out.append('g')
            elif n == 3:
                out.append('magenta')
            else:
                out.append('b')
    elif stripe == 85:
        for n in c:
            if n == 0:
                out.append('k')
            elif n == 1:
                out.append('r')
            elif n == 2:
                out.append('gray')
            elif n == 3:
                out.append('g')
            else:
                out.append('b')
    else:
        for n in c:
            if n == 0:
                out.append('k')
            elif n == 1:
                out.append('b')
            elif n == 2:
                out.append('orange')
            elif n == 3:
                out.append('black')
            else:
                out.append('purple')

    return out

#=======================================
#FUNCTION DEFINITIONS
#=======================================

def readStarFile(f):
    print('reading in star file ',f)
    data = []

    f = open(f)
    f.readline() #first line is number of stars in that stripe

    #turn data in file into separated floats
    for line in f:
        line = line.split(' ')

        for i in range(3):
            line[i] = float(line[i].strip())

        data.append(line)

    return data

def readColorList(f, stripe):
    print('reading in output file ',f)
    data = []

    f = open(f)

    #turn data in file into separated floats
    for line in f:
        line = line.split(' ')

        data.append(int(line[0]))

    data = numtocolor(data, stripe)
    data = np.array(data)

    return data

def plotStarData(star_data, c_list, mode="radec", wrapra=False):

    ra = np.array([])
    dec = np.array([])
    r = np.array([])
    c_list_stripe = []
    n = 0

    for i in range(len(stripes)):
        wedge = stripes[i]

        stripe = star_data[i]
        stripe_T = np.transpose(stripe)

        l = stripe_T[0]
        b = stripe_T[1]
        r_tmp = stripe_T[2]

        #convert l, b to ra, dec
        co = SkyCoord(l*u.degree, b*u.degree, frame='galactic')
        co = co.transform_to('icrs')
        ra_tmp = co.ra.value
        dec_tmp = co.dec.value

        #remove all ra/dec from array if they overlap the next stripe
        if i < len(stripes) - 1: #do not perform on last stripe, there is nothing below it
            mu, nu = EqToGC(ra_tmp, dec_tmp, wedge+1, wrapra=wrapra)

            for j in range(len(mu)):
                if nu[j] < 1.25: #point is overlapping the below stripe
                    ra_tmp[j] = -9999
                    dec_tmp[j] = -9999
                    r_tmp[j] = -9999
                else:
                    c_list_stripe.append(c_list[n])
                n+=1
        else:
            for j in range(len(ra_tmp)):
                c_list_stripe.append(c_list[n])
                n+=1

        index = np.argwhere(ra_tmp==-9999)
        ra_tmp = np.delete(ra_tmp, index)
        dec_tmp = np.delete(dec_tmp, index)
        r_tmp = np.delete(r_tmp, index)

        ra = np.append(ra, ra_tmp)
        dec = np.append(dec, dec_tmp)
        r = np.append(r, r_tmp)
        #-------------------------------------------------------

    ra, dec = angle_bounds(ra, dec, wrapra=wrapra)

    c_list_stripe = np.array(c_list_stripe)

    index = np.argwhere(c_list_stripe != 'r')
    #index = np.append(index, np.argwhere(c_list != 'b'))
    ra = np.delete(ra, index)
    dec = np.delete(dec, index)
    r = np.delete(r, index)
    c_list_stripe = np.delete(c_list_stripe, index)

    if mode == "radec":
        #plt.scatter(ra, dec, c=c_list_stripe, marker='.', s=1)
        plt.hist2d(ra, dec, bins=[200, 100], cmap='binary')
    elif mode == "radist":
        #plt.scatter(ra, r, c=c_list_stripe, marker='.', s=1)
        plt.hist2d(ra, r, bins=[200, 100], cmap='binary')
    else:
        print("error <plotStarData()>: mode must be 'radec' or 'radist'")
#=======================================
#RUNTIME
#=======================================

#read in the star data
star_data = [] #array of the data for every stripe
c_list = []
"""
[ <overall list>
    [ <one stripe>
        [ <one star>
            mu, nu, r (at least I think it's mu, nu, r, might be wrong but I'll check that later)
        ] ...
    ] ...
]
"""
for stripe in stripes:
    star_data.append(readStarFile('/home/donlot/Desktop/separation/separation_data/StarsSouth/stars-'+str(stripe)+'_2.txt'))
    c_list = np.append(c_list, readColorList('/home/donlot/mwah_nbody/build_sep/bin/sep_output/separation_output_'+str(stripe)+'.out', stripe))

star_data = np.array(star_data)

fig = plt.figure(figsize=(figwidth, figheight))

plotStarData(star_data, c_list, mode='radec', wrapra=True)
plt.title('Separation Application: South Redo 3')
plt.xlabel('RA')
plt.ylabel('Dec')
plt.xlim((45, -15))
plt.ylim((-12.5, 7.5))
plt.savefig("/home/donlot/Desktop/separation/figures/South_Redo_3_separation_ra_dec_bif_density.png")

fig = plt.figure(figsize=(figwidth, figheight))

plotStarData(star_data, c_list, mode='radist', wrapra=True)
plt.title('Separation Application: South Redo 3')
plt.xlabel('RA')
plt.ylabel('Dist')
plt.xlim((45, -15))
plt.ylim((0,45))
plt.savefig("/home/donlot/Desktop/separation/figures/South_Redo_3_separation_ra_dist_bif_density.png")
