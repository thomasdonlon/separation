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

#place star files here
star_files = [ \
            "/home/donlot/Desktop/separation/separation_data/StarsSouth/stars-80_2.txt", \
            "/home/donlot/Desktop/separation/separation_data/StarsSouth/stars-81_2.txt", \
            "/home/donlot/Desktop/separation/separation_data/StarsSouth/stars-82_2.txt", \
            "/home/donlot/Desktop/separation/separation_data/StarsSouth/stars-83_2.txt", \
            "/home/donlot/Desktop/separation/separation_data/StarsSouth/stars-84_2.txt", \
            "/home/donlot/Desktop/separation/separation_data/StarsSouth/stars-85_2.txt", \
            "/home/donlot/Desktop/separation/separation_data/StarsSouth/stars-86_2.txt" \
            ]

#place stream files here
#enter the stripe numbers here in the same order as they appear in star_files
stream_files = [ \
            "/home/donlot/Desktop/separation/separation_data/South_Redo_3/stream_opt_80_bundle_4_south4s_1.txt", \
            "/home/donlot/Desktop/separation/separation_data/South_Redo_3/stream_opt_81_bundle_4_south4s_1.txt", \
            "/home/donlot/Desktop/separation/separation_data/South_Redo_3/stream_opt_82_bundle_4_south4s_1.txt", \
            "/home/donlot/Desktop/separation/separation_data/South_Redo_3/stream_opt_83_bundle_4_south4s_1.txt", \
            "/home/donlot/Desktop/separation/separation_data/South_Redo_3/stream_opt_84_bundle_4_south4s_1.txt", \
            "/home/donlot/Desktop/separation/separation_data/South_Redo_3/stream_opt_85_bundle_4_south4s_1.txt", \
            "/home/donlot/Desktop/separation/separation_data/South_Redo_3/stream_opt_86_bundle_4_south4s_1.txt" \
                ]

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

def lbToEq (l_deg, b_deg, wrapra=False):
    """ Converts galactic l,b in to Equatorial ra, dec; from Binney and Merrifield, p. 31;
    l, b must be arrays of same shape"""
    l, b = (l_deg*rad), (b_deg*rad)
    # Conversion Code
    t = lCP - l
    dec = np.arcsin(np.sin(decGP)*np.sin(b) + np.cos(decGP)*np.cos(b)*np.cos(t) )
    ra = np.arctan2( (np.cos(b)*np.sin(t)),
                    ( (np.cos(decGP)*np.sin(b)) - (np.sin(decGP)*np.cos(b)*np.cos(t)))  )
    ra = ra + raGP
    ra, dec = angle_bounds(ra*deg, dec*deg, wrapra=wrapra)
    return ra, dec

def EqTolb (ra_deg, dec_deg):
    """ Converts equatorial ra, dec, into galactic l,b;  from Binney and Merrifield, p. 31
    following the method of http://star-www.st-and.ac.uk/~spd3/Teaching/AS3013/programs/radec2lb.f90
    NOT QUITE - I use arctan2 method instead"""
    ra, dec = (ra_deg*rad), (dec_deg*rad)
    ra, dec = angle_bounds(ra*deg, dec*deg) #removes wrapping
    # Conversion Code
    r = (ra - raGP)
    b = np.arcsin( np.sin(decGP)*np.sin(dec) + np.cos(decGP)*np.cos(dec)*np.cos(r) )
    t = np.arctan2((np.cos(dec)*np.sin(r)),
                   (np.cos(decGP)*np.sin(dec) - np.sin(decGP)*np.cos(dec)*np.cos(r)) )
    l = (lCP - t)
    b, l = angle_bounds((b*deg), (l*deg))
    return l, b

def thetaphi2draddec (theta, phi, ra, dec, r):
    #converts a theta, phi orientation and ra, dec, r location into a change in ra, dec that can be used for plotting directional arrows
    #theta, phi should be in rad, everything else should be in degrees
    #get xyz
    c = SkyCoord(ra*u.degree, dec*u.degree, frame='icrs')
    c = c.transform_to('galactic')
    l = c.l.value
    b = c.b.value
    x = r * np.cos(l*np.pi/180) * np.cos(b*np.pi/180) - 8
    y = r * np.sin(l*np.pi/180) * np.cos(b*np.pi/180)
    z = r * np.sin(b*np.pi/180)

    #move along xyz a short distance in the theta, phi direction
    dx = 0.1 #kpc, this can be adjusted to make the arrows more accurate if necessary
    x += dx * np.cos(phi) * np.sin(theta)
    y += dx * np.sin(phi) * np.sin(theta)
    z += dx * np.cos(theta)

    new_r = ((x+8)**2 + y**2 + z**2)**0.5

    #then translate that shifted point back into ra, dec
    l = np.arctan2(y, x + 8) * 180/np.pi
    b = np.arcsin(z/new_r) * 180/np.pi

    c = SkyCoord(l*u.degree, b*u.degree, frame='galactic')
    c = c.transform_to('icrs')
    new_ra = c.ra.value
    new_dec = c.dec.value
    dra = new_ra - ra
    ddec = new_dec - dec
    dr = new_r - r

    #make all the arrows the same length on the sky
    dlength = 1 #degrees
    norm = (dra**2 + ddec**2)**0.5
    dra *= (dlength/norm)
    ddec *= (dlength/norm)
    dr *= dlength/(dra**2 + (0.5*dr)**2)**0.5

    return dra, ddec, dr

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
            line[i] = float(line[i])

        data.append(line)

    return data

def readStreamFile(f):
    print('reading in stream file ',f)
    data = []

    f = open(f)
    #turn data in file into separated floats
    line = f.readline()
    line = line.split(',')
    line = line[4:len(line)-1] #omit first two elements and the last element of line

    n = 0
    stream = []
    for i in range(len(line)):
        line[i] = line[i].strip(' ').strip('[').strip(']') #clean up format
        line[i] = float(line[i]) #convert to float
        stream.append(line[i])

        if n == 5: #6 parameters to a stream
            stream = np.array(stream)
            data.append(stream) #does not alias

            stream = [] #reset counter
            n = -1

        n += 1

    return data

def plotStarData(star_data, mode="radec", wrapra=False):

    ra = np.array([])
    dec = np.array([])
    r = np.array([])

    for i in range(len(stripes)):
        wedge = stripes[i]

        stripe = star_data[i]
        stripe_T = np.transpose(stripe)

        l = stripe_T[0]
        b = stripe_T[1]
        r_tmp = stripe_T[2]

        #convert l, b to ra, dec
        ra_tmp, dec_tmp = lbToEq(l, b, wrapra=wrapra)

        #remove all ra/dec from array if they overlap the next stripe
        if i < len(stripes) - 1: #do not perform on last stripe, there is nothing below it
            mu, nu = EqToGC(ra_tmp, dec_tmp, wedge+1, wrapra=wrapra)

            for i in range(len(mu)):
                if nu[i] < 1.25: #point is overlapping the below stripe
                    ra_tmp[i] = -9999
                    dec_tmp[i] = -9999
                    r_tmp[i] = -9999

        index = np.argwhere(ra_tmp==-9999)
        ra_tmp = np.delete(ra_tmp, index)
        dec_tmp = np.delete(dec_tmp, index)
        r_tmp = np.delete(r_tmp, index)

        ra = np.append(ra, ra_tmp)
        dec = np.append(dec, dec_tmp)
        r = np.append(r, r_tmp)

    if mode == "radec":
        #plt.scatter(ra, dec, marker='.', c='k', s=1)
        plt.hist2d(ra, dec, bins=[200, 100], cmap='binary')
    elif mode == "radist":
        #plt.scatter(ra, r, marker='.', c='k', s=1)
        plt.hist2d(ra, r, bins=[200, 100], cmap='binary')
    else:
        print("error <plotStarData()>: mode must be 'radec' or 'radist'")

def plotStreamData(stream_data, mode="radec", wrapra=False, outline_streams=False):

    for i in range(len(stripes)):
        wedge = stripes[i]
        stripe = stream_data[i]

        for j in range(len(stripe)):
            epsilon = stripe[j][0]
            mu = stripe[j][1]
            r = stripe[j][2]
            theta = stripe[j][3]
            phi = stripe[j][4]
            sigma = stripe[j][5]

            ra, dec = GCToEq(mu, 0, wedge, wrapra=wrapra)
            dra, ddec, dr = thetaphi2draddec(theta, phi, ra, dec, r)

            if mode == "radec":
                if outline_streams and j == len(stripe)-1:
                    '''
                    #100 kpc limit
                    plt.plot([18.67, 16.47, 17.12, 16.00], [2.43, 0.0, -2.44, -4.91], c='r')
                    plt.plot([25.39, 24.15, 27.51, 23.14, 15.79, 15.87], [4.69, 0.0, -2.31, -4.75, -7.37, -9.82], c='b')
                    plt.plot([34.87, 33.98, 39.10], [2.17, 0.0, -2.07], c='g')
                    plt.plot([35.09, 37.04, 32.70, 29.31, 26.87, 21.90, 22.01], [4.33, 2.12, 0.0, -2.28, -4.64, -7.18, -9.57], c='magenta')
                    '''
                    #60 kpc limit
                    plt.plot([25.87, 17.54, 9.76, 354.26-360], [4.67, 2.44, 0.0, -7.37], c='r')
                    #plt.plot([34.58, 26.90, 5.80], [-4.35, -2.31, 4.99], c='teal')
                    plt.plot([37.43, 30.59, 27.48, 21.28, 14.58, 8.60], [2.11, 0.0, -2.27, -4.80, -7.39, -9.98], c='b')
                    plt.plot([36.89, 33.17, 29.34, 26.06, 22.8, 14.66], [2.12, 0.0, -2.27, -4.66, -7.14, -9.86], c='magenta')
                    plt.plot([18.93, 16.44, 14.23], [2.42, 0.0, -7.40], c='orange')
                plt.scatter(ra, dec, marker='o', c='navy', s=20)
                plt.arrow(ra, dec, dra, ddec, color='navy', head_width=0.2)
                #plt.annotate(str(j+(wedge-80)*4), (ra, dec), color='k', fontsize=25)
                plt.annotate(str(wedge) + '.' + str(j+1), (ra, dec), color='k', fontsize=10)
            elif mode == "radist":
                if outline_streams and j == len(stripe)-1:
                    '''
                    #100 kpc limit
                    plt.plot([18.67, 16.47, 17.12, 16.00], [25.28, 23.56, 30.84, 20.75], c='r')
                    plt.plot([25.39, 24.15, 27.51, 23.14, 15.79, 15.87], [23.87, 18.26, 20.43, 20.63, 20.66, 20.66], c='b')
                    plt.plot([34.87, 33.98, 39.10], [42.09, 46.32, 49.05], c='g')
                    plt.plot([35.09, 37.04, 32.70, 29.31, 26.87, 21.90, 22.01], [26.47, 27.53, 27.11, 26.11, 27.34, 28.04, 28.04], c='magenta')
                    '''
                    #60 kpc limit
                    plt.plot([25.87, 17.54, 9.76, 354.26-360], [23.67, 31.07, 40.16, 43.15], c='r')
                    #plt.plot([34.58, 26.90, 5.80], [59.19, 47.73, 45.65], c='teal')
                    plt.plot([37.43, 30.59, 27.48, 21.28, 14.58, 8.60], [25.80, 22.58, 21.09, 19.69, 20.32, 20.02], c='b')
                    plt.plot([36.89, 33.17, 29.34, 26.06, 22.8, 14.66], [27.23, 27.46, 25.48, 27.14, 26.97, 26.13], c='magenta')
                    plt.plot([18.93, 16.44, 14.23], [25.98, 26.84, 24.33], c='orange')
                plt.scatter(ra, r, marker='o', c='navy', s=20)
                plt.arrow(ra, r, dra, dr, color='navy', head_width=0.2)
                #plt.annotate(str(j+(wedge-80)*4), (ra, r), color='k', fontsize=25)
                plt.annotate(str(wedge) + '.' + str(j+1), (ra, r), color='k', fontsize=10)
            elif mode == "rasigma":
                if outline_streams and j == len(stripe)-1:
                    '''
                    #100 kpc limit
                    plt.plot([18.67, 16.47, 17.12, 16.00], [1.49, 1.50, 1.41, 0.19], c='r')
                    plt.plot([25.39, 24.15, 27.51, 23.14, 15.79, 15.87], [3.77, 5.75, 4.99, 5.40, 3.60, 3.60], c='b')
                    plt.plot([34.87, 33.98, 39.10], [24.97, 10.37, 24.77], c='g')
                    plt.plot([35.09, 37.04, 32.70, 29.31, 26.87, 21.90, 22.01], [1.44, 1.18, 1.43, 1.44, 1.78, 0.80, 0.80], c='magenta')
                    '''
                    #60 kpc limit
                    plt.plot([25.87, 17.54, 9.76, 354.26-360], [4.03, 24.01, 22.32, 23.51], c='r')
                    #plt.plot([34.58, 26.90, 5.80], [18.71, 23.22, 2.85], c='teal')
                    plt.plot([37.43, 30.59, 27.48, 21.28, 14.58, 8.60], [4.65, 4.72, 5.18, 4.72, 3.84, 4.59], c='b')
                    plt.plot([36.89, 33.17, 29.34, 26.06, 22.8, 14.66], [0.95, 0.85, 1.21, 1.27, 0.91, 0.85], c='magenta')
                    plt.plot([18.93, 16.44, 14.23], [1.51, 1.40, 0.22], c='orange')
                plt.scatter(ra, sigma, marker='o', c='navy', s=20)
                plt.annotate(str(wedge) + '.' + str(j+1), (ra, sigma), color='k', fontsize=10)
            else:
                print("error <plotStreamData()>: mode must be 'radec', 'radist', or 'rasigma'")

def printData(filename, stream_data, star_data, mode="radec", wrapra=False):

    print('---------------------')
    print('Printing Data Nicely')
    print('---------------------')

    f = open(filename, 'w')
    fig = plt.figure(figsize=(figwidth, figheight))

    for i in range(len(stripes)):
        wedge = stripes[i]
        stripe = stream_data[i]

        f.write('Stripe ' + str(wedge)+'\n')
        f.write('Stream\tStars\tR.A. (deg)\t\tDec. (deg)\t\tR (kpc)\t\t\tsigma (kpc)\n')
        f.write('-------------------------------------------------------\n')

        epsilon = []
        mu = []
        r = []
        theta = []
        phi = []
        sigma = []
        ra = []
        dec = []
        for j in range(len(stripe)):
            epsilon.append(stripe[j][0])
            mu = stripe[j][1]
            r.append(stripe[j][2])
            theta.append(stripe[j][3])
            phi.append(stripe[j][4])
            sigma.append(stripe[j][5])

            ra_tmp, dec_tmp = GCToEq(mu, 0, wedge, wrapra=True)
            ra.append(ra_tmp)
            dec.append(dec_tmp)

        #convert epsilon into star counts
        total_star_counts = len(np.transpose(star_data[i])[0])
        star_counts = []

        divisor = 1
        for j in range(len(epsilon)):
            divisor += np.exp(epsilon[j])

        star_counts = [int(np.exp(e)/divisor*total_star_counts) for e in epsilon]

        for j in range(len(star_counts)):
            f.write(str(wedge)+'.'+str(j+1)+'\t'+str(star_counts[j])+'\t'+str(ra[j])+'\t'+str(dec[j])+'\t'+str(r[j])+'\t'+str(sigma[j])+'\n')

        f.write('-------------------------------------------------------\n')

        plt.scatter(ra, star_counts, marker='o', c='navy', s=20)

        for j in range(len(ra)):
            plt.annotate(str(wedge) + '.' + str(j+1), (ra[j], star_counts[j]), color='k', fontsize=10)

    print('---------------------')
    print('Plotting Ra/StarCount Plot')
    print('---------------------')
    '''
    plt.plot([18.67, 16.47, 17.12, 16.00], [1106, 1484, 501, 448], c='r')
    plt.plot([25.39, 24.15, 27.51, 23.14, 15.79, 15.87], [6412, 17184, 16967, 15261, 13923, 18727], c='b')
    plt.plot([34.87, 33.98, 39.10], [8715, 5596, 7800], c='g')
    plt.plot([35.09, 37.04, 32.70, 29.31, 26.87, 21.90, 22.01], [1246, 3360, 4690, 3285, 3052, 2171, 2921], c='magenta')
    '''
    #60 kpc limit
    plt.plot([25.87, 17.54, 9.76, 354.26-360], [6880, 9339, 7582, 4865], c='r')
    #plt.plot([34.58, 26.90, 5.80], [5371, 7162, 1298], c='teal')
    plt.plot([37.43, 30.59, 27.48, 21.28, 14.58, 8.60], [12228, 16902, 18184, 14621, 14086, 21055], c='b')
    plt.plot([36.89, 33.17, 29.34, 26.06, 22.8, 14.66], [2949, 2602, 2713, 3112, 1959, 592], c='magenta')
    plt.plot([18.93, 16.44, 14.23], [943, 784, 20], c='orange')

    plt.title('Separation Application: South Redo 3')
    plt.xlabel('RA')
    plt.ylabel('Star Counts')
    plt.xlim((45, -15))
    plt.ylim((0, 22000))

    plt.savefig("/home/donlot/Desktop/separation/figures/South_Redo_3_ra_star_counts")


#=======================================
#RUNTIME
#=======================================

#read in the star data
star_data = [] #array of the data for every stripe
"""
[ <overall list>
    [ <one stripe>
        [ <one star>
            mu, nu, r (at least I think it's mu, nu, r, might be wrong but I'll check that later)
        ] ...
    ] ...
]
"""

for file in star_files:
    star_data.append(readStarFile(file))
star_data = np.array(star_data)


#read in the optimized stream data
stream_data = [] #will probably change the data structure here, but for now it is an array
"""
[ <overall list>
    [ <one stripe>
        [ <one stream center>
            epsilon, mu, r, theta, phi, sigma (need to figure out what these parameters are)
            epsilon = stream weight (fraction of stars in the wedge that belong to the stream)
            f = exp(epsilon)/(1+sum(exp(epsilon_j) over j=1 to j=k for k streams))

            mu = stream location in sdss GC coords
            r = stream distance
            theta = stream position/orientation
            phi = stream position/orientation
            sigma = width parameter

        ] ...
    ] ...
]
"""

for file in stream_files:
    stream_data.append(readStreamFile(file))
stream_data = np.array(stream_data)


#plot the star data for ra/dec
fig = plt.figure(figsize=(figwidth, figheight)) #generate the figure for the plotting routines to use

print('---------------------')
print('Plotting Ra/Dec Plot')
print('---------------------')

plotStarData(star_data, mode="radec", wrapra=True)

#plot the optimized stream data
plotStreamData(stream_data, mode="radec", wrapra=True)

#any extra functionality that we want within the separation application plotter
plt.title('Separation Application: South Redo 3')
plt.xlabel('RA')
plt.ylabel('Dec')
plt.xlim((45, -15))
plt.ylim((-12.5, 7.5))

#print out plot
plt.savefig("/home/donlot/Desktop/separation/figures/South_Redo_3_ra_dec")

#---------------------------------------

#plot the star data for ra/dec
fig = plt.figure(figsize=(figwidth, figheight)) #generate the figure for the plotting routines to use

print('---------------------')
print('Plotting Ra/Dist Plot')
print('---------------------')

plotStarData(star_data, mode="radist", wrapra=True)

#plot the optimized stream data
plotStreamData(stream_data, mode="radist", wrapra=True)

#any extra functionality that we want within the separation application plotter
plt.title('Separation Application: South Redo 3')
plt.xlabel('RA')
plt.ylabel('Dist')
plt.xlim((45, -15))
plt.ylim((0,60))

#print out plot
plt.savefig("/home/donlot/Desktop/separation/figures/South_Redo_3_ra_dist")

#---------------------------------------

fig = plt.figure(figsize=(figwidth, figheight)) #generate the figure for the plotting routines to use

print('---------------------')
print('Plotting Ra/Dec Plot')
print('---------------------')

plotStarData(star_data, mode="radec", wrapra=True)

#plot the optimized stream data
plotStreamData(stream_data, mode="radec", wrapra=True, outline_streams=True)

#any extra functionality that we want within the separation application plotter
plt.title('Separation Application: South Redo 3')
plt.xlabel('RA')
plt.ylabel('Dec')
plt.xlim((45, -15))
plt.ylim((-12.5, 7.5))

#print out plot
plt.savefig("/home/donlot/Desktop/separation/figures/South_Redo_3_ra_dec_with_streams")

#---------------------------------------

#plot the star data for ra/dec
fig = plt.figure(figsize=(figwidth, figheight)) #generate the figure for the plotting routines to use

print('---------------------')
print('Plotting Ra/Dist Plot')
print('---------------------')

plotStarData(star_data, mode="radist", wrapra=True)

#plot the optimized stream data
plotStreamData(stream_data, mode="radist", wrapra=True, outline_streams=True)

#any extra functionality that we want within the separation application plotter
plt.title('Separation Application: South Redo 3')
plt.xlabel('RA')
plt.ylabel('Dist')
plt.xlim((45, -15))
plt.ylim((0,60))

#print out plot
plt.savefig("/home/donlot/Desktop/separation/figures/South_Redo_3_ra_dist_with_streams")

fig = plt.figure(figsize=(figwidth, figheight)) #generate the figure for the plotting routines to use

print('---------------------')
print('Plotting Ra/Sigma Plot')
print('---------------------')

#plot the optimized stream data
plotStreamData(stream_data, mode="rasigma", wrapra=True, outline_streams=True)

#any extra functionality that we want within the separation application plotter
plt.title('Separation Application: South Redo 3')
plt.xlabel('RA')
plt.ylabel('Sigma')
plt.xlim((45, -15))
plt.ylim((0, 30))

#print out plot
plt.savefig("/home/donlot/Desktop/separation/figures/South_Redo_3_ra_sigma_with_streams")

printData('/home/donlot/Desktop/separation/separation_data/South_Redo_3/easy_to_read.txt', stream_data, star_data)
