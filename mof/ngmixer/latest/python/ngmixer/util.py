from __future__ import print_function
import os
import numpy
from numpy import log
import fitsio
import time
import sys

import ngmix
from ngmix import srandu, GMixRangeError
from ngmix.priors import LOWVAL
from .defaults import VERBOSITY

# coordinates
# ra = -u
# ra = -phi
# v = dec
# theta = 90 - dec

# unit vector directions
# u = -ra on sphere = +phi on sphere
# v = dec on sphere = -theta on sphere

def radec_to_unitvecs_ruv(ra,dec):
    theta,phi = radec_to_thetaphi(ra,dec)
    return thetaphi_to_unitvecs_ruv(theta,phi)

def radec_to_thetaphi(ra,dec):
    return (90.0-dec)/180.0*numpy.pi,-1.0*ra/180.0*numpy.pi

def thetaphi_to_unitvecs_ruv(theta,phi):
    sint = numpy.sin(theta)
    cost = numpy.cos(theta)
    sinp = numpy.sin(phi)
    cosp = numpy.cos(phi)

    rhat = numpy.array([sint*cosp,sint*sinp,cost])
    that = numpy.array([cost*cosp,cost*sinp,-1.0*sint])
    phat = numpy.array([-1.0*sinp,cosp,0.0])

    return rhat,phat,-1.0*that


def interpolate_image(rowcen1, colcen1, jacob1, im1, 
                      rowcen2, colcen2, jacob2):
    """
    interpolate from im1 to im2 using jacobs
    
    assumes im1 and im2 are the same size
    
    returns 
        
        interpolated image, pixels in im2 off of im1, pixels in im2 on im1
    """
    
    im2 = numpy.zeros_like(im1)
        
    rows2,cols2 = numpy.mgrid[0:im2.shape[0], 0:im2.shape[1]]
    rows2 = rows2 - rowcen2
    cols2 = cols2 - colcen2

    jinv1 = jacob1.getI()
            
    # convert pixel coords in second cutout to u,v
    u = rows2*jacob2[0,0] + cols2*jacob2[0,1]
    v = rows2*jacob2[1,0] + cols2*jacob2[1,1]
    
    # now convert into pixels for first image
    row1 = rowcen1 + u*jinv1[0,0] + v*jinv1[0,1]
    col1 = colcen1 + u*jinv1[1,0] + v*jinv1[1,1]
    
    row1 = row1.round().astype('i8')
    col1 = col1.round().astype('i8')
    
    wbad = numpy.where((row1 < 0)             |
                       (row1 >= im1.shape[0]) |
                       (col1 < 0)             |
                       (col1 >= im1.shape[1]))

    wgood = numpy.where((row1 >= 0)            &
                        (row1 < im1.shape[0])  &
                        (col1 >= 0)            &
                        (col1 < im1.shape[1]))
    
    # clipping makes the notation easier
    row1 = row1.clip(0,im1.shape[0]-1)
    col1 = col1.clip(0,im1.shape[1]-1)
        
    rows2,cols2 = numpy.mgrid[0:im2.shape[0], 0:im2.shape[1]]
    
    # fill the image
    im2[rows2[wgood],cols2[wgood]] = im1[row1[wgood],col1[wgood]]

    return im2,wbad,wgood

def interpolate_image_diffsize(rowcen1, colcen1, jacob1, im1, 
                               rowcen2, colcen2, jacob2, im2):
    """
    interpolate from im1 to im2 using jacobs
    
    assumes im1 and im2 are the same size
    
    returns 
        
        interpolated image, pixels in im2 off of im1, pixels in im2 on im1
    """

    rows2abs,cols2abs = numpy.mgrid[0:im2.shape[0], 0:im2.shape[1]]
    rows2 = rows2abs - rowcen2
    cols2 = cols2abs - colcen2

    jinv1 = jacob1.getI()
            
    # convert pixel coords in second cutout to u,v
    u = rows2*jacob2[0,0] + cols2*jacob2[0,1]
    v = rows2*jacob2[1,0] + cols2*jacob2[1,1]
    
    # now convert into pixels for first image
    row1 = rowcen1 + u*jinv1[0,0] + v*jinv1[0,1]
    col1 = colcen1 + u*jinv1[1,0] + v*jinv1[1,1]
    
    row1 = row1.round().astype('i8')
    col1 = col1.round().astype('i8')

    wgood = numpy.where((row1 >= 0)            &
                        (row1 < im1.shape[0])  &
                        (col1 >= 0)            &
                        (col1 < im1.shape[1]))
    
    # clipping makes the notation easier
    row1 = row1.clip(0,im1.shape[0]-1)
    col1 = col1.clip(0,im1.shape[1]-1)
    
    # fill the image
    im2[rows2abs[wgood],cols2abs[wgood]] = im1[row1[wgood],col1[wgood]]



def print_with_verbosity(*args,**kwargs):
    """
    print with verbosity=XXX keyword
    """
    verbosity=kwargs.pop('verbosity',0)
    if verbosity <= VERBOSITY():
        print(*args,**kwargs)


def clip_element_wise(arr, minvals, maxvals):
    """
    Clip each element of an array separately
    """
    for i in xrange(arr.size):
        arr[i] = arr[i].clip(min=minvals[i],max=maxvals[i])

class UtterFailure(Exception):
    """
    could not make a good guess
    """
    def __init__(self, value):
        self.value = value
    def __str__(self):
        return repr(self.value)

class MissingDataError(Exception):
    """
    could not make a good guess
    """
    def __init__(self, value):
        self.value = value
    def __str__(self):
        return repr(self.value)


class Namer(object):
    """
    create strings with a specified front prefix
    """
    def __init__(self, front=None, back=None):
        if front=='':
            front=None
        if back=='':
            back=None

        self.front=front
        self.back=back

        if self.front is None and self.back is None:
            self.nomod=True
        else:
            self.nomod=False



    def __call__(self, name):
        if self.nomod:
            return name
        else:
            n=name
            if self.front is not None:
                n = '%s_%s' % (self.front, n)
            if self.back is not None:
                n = '%s_%s' % (n, self.back)
            return n

def print_pars(pars, fmt='%8.3g', front=None, verbosity=0):
    """
    print the parameters with a uniform width
    """
    from sys import stdout
    msg = ""
    if front is not None:
        msg += front
        msg += " "

    allfmt = ' '.join( [fmt+' ']*len(pars) )
    msg += allfmt % tuple(pars)

    print_with_verbosity(msg,verbosity=verbosity)

def print_pars_and_logl(pars, logl, fmt='%8.3g', front=None):
    """
    print the parameters with a uniform width
    """
    from sys import stdout
    msg = ""
    if front is not None:
        msg += front
        msg += " "

    allfmt = ' '.join( [fmt+' ']*len(pars) )
    msg += allfmt % tuple(pars)
    msg += " logl: %g" % logl
    print_with_verbosity(msg,verbosity=verbosity)

def plot_autocorr(trials, window=100, show=False, **kw):
    import biggles
    import emcee

    arr=biggles.FramedArray(trials.shape[1], 1)
    arr.uniform_limits=True

    func=emcee.autocorr.function(trials)
    tau2 = emcee.autocorr.integrated_time(trials, window=window)

    xvals=numpy.arange(func.shape[0])
    zc=biggles.Curve( [0,func.shape[0]-1],[0,0] )

    for i in xrange(trials.shape[1]):
        pts=biggles.Curve(xvals,func[:,i],color='blue')

        lab=biggles.PlotLabel(0.9,0.9,
                              r'$%s tau\times 2: %s$' % (i,tau2[i]),
                              halign='right')
        arr[i,0].add(pts,zc,lab)

    if show:
        arr.show(**kw)

    return arr

class CombinedImageFlags(object):
    """
    replacement astrometry flags
    """
    def __init__(self, filename):
        self.filename=filename

        self._load()

    def _load(self):
        import json
        with open(self.filename) as fobj:
            self.data=json.load(fobj)

    def get_key(self, filename):
        bname=os.path.basename(filename)
        bs=bname.split('_')

        expid = int(bs[1])
        ccdid = int( bs[2].split('.')[0] )

        key='%s-%02d' % (expid, ccdid)
        return key

    def get_flags_multi(self, image_names, default=0):

        flaglist=numpy.zeros(len(image_names))
        for i,image_name in enumerate(image_names):
            flags=self.get_flags(image_name,default=default)

            flaglist[i] = flags

        return flaglist

    def get_flags(self, image_name, default=0):
        """
        match based on image id
        """

        key=self.get_key(image_name)

        flags = self.data.get(key,default)
        return flags


class AstromFlags(object):
    """
    replacement astrometry flags
    """
    def __init__(self, filename):
        self.data=fitsio.read(filename,lower=True)

    def get_flags(self, image_ids):
        """
        match based on image id
        """
        import esutil as eu

        image_ids=numpy.array(image_ids, ndmin=1, dtype='i8',copy=False)

        # default is flagged, to indicated not found
        flags=numpy.ones(image_ids.size,dtype='i8')
        minput, mastro = eu.numpy_util.match(image_ids, self.data['imageid'])

        nmiss=image_ids.size - minput.size
        if nmiss > 0:
            print("        %d/%d did not "
                  "match astrom flags" % (nmiss,image_ids.size))
        else:
            print("    all matched")



        if minput.size > 0:
            flags[minput] = self.data['astrom_flag'][mastro]

        return flags

def seed_numpy(random_seed):
    """
    set up random number generation with the input seed
    """
    if random_seed is not None:
        numpy.random.seed(random_seed)


