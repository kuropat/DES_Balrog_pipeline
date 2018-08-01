from __future__ import print_function
import os
import sys
import numpy
import fitsio
from .. import files
from ..util import Namer

from .concat import Concat,ConcatError
from .desconcat import DESConcat

class FPConcat(DESConcat):
    def __init__(self,*args,**kwargs):
        assert 'bands' in kwargs,"band names must be supplied to FPConcat"
        self.bands = kwargs.pop('bands')
        self.nbands = len(self.bands)

        super(DESConcat,self).__init__(*args,**kwargs)

    def read_chunk(self, fname):
        """
        Read the chunk data
        """
        return super(DESConcat,self).read_chunk(fname)

    def pick_fields(self, data0, meta):
        """
        pick out some fields, add some fields, rename some fields
        """
        import esutil as eu

        nbands=self.nbands

        names=list( data0.dtype.names )
        dt=[tdt for tdt in data0.dtype.descr]

        n=Namer('fp')

        # make sure the mag is associated with the fp values
        ind = names.index( n('flux_s2n') )
        if ind == len(names)-1:
            dt.append( (n('mag'),'f8',nbands) )
        else:
            dt.insert(ind, (n('mag'),'f8',nbands) )

        data=numpy.zeros(data0.size, dtype=dt)
        eu.numpy_util.copy_fields(data0, data)

        for band in xrange(nbands):
            self.calc_mag(data, meta, band)

        return data

    def calc_mag(self, data, meta, band):

        n=Namer('fp')

        nband=self.nbands
        if nband == 1:
            data[n('mag')] = -9999.
            flux = data[n('flux')].clip(min=0.001)

        else:
            data[n('mag')][:,band] = -9999.
            flux = data[n('flux')][:,band].clip(min=0.001)

        magzero=meta['magzp_ref'][band]

        if nband==1:
            data[n('mag')] = magzero - 2.5*numpy.log10( flux )
        else:
            data[n('mag')][:,band] = magzero - 2.5*numpy.log10( flux )


