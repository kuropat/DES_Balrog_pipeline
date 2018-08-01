#!/usr/bin/env python
from __future__ import print_function
import os
import sys
import copy
import numpy
import fitsio
import meds

from .pbar import prange

def get_dummy_fofs(numbers):
    nobj = len(numbers)
    fofs = numpy.zeros(nobj,dtype=[('fofid','i8'),('number','i8')])
    fofs['fofid'][:] = numpy.arange(nobj)
    fofs['number'][:] = numbers[:] #subscript should make a copy
    return fofs

class MedsNbrs(object):
    """
    Gets nbrs of any postage stamp in the MEDS.

    A nbr is defined as any stamp which overlaps the stamp under consideration
    given a buffer or is in the seg map. See the code below.

    Options:
        buff_type - how to compute buffer length for stamp overlap
            'min': minimum of two stamps
            'max': max of two stamps
            'tot': sum of two stamps

        buff_frac - fraction by whch to multiply the buffer

        maxsize_to_replace - postage stamp size to replace with maxsize
        maxsize - size ot use instead of maxsize_to_replace to compute overlap

        check_seg - use object's seg map to get nbrs in addition to postage stamp overlap
    """

    def __init__(self,meds_list,conf):

        if isinstance(meds_list, meds.MEDS):
            meds_list=[meds_list]

        self.meds_list = meds_list
        self.conf = conf

        self._init_bounds()

    def _init_bounds(self):
        self.l = {}
        self.r = {}
        self.t = {}
        self.b = {}
        self.sze = {}

        for band,m in enumerate(self.meds_list):
            #expand the stamps and get edges
            dsize = (self.conf['new_maxsize']-self.conf['maxsize_to_replace'])/2
            self.sze[band] = m['box_size'].copy()
            self.l[band] = m['orig_start_row'][:,0].copy()
            self.r[band] = m['orig_start_row'][:,0].copy()
            self.b[band] = m['orig_start_col'][:,0].copy()
            self.t[band] = m['orig_start_col'][:,0].copy()

            q, = numpy.where(self.sze[band] == self.conf['maxsize_to_replace'])
            if q.size > 0:
                self.sze[band][q[:]] = self.conf['new_maxsize']
                self.l[band][q[:]] -= dsize
                self.b[band][q[:]] -= dsize

            self.r[band] += self.sze[band]
            self.t[band] += self.sze[band]

    def get_nbrs(self,verbose=True):
        #data types
        nbrs_data = []
        dtype = [('number','i8'),('nbr_number','i8')]
        print("config:",self.conf)

        for mindex in prange(self.meds_list[0].size):
            nbrs = []
            for band,m in enumerate(self.meds_list):
                #make sure MEDS lists have the same objects!
                assert m['number'][mindex] == self.meds_list[0]['number'][mindex]
                assert m['id'][mindex] == self.meds_list[0]['id'][mindex]
                assert m['number'][mindex] == mindex+1

                #add on the nbrs
                nbrs.extend(list(self.check_mindex(mindex,band)))

            #only keep unique nbrs
            nbrs = numpy.unique(numpy.array(nbrs))

            #add to final list
            for nbr in nbrs:
                nbrs_data.append((m['number'][mindex],nbr))

        #return array sorted by number
        nbrs_data = numpy.array(nbrs_data,dtype=dtype)
        i = numpy.argsort(nbrs_data['number'])
        nbrs_data = nbrs_data[i]

        return nbrs_data

    def check_mindex(self,mindex,band):
        m = self.meds_list[band]

        #check that current gal has OK stamp, or return bad crap
        if m['orig_start_row'][mindex,0] == -9999 or m['orig_start_col'][mindex,0] == -9999:
            nbr_numbers = numpy.array([-1],dtype=int)
            return nbr_numbers

        #get the nbrs from two sources
        # 1) intersection of postage stamps
        # 2) seg map vals
        nbr_numbers = []

        #box intersection test and exclude yourself
        #use buffer of 1/4 of smaller of pair of stamps
        buff = self.sze[band].copy()
        if self.conf['buff_type'] == 'min':
            q, = numpy.where(buff[mindex] < buff)
            if len(q) > 0:
                buff[q[:]] = buff[mindex]
        elif self.conf['buff_type'] == 'max':
            q, = numpy.where(buff[mindex] > buff)
            if len(q) > 0:
                buff[q[:]] = buff[mindex]
        elif self.conf['buff_type'] == 'tot':
            buff = buff[mindex] + buff
        else:
            assert False, "buff_type '%s' not supported!" % self.conf['buff_type']
        buff = buff*self.conf['buff_frac']
        q, = numpy.where((~((self.l[band][mindex] > self.r[band]-buff) | (self.r[band][mindex] < self.l[band]+buff) |
                            (self.t[band][mindex] < self.b[band]+buff) | (self.b[band][mindex] > self.t[band]-buff))) &
                         (m['number'][mindex] != m['number']) &
                         (m['orig_start_row'][:,0] != -9999) & (m['orig_start_col'][:,0] != -9999))

        if len(q) > 0:
            nbr_numbers.extend(list(m['number'][q]))

        #check coadd seg maps
        if self.conf['check_seg']:
            try:
                segmap = m.get_cutout(mindex,0,type='seg')
                q = numpy.where((segmap > 0) & (segmap != m['number'][mindex]))
                if len(q) > 0:
                    nbr_numbers.extend(list(numpy.unique(segmap[q])))
            except:
                pass

        #cut weird crap
        if len(nbr_numbers) > 0:
            nbr_numbers = numpy.array(nbr_numbers,dtype=int)
            nbr_numbers = numpy.unique(nbr_numbers)
            inds = nbr_numbers-1
            q, = numpy.where((m['orig_start_row'][inds,0] != -9999) & (m['orig_start_col'][inds,0] != -9999))
            if len(q) > 0:
                nbr_numbers = list(nbr_numbers[q])
            else:
                nbr_numbers = []

        #if have stuff return unique else return -1
        if len(nbr_numbers) == 0:
            nbr_numbers = numpy.array([-1],dtype=int)
        else:
            nbr_numbers = numpy.array(nbr_numbers,dtype=int)
            nbr_numbers = numpy.unique(nbr_numbers)

        return nbr_numbers

class NbrsFoF(object):
    def __init__(self,nbrs_data):
        self.nbrs_data = nbrs_data
        self.Nobj = len(numpy.unique(nbrs_data['number']))

        #records fofid of entry
        self.linked = numpy.zeros(self.Nobj,dtype='i8')
        self.fofs = {}

        self._fof_data = None

    def get_fofs(self,verbose=True):
        self._make_fofs(verbose=verbose)
        return self._fof_data

    def _make_fofs(self,verbose=True):
        #init
        self._init_fofs()


        for i in prange(self.Nobj):
            self._link_fof(i)

        for fofid,k in enumerate(self.fofs):
            inds = numpy.array(list(self.fofs[k]),dtype=int)
            self.linked[inds[:]] = fofid
        self.fofs = {}

        self._make_fof_data()

    def _link_fof(self,mind):
        #get nbrs for this object
        nbrs = set(self._get_nbrs_index(mind))

        #always make a base fof
        if self.linked[mind] == -1:
            fofid = copy.copy(mind)
            self.fofs[fofid] = set([mind])
            self.linked[mind] = fofid
        else:
            fofid = copy.copy(self.linked[mind])

        #loop through nbrs
        for nbr in nbrs:
            if self.linked[nbr] == -1 or self.linked[nbr] == fofid:
                #not linked so add to current
                self.fofs[fofid].add(nbr)
                self.linked[nbr] = fofid
            else:
                #join!
                self.fofs[self.linked[nbr]] |= self.fofs[fofid]
                del self.fofs[fofid]
                fofid = copy.copy(self.linked[nbr])
                inds = numpy.array(list(self.fofs[fofid]),dtype=int)
                self.linked[inds[:]] = fofid

    def _make_fof_data(self):
        self._fof_data = []
        for i in xrange(self.Nobj):
            self._fof_data.append((self.linked[i],i+1))
        self._fof_data = numpy.array(self._fof_data,dtype=[('fofid','i8'),('number','i8')])
        i = numpy.argsort(self._fof_data['number'])
        self._fof_data = self._fof_data[i]
        assert numpy.all(self._fof_data['fofid'] >= 0)

    def _init_fofs(self):
        self.linked[:] = -1
        self.fofs = {}

    def _get_nbrs_index(self,mind):
        q, = numpy.where((self.nbrs_data['number'] == mind+1) & (self.nbrs_data['nbr_number'] > 0))
        if len(q) > 0:
            return list(self.nbrs_data['nbr_number'][q]-1)
        else:
            return []

class NbrsFoFExtractor(object):
    """
    Class to extract subet set of FoF file and destroy on exit if wanted.
    """

    def __init__(self, fof_file, start, end, sub_file, cleanup=False):
        self.fof_file = fof_file
        self.start = start
        self.end = end
        self.sub_file = sub_file
        self.cleanup = cleanup
        self._check_inputs()

        self._extract()

    def __enter__(self):
        return self

    def __exit__(self, exception_type, exception_value, traceback):
        self.close()

    def __del__(self):
        self.close()

    def close(self):
        if self.cleanup:
            if os.path.exists(self.sub_file):
                print('removing sub file:',self.sub_file)
                os.remove(self.sub_file)

    def _get_inds(self, data):
        inds = []
        for fofid in range(self.start,self.end+1):
            q, = numpy.where(data['fofid'] == fofid)
            if len(q) > 0:
                inds.extend(list(q))
        inds = numpy.array(inds,dtype=int)
        #always write this sorted!
        q = numpy.argsort(data['number'][inds])
        inds = inds[q]
        self.numbers = data['number'][inds]
        return inds

    def _extract(self):

        with fitsio.FITS(self.fof_file) as infits:
            print('opening sub file:',self.sub_file)
            with fitsio.FITS(self.sub_file,'rw',clobber=True) as outfits:
                old_data = infits[1][:]
                inds = self._get_inds(old_data)
                obj_data = old_data[inds]
                outfits.write(obj_data)

    def _check_inputs(self):
        if self.fof_file==self.sub_file:
            raise ValueError("output file name equals input")

        if self.start > self.end:
            raise ValueError("one must extract at least one object")
