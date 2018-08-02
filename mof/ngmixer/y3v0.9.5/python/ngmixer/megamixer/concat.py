from __future__ import print_function
import os
import sys
import numpy
import fitsio
from .. import files
from ..util import Namer

import esutil as eu

class ConcatError(Exception):
    def __init__(self, value):
        self.value = value
    def __str__(self):
        return repr(self.value)

class Concat(object):
    """
    Concatenate split files
    """
    def __init__(self,
                 config_file,
                 chunk_list,
                 output_file,
                 clobber=False,
                 skip_errors=False,
                 **kwargs):

        self.config_file=config_file

        self.output_file = output_file
        self.clobber = clobber
        self.skip_errors = skip_errors
        self.chunk_file_list = chunk_list

        self.config = files.read_yaml(config_file)

        self.make_collated_dir()
        self.tmpdir = files.get_temp_dir() 

    def pick_fields(self,data0):
        """
        modify and/or add fields to the data array
        """
        return data0

    def pick_epoch_fields(self,epoch_data0):
        """
        modify and/or add fields to the epoch data array
        """
        return epoch_data0

    def pick_nbrs_fields(self,nbrs_data0):
        """
        modify and/or add fields to the nbrs data array
        """
        return nbrs_data0

    def make_collated_dir(self):
        """
        set collated file output dir
        """
        files.makedirs_fromfile(self.output_file)


    def read_chunk(self, fname):
        """
        Read the chunk data
        """
        if not os.path.exists(fname):
            raise ConcatError("file not found: %s" % fname)

        try:
            with fitsio.FITS(fname) as fobj:
                data0 = fobj['model_fits'][:]
                epoch_data0 = fobj['epoch_data'][:]
                if 'nbrs_data' in fobj:
                    nbrs_data0 = fobj['nbrs_data'][:]
                else:
                    nbrs_data0 = None
                meta = fobj['meta_data'][:]
        except IOError as err:
            raise ConcatError(str(err))

        if 'magzp_ref' not in meta.dtype.names:
            meta = self._add_magzp_ref(meta)

        data = self.pick_fields(data0,meta)
        epoch_data = self.pick_epoch_fields(epoch_data0)
        if nbrs_data0 is not None:
            nbrs_data = self.pick_nbrs_fields(nbrs_data0)
        else:
            nbrs_data = None

        return data, epoch_data, nbrs_data, meta

    def _add_magzp_ref(self, meta):
        if 'magzp_ref' not in self.config:
            raise ValueError("You must have magzp_ref in the "
                             "meta_data extension or in the config file")
        #print("adding magzp_ref from config file")
        add_dt=[('magzp_ref','f8')]
        new_meta = eu.numpy_util.add_fields(meta, add_dt)
        new_meta['magzp_ref'] = self.config['magzp_ref']

        return new_meta

    def verify(self):
        """
        just run through and read the data, verifying we can read it
        """

        dlist = []

        nchunk=len(self.chunk_file_list)
        for i,fname in enumerate(self.chunk_file_list):

            print('\t%d/%d %s' %(i+1,nchunk,fname))
            sys.stdout.flush()
            try:
                data, epoch_data, nbrs_data, meta = self.read_chunk(fname)
                dlist.extend(list(data))
            except ConcatError as err:
                print("error found: %s" % str(err))

        if len(dlist) == 0:
            print("found no data! could not make or verify collated file %s!" % self.output_file)
        else:
            data = numpy.array(dlist,dtype=data.dtype.descr)
            if not numpy.array_equal(numpy.sort(numpy.unique(data['number'])),numpy.sort(data['number'])):
                print("object 'number' field is not unique!")

    def concat(self):
        """
        actually concatenate the data, and add any new fields
        """
        print('writing:',self.output_file)

        if os.path.exists(self.output_file) and not self.clobber:
            print('file already exists, skipping')
            return

        dlist=[]
        elist=[]
        nlist=[]
        ndtype = None
        
        nchunk=len(self.chunk_file_list)
        for i,fname in enumerate(self.chunk_file_list):

            print('\t%d/%d %s' %(i+1,nchunk,fname))
            sys.stdout.flush()
            try:
                data, epoch_data, nbrs_data, meta = self.read_chunk(fname)
                dlist.extend(list(data))
                elist.extend(list(epoch_data))
                if nbrs_data is not None:
                    nlist.extend(list(nbrs_data))
                    ndtype = nbrs_data.dtype.descr
            except ConcatError as err:
                if not self.skip_errors:
                    raise err
                print("\tskipping problematic chunk")

        if len(dlist) == 0 or len(elist)==0:
            print("\tNo good chunks found, skipping entire data set")
            return

        data = numpy.array(
            dlist,
            dtype=dlist[0].dtype.descr,
        )
        epoch_data = numpy.array(
            elist,
            dtype=elist[0].dtype.descr,
        )
        if len(nlist) > 0:
            nbrs_data = numpy.array(nlist,dtype=ndtype)
        else:
            nbrs_data = None

        if not numpy.array_equal(numpy.sort(numpy.unique(data['number'])),numpy.sort(data['number'])):
            print("object 'number' field is not unique!")

        # note using meta from last file
        self._write_data(data, epoch_data, nbrs_data, meta)

    def _write_data(self, data, epoch_data, nbrs_data, meta):
        """
        write the data, first to a local file then staging out
        the the final location
        """
        with files.StagedOutFile(self.output_file, tmpdir=self.tmpdir) as sf:
            with fitsio.FITS(sf.path,'rw',clobber=True) as fits:
                fits.write(data,extname="model_fits")
                fits.write(epoch_data,extname='epoch_data')
                if nbrs_data is not None:
                    fits.write(nbrs_data,extname='nbrs_data')
                fits.write(meta,extname="meta_data")

        print('output is in:',self.output_file)
