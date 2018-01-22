"""
code for image i/o
"""

class ImageIO(object):
    """
    abstract base class for reading images

    Basic Intro
    -----------

    To implement an image reader, you need to fill out the abstract methods below
    and make a __next__ and next method like this.

    class MyImageIO(ImageIO):
        def __init__(self,*args,**kwargs):
            super(MEDSImageIO, self).__init__(*args,**kwargs)
            self.num_fofs = # of mb obs lists to return

        def __next__(self):
            if self.fofindex >= self.num_fofs:
                raise StopIteration
            else:
                coadd_mb_obs_lists,se_mb_obs_lists = # get mb obs lists for fofindex
                self.fofindex += 1
                return coadd_mb_obs_lists,se_mb_obs_lists

        next = __next__

        ...

    Meta Data
    ---------

    The coadd_mb_obs_lists and se_mb_obs_lists need to have their meta data set with the field

        'id': unique id for object in obs
        'meta_data': numpy array with meta data (same dtype as returned by get_meta_data_dtype)
        'obj_flags': non-zero if whole object should be ignored

    Each list of obs per band needs to have the following meta data set
    
        'band_num': zero-indexed band 

    Each individual observation in the lists needs to have the following meta data fields set

        'id': unique id for object in obs
        'band_ind': id for object's cutout in the band
        'flags': non-zero if the observation should be ignored
        'meta_data': numpy array with epoch meta data (same dtype as returned by get_epoch_meta_data_dtype)
        'cutout_index': index of cutout in the band (must be unique at least within this band)

    Each psf of each observation needs to have the meta data field

        'Tguess': guess for size of PSF in arcsec

    You can use the update_meta_data method of these objects to set the data.

    Sub-files and extra data
    ------------------------

    The image io class will be passed three extra variables

        fof_range - range of FoF's (or objects if each object is a FoF) to yeild
        fof_file - a file (user defined) that possibly breaks the set of objects into subsets to
                   sent back from __next__ all at once
        mof_file - a MOF output file
        extra_data - a dict with extra data fields that can be used by the image i/o class

    The basic idea here is the user can define subsets of objects to all be dealt with by the calling class at once.

    This can enable load balancing or other more complicated fitting patterns.

    Note that the FoF file name is just pass to the image i/o class, so there is no restriction on
    the format or tags.

    Nbrs Modeling
    -------------
    To enable the nbrs modeling, each object has to mark which nbrs must be modeled.

    Each mb_obs_list in mb_obs_lists needs five meta data fields that are used to do this

        nbrs_inds: index into mb_obs_lists of the nbr obs
        nbrs_psfs: psf obs for each nbr
        nbrs_jacs: jaobians for each nbr
        nbrs_flags: only render nbrs with flags == 0
        cen_ind: index into mb_obs_lists of the central object
        nbrs_ids: id of the nbr - useful for double checking things

    """

    def __init__(self,*args,**kwargs):
        self.set_fof_start(0)
        self.num_fofs = 0

        if 'fof_range' in kwargs:
            self.fof_range = kwargs['fof_range']
        else:
            self.fof_range = None

        if 'fof_file' in kwargs:
            self.fof_file = kwargs['fof_file']
        else:
            self.fof_file = None

        self.mof_file=kwargs.get('mof_file',None)

        if 'extra_data' in kwargs:
            self.extra_data = kwargs['extra_data']
        else:
            self.extra_data = {}

    def get_file_meta_data(self):
        """
        returns file meta data

        if there is none, return None
        """
        raise NotImplementedError("get_file_meta_data method of ImageIO must be defined in subclass.")

    def get_num_bands(self):
        """
        returns number of bands for galaxy images
        """
        raise NotImplementedError("get_num_bands method of ImageIO must be defined in subclass.")

    def get_meta_data_dtype(self):
        """
        returns a numpy dtype for the galaxy meta data as a list

        For example

            return [('id','i8'),('number','i4')]
        """
        raise NotImplementedError("get_meta_data_dtype method of ImageIO must be defined in subclass.")

    def get_epoch_meta_data_dtype(self):
        """
        returns a numpy dtype for the galaxy per epoch meta data as a list

        For example

            return [('id','i8'),('number','i4'),('icut','i4')]
        """
        raise NotImplementedError("get_epoch_meta_data_dtype method of ImageIO must be defined in subclass.")

    def set_fof_start(self,start):
        self.fof_start = start

    def get_num_fofs(self):
        """
        returns # of fofs code will yeild
        """
        raise NotImplementedError("get_num_fofs method of ImageIO must be defined in subclass.")

    def __iter__(self):
        self.fofindex = self.fof_start
        return self
