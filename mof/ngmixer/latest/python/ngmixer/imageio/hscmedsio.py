from __future__ import print_function
import numpy

from .desmedsio import (
    SVDESMEDSImageIO,
    Y3DESMEDSImageIO,
)


HSC_BADPIX_MAP = {
    'BAD':                1,
    'SAT':                2,
    'INTRP':              4,
    'CR':                 8,
    'EDGE':              16,

    'DETECTED':          32, # should not mask
    'DETECTED_NEGATIVE': 64, # should not mask

    'SUSPECT':          128,
    'NO_DATA':          256,
    'BRIGHT_OBJECT':    512,
    'CLIPPED':         1024,
    'CROSSTALK':       2048,
    'NOT_DEBLENDED':   4096,
    'UNMASKEDNAN':     8192,
}

HSC_DONOTMASK = ['DETECTED','DETECTED_NEGATIVE']

class HSCDR1MedsImageIO(Y3DESMEDSImageIO):
    def _get_image_flags(self, band, mindex):
        """
        currently don't have image flags
        """
        meds=self.meds_list[band]
        ncut = meds['ncutout'][mindex]
        return numpy.zeros(ncut, dtype='i4')

    def _fill_obs_meta_data(self,obs, band, mindex, icut):
        """
        fill meta data to be included in output files
        """
        super(SVDESMEDSImageIO, self)._fill_obs_meta_data(obs, band, mindex, icut)
        meds=self.meds_list[band]
        image_id=-9999
        obs.meta['meta_data']['image_id'][0]  = image_id

    def _set_extra_mask_flags(self):
        flag_list = self.conf.get('extra_mask_flag_list',None)

        if flag_list is not None:

            flags = 0
            for flagname in flag_list:
                assert flagname not in HSC_DONOTMASK,"don't mask %s" % HSC_DONOTMASK
                flags += HSC_BADPIX_MAP[flagname]

            self.conf['extra_mask_flags'] = flags


