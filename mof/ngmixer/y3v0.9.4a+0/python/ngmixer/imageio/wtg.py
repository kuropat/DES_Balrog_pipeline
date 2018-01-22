from __future__ import print_function
import os
from .desmedsio import Y3DESMEDSImageIO

class WTGMEDSImageIO(Y3DESMEDSImageIO):
    """
    Class for weighing the giants

    differs from DES Y3 in the psf map structure
    """

    def _load_psf_map(self, **kw):
        """
        we fake the coadd psf
        """
        extra_data=kw.get('extra_data',{})

        map_file=extra_data.get('psf_map',None)
        if map_file is None:
            raise RuntimeError("for Y3 you must send a map file")

        print("reading psf map:",map_file)
        psf_map={}
        with open(map_file) as fobj:
            for line in fobj:

                ls=line.strip().split()

                if len(ls) != 2:
                    raise ValueError("badly formatted psf map line: '%s'" % line)

                key, psfpath = ls

                psf_map[key] = psfpath

        self._psf_map=psf_map


    def _psfex_path_from_image_path(self, meds, image_path):
        """
        infer the psfex path from the image path.
        """

        bname = os.path.basename(image_path)
        key = bname.replace('.sub.fits','')

        psf_path = self._psf_map[key]
        psf_path = os.path.expandvars(psf_path)
        return psf_path


