from __future__ import print_function
import os
import sys
import numpy
import fitsio
from .. import files
from ..util import Namer

from .concat import Concat,ConcatError

SHAPENOISE=0.20
SHAPENOISE2=SHAPENOISE**2

class DESConcat(Concat):
    """
    Concat for DES, database prep and blinding
    """
    def __init__(self,*args,**kwargs):
        assert 'bands' in kwargs,"band names must be supplied to DESConcat"
        self.bands = kwargs.pop('bands')
        self.nbands = len(self.bands)
        self.blind = kwargs.pop('blind',True)

        super(DESConcat,self).__init__(*args,**kwargs)
        #self.config['fit_models'] = list(self.config['model_pars'].keys())

    def read_chunk(self, fname):
        """
        Read the chunk data
        """
        d,ed,nd,m = super(DESConcat,self).read_chunk(fname)

        if self.blind:
            self.blind_data(d)

        return d,ed,nd,m

    def blind_data(self,data):
        """
        multiply all shear type values by the blinding factor

        This must be run after copying g values out of pars into
        the model_g fields

        This also includes the Q values from B&A
        """

        blinding = self.config['collate']['blinding']
        if blinding =='y3':
            from blind_des_catalog_y3 import blind_arrays
        else:
            raise ValueError("bad blinding scheme: '%s'" % blinding)

        models=self.get_models(data)

        names=data.dtype.names

        if 'mcal_g' not in names:
            raise RuntimeError("you requested blinding but "
                               "mcal_g field is not present in the data")

        w,=numpy.where(data['mcal_flags'] == 0)
        types=[
            ('mcal_pars',   2,3),
            ('mcal_g',      0,1),
        ]

        for name,i1,i2 in types:

            bg1,bg2 = blind_arrays(
                data[name][w,i1],
                data[name][w,i2],
            )
            data[name][w,i1] = bg1
            data[name][w,i2] = bg2

    def pick_epoch_fields(self, epoch_data0):
        """
        pick out some fields, add some fields, rename some fields
        """

        wkeep,=numpy.where(epoch_data0['cutout_index'] >= 0)
        if wkeep.size==0:
            print("None found with cutout_index >= 0")
            print(epoch_data0['cutout_index'])
            #return numpy.zeros(1)
            return []

        epoch_data0=epoch_data0[wkeep]

        dt=epoch_data0.dtype.descr

        names=epoch_data0.dtype.names
        ind=names.index('band_num')
        dt.insert( ind, ('band','S1') )

        epoch_data = numpy.zeros(epoch_data0.size, dtype=dt)
        for nm in epoch_data0.dtype.names:
            if nm in epoch_data.dtype.names:
                epoch_data[nm] = epoch_data0[nm]

        for band_num in xrange(self.nbands):
            w,=numpy.where(epoch_data['band_num'] == band_num)
            if w.size > 0:
                epoch_data['band'][w] = self.bands[band_num]

        return epoch_data

    def pick_nbrs_fields(self, nbrs_data0):
        """
        pick out some fields, add some fields, rename some fields
        """

        dt=nbrs_data0.dtype.descr

        names=nbrs_data0.dtype.names
        ind=names.index('band_num')
        dt.insert( ind, ('band','S1') )

        nbrs_data = numpy.zeros(nbrs_data0.size, dtype=dt)
        for nm in nbrs_data0.dtype.names:
            if nm in nbrs_data.dtype.names:
                nbrs_data[nm] = nbrs_data0[nm]

        for band_num in xrange(self.nbands):
            w,=numpy.where(nbrs_data['band_num'] == band_num)
            if w.size > 0:
                nbrs_data['band'][w] = self.bands[band_num]

        return nbrs_data

    def get_models(self, data):
        models=[]

        model_names = self.config['model_pars'].keys()
        model_names = ['%s_max' % mod for mod in model_names] + model_names
        model_names = ['coadd_%s' % mod for mod in model_names] + model_names

        names=list( data.dtype.names )
        for model in model_names:
            n=Namer(model)

            if n('flux') in names or n('pars') in names:
                models.append(model)

        return models

    def make_flux_tuple(self, name, dtype, nbands):
        if nbands==1:
            tup=(name, dtype)
        else:
            if 'cov' in name:
                tup=(name,dtype,(nbands,nbands))
            else:
                tup=(name,dtype,nbands)
        return tup

    def pick_fields(self, data0, meta):
        """
        pick out some fields, add some fields, rename some fields
        """
        import esutil as eu

        nbands=self.nbands

        models=self.get_models(data0)

        names=list( data0.dtype.names )
        dt=[tdt for tdt in data0.dtype.descr]

        if 'coadd_psf_flux_err' in names:
            flux_ind = names.index('coadd_psf_flux_err')
            dt.insert(flux_ind+1, ('coadd_psf_flux_s2n','f8',nbands) )
            names.insert(flux_ind+1,'coadd_psf_flux_s2n')

            dt.insert(flux_ind+2, ('coadd_psf_mag','f8',nbands) )
            names.insert(flux_ind+2,'coadd_psf_mag')


        tmpind = names.index('psf_flux_s2n')
        dt.insert(tmpind, ('psf_mag','f8',nbands) )
        names.insert(tmpind,'psf_mag')

        T_namers=[]
        do_flux=False
        for ft in models:

            n=Namer(ft)

            if '%s_g_cov' % ft in names:
                gcind = names.index('%s_g_cov' % ft)
                wtf = (n('weight'), 'f8')
                dt.insert(gcind+1, wtf)
                names.insert(gcind+1, n('weight'))

            if n('flux') in names:
                # we can just copy it over
                flux_ind = names.index(n('flux'))
            else:
                # we need to get flux from the pars vector
                do_flux = True
                pars_cov_ind = names.index(n('pars_cov'))

                offset = 1

                for fi,name in enumerate(['flux','flux_cov']):
                    tup=self.make_flux_tuple(n(name),'f8',nbands)
                    dt.insert(pars_cov_ind+offset, tup)
                    names.insert(pars_cov_ind+offset,n(name))

                flux_ind = names.index(n('flux'))

            offset=1

            for fi,name in enumerate(['flux_s2n','mag','logsb']):
                tup=self.make_flux_tuple(n(name),'f8',nbands)

                dt.insert(flux_ind+offset, tup)
                names.insert(flux_ind+offset,n(name))

                offset += 1

            if n('T') not in data0.dtype.names:
                T_namers.append(n)
                fadd=[(n('T'),'f8'),
                      (n('T_err'),'f8'),
                      (n('T_s2n'),'f8')]
                ind = names.index('%s_pars_cov' % ft)
                for f in fadd:
                    dt.insert(ind+1, f)
                    names.insert(ind+1, f[0])
                    ind += 1

        data=numpy.zeros(data0.size, dtype=dt)
        eu.numpy_util.copy_fields(data0, data)


        all_models=models
        if 'coadd_psf_flux' in names:
            all_models=all_models + ['coadd_psf']
        if 'psf_flux' in names:
            all_models=all_models + ['psf']

        if len(T_namers) > 0:
            self.add_T_info(T_namers, data)

        for ft in all_models:
            for band in xrange(nbands):
                self.calc_mag_and_flux_stuff(
                    data, meta, ft, band,
                    do_flux=do_flux,
                )

        #self.add_weight(data, models)

        return data

    def add_T_info(self, T_namers, data):
        """
        Add T S/N etc.
        """
        for n in T_namers:

            if n('T_s2n') in data.dtype.names:

                data[n('T')][:]   = -9999.0
                data[n('T_err')][:]  =  9999.0
                data[n('T_s2n')][:] = -9999.0

                Tcov=data[n('pars_cov')][:,4,4]
                w,=numpy.where( (data[n('flags')] == 0) & (Tcov > 0.0) )
                if w.size > 0:
                    data[n('T')][w]   = data[n('pars')][w, 4]
                    data[n('T_err')][w]  =  numpy.sqrt(Tcov[w])
                    data[n('T_s2n')][w] = data[n('T')][w]/data[n('T_err')][w]


    def add_weight(self, data, models):
        """
        Add weight for each model
        """
        for model in models:
            n=Namer(model)

            w,=numpy.where(data[n('flags')]==0)
            if w.size > 0:
                if n('g_cov') in data.dtype.names:
                    c=data[n('g_cov')]
                    Csum = c[w,0,0] + c[w,1,1]
                    weight=1.0/(2.*SHAPENOISE2 + Csum)
                    data[n('weight')][w] = weight

    def calc_mag_and_flux_stuff(self, data, meta, model, band, do_flux=False):
        """
        Get magnitudes

        if do_flux, we will get flux from pars
        """

        do_flux_not_psf=(do_flux and 'psf' not in model)

        names = data.dtype.names

        n=Namer(model)

        nband=self.nbands

        if nband == 1:
            data[n('mag')] = -9999.
            data[n('flux_s2n')] = 0.0
            if do_flux_not_psf:
                data[n('flux')] = -9999.
                data[n('flux_cov')] = -9999.
        else:
            data[n('mag')][:,band] = -9999.
            data[n('flux_s2n')][:,band] = 0.0
            if do_flux_not_psf:
                data[n('flux')][:,band] = -9999.
                data[n('flux_cov')][:,:,band] = -9999.

        if 'psf' not in model:
            if nband == 1:
                data[n('logsb')] = -9999.0
            else:
                data[n('logsb')][:,band] = -9999.0

        if model in ['coadd_psf','psf']:
            if nband == 1:
                w,=numpy.where(data[n('flags')] == 0)
            else:
                w,=numpy.where(data[n('flags')][:,band] == 0)
        else:
            w,=numpy.where(data[n('flags')] == 0)

        if w.size > 0:
            if do_flux_not_psf:
                if nband == 1:
                    data[n('flux')][w] = data[n('pars')][w,5]
                    data[n('flux_cov')][w] = data[n('pars_cov')][w,5,5]
                else:
                    data[n('flux')][w,band] = data[n('pars')][w,5+band]
                    data[n('flux_cov')][w,band,band] = \
                        data[n('pars_cov')][w,5+band,5+band]

            if nband == 1:
                flux = (data[n('flux')][w]).clip(min=0.001)
            else:
                flux = (data[n('flux')][w,band]).clip(min=0.001)

            magzero=meta['magzp_ref'][band]

            if nband==1:
                data[n('mag')][w] = magzero - 2.5*numpy.log10( flux )
            else:
                data[n('mag')][w,band] = magzero - 2.5*numpy.log10( flux )

            if 'psf' not in model:
                if nband==1:
                    sb = data[n('flux')][w]/data[n('T')][w]
                    data[n('logsb')][w] = numpy.log10(numpy.abs(sb))
                else:
                    sb = data[n('flux')][w,band]/data[n('T')][w]
                    data[n('logsb')][w,band] = numpy.log10(numpy.abs(sb))

            if model in ['coadd_psf','psf']:
                if nband==1:
                    flux=data[n('flux')][w]
                    flux_err=data[n('flux_err')][w]
                    w2,=numpy.where(flux_err > 0)
                    if w2.size > 0:
                        flux=flux[w2]
                        flux_err=flux_err[w2]
                        data[n('flux_s2n')][w[w2]] = flux/flux_err
                else:
                    flux=data[n('flux')][w,band]
                    flux_err=data[n('flux_err')][w,band]
                    w2,=numpy.where(flux_err > 0)
                    if w2.size > 0:
                        flux=flux[w2]
                        flux_err=flux_err[w2]
                        data[n('flux_s2n')][w[w2],band] = flux/flux_err
            else:
                if nband==1:
                    flux=data[n('flux')][w]
                    flux_var=data[n('flux_cov')][w]

                    w2,=numpy.where(flux_var > 0)
                    if w.size > 0:
                        flux=flux[w2]
                        flux_err=numpy.sqrt(flux_var[w2])
                        data[n('flux_s2n')][w[w2]] = flux/flux_err
                else:
                    flux=data[n('flux')][w,band]
                    flux_var=data[n('flux_cov')][w,band,band]

                    w2,=numpy.where(flux_var > 0)
                    if w.size > 0:
                        flux=flux[w2]
                        flux_err=numpy.sqrt(flux_var[w2])
                        data[n('flux_s2n')][w[w2], band] = flux/flux_err

class DESAdmomMetacalConcat(DESConcat):
    def get_models(self, data):
        return ['gauss','mcal']

    def pick_fields(self, data0, meta):
        return data0

    '''
    def pick_fields(self, data0, meta):
        """
        pick out some fields, add some fields, rename some fields
        """
        import esutil as eu

        nbands=self.nbands

        models=self.get_models(data0)

        names=list( data0.dtype.names )
        dt=[tdt for tdt in data0.dtype.descr]

        tmpind = names.index('psf_flux_s2n')
        dt.insert(tmpind, ('psf_mag','f8',nbands) )
        names.insert(tmpind,'psf_mag')

        for ft in models:

            n=Namer(ft)

            gcind = names.index('%s_g_cov' % ft)
            wtf = (n('weight'), 'f8')
            dt.insert(gcind+1, wtf)
            names.insert(gcind+1, n('weight'))

            flux_ind = names.index(n('flux'))

            offset=1

            for fi,name in enumerate(['mag','logsb']):
                tup=self.make_flux_tuple(n(name),'f8',nbands)

                dt.insert(flux_ind+offset, tup)
                names.insert(flux_ind+offset,n(name))

                offset += 1

        data=numpy.zeros(data0.size, dtype=dt)
        eu.numpy_util.copy_fields(data0, data)

        all_models=models
        if 'psf_flux' in names:
            all_models=all_models + ['psf']

        for ft in all_models:
            for band in xrange(nbands):
                self.calc_mag_and_flux_stuff(
                    data, meta, ft, band,
                )

        #self.add_weight(data, models)

        return data

    def calc_mag_and_flux_stuff(self, data, meta, model, band):
        """
        Get magnitudes
        """

        names = data.dtype.names

        n=Namer(model)

        nband=self.nbands

        if nband == 1:
            data[n('mag')] = -9999.
        else:
            data[n('mag')][:,band] = -9999.

        if 'psf' not in model:
            if nband == 1:
                data[n('logsb')] = -9999.0
            else:
                data[n('logsb')][:,band] = -9999.0

        if model in ['psf']:
            if nband == 1:
                w,=numpy.where(data[n('flags')] == 0)
            else:
                w,=numpy.where(data[n('flags')][:,band] == 0)
        else:
            w,=numpy.where(data[n('flags')] == 0)

        if w.size > 0:

            if nband == 1:
                flux = (data[n('flux')][w]).clip(min=0.001)
            else:
                flux = (data[n('flux')][w,band]).clip(min=0.001)

            magzero=meta['magzp_ref'][band]

            if nband==1:
                data[n('mag')][w] = magzero - 2.5*numpy.log10( flux )
            else:
                data[n('mag')][w,band] = magzero - 2.5*numpy.log10( flux )

            if 'psf' not in model:
                if nband==1:
                    sb = data[n('flux')][w]/data[n('T')][w]
                    data[n('logsb')][w] = numpy.log10(numpy.abs(sb))
                else:
                    sb = data[n('flux')][w,band]/data[n('T')][w]
                    data[n('logsb')][w,band] = numpy.log10(numpy.abs(sb))
    '''
class DESMetacalConcat(DESConcat):
    #def __init__(self,*args,**kwargs):
    #    super(DESMetacalConcat,self).__init__(*args, **kwargs)
    #    self.config['fit_models'] += ['mcal']

    def get_models(self, data):
        models=super(DESMetacalConcat,self).get_models(data)

        # remove max models for metacal
        models = [m for m in models if 'max' not in m]
        models += ['mcal']
        return models
