#######################################################
# imports
from . import files
from . import ngmixing
from . import priors
from . import util
from . import defaults
from . import bootfit
from . import mofngmixing
from . import megamixer
from . import render_ngmix_nbrs
from .render_ngmix_nbrs import RenderNGmixNbrs, DESRenderNGmixNbrs

#######################################################
# setup fitter dict
# each time you create a new fitter class, add it to this dict
FITTERS = {}

# boot fitters
from . import bootfit
FITTERS['max-ngmix-boot'] = bootfit.MaxNGMixBootFitter
FITTERS['metacal-ngmix-boot'] = bootfit.MetacalNGMixBootFitter
FITTERS['isamp-ngmix-boot'] = bootfit.ISampNGMixBootFitter

from . import deconvolvers
FITTERS['deconv'] = deconvolvers.Deconvolver
FITTERS['metacal-deconv'] = deconvolvers.MetacalDeconvolver
FITTERS['metacal-deconv-psfbase'] = deconvolvers.MetacalDeconvolverPSFBase


try:
    from .githash import hash as __gitrepohash__
except:
    __gitrepohash__ = None
