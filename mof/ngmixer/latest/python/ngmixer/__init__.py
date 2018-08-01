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

from . import forcephot

#######################################################
# setup fitter dict
# each time you create a new fitter class, add it to this dict
FITTERS = {}

# boot fitters
from . import bootfit
FITTERS['max-ngmix-boot'] = bootfit.MaxNGMixBootFitter
FITTERS['metacal-ngmix-boot'] = bootfit.MetacalNGMixBootFitter
FITTERS['metacal-admom-boot'] = bootfit.MetacalAdmomBootFitter
FITTERS['isamp-ngmix-boot'] = bootfit.ISampNGMixBootFitter

# deprecated
from . import deconvolvers
FITTERS['deconv'] = deconvolvers.Deconvolver
FITTERS['metacal-deconv'] = deconvolvers.MetacalDeconvolver
FITTERS['metacal-deconv-psfbase'] = \
        deconvolvers.MetacalDeconvolverPSFBase

# forced photometry fitter
from . import forcephot
FITTERS['galsim-forcephot'] = forcephot.ForcedPhotometryFitter


try:
    from .githash import hash as __gitrepohash__
except:
    __gitrepohash__ = None
