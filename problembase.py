__author__ = "Erich L Foster <erichlf@gmail.com>"
__date__ = "2018-12-26"
#
#   adapted from problembase.py in nsbench originally developed by
#   Anders Logg <logg@simula.no>
#

from dolfin import *
try:
    from dolfin_adjoint import *

    dolfin.parameters['adjoint']['record_all'] = True
    adjointer = True
except:
    print('WARNING: Could not import DOLFIN-Adjoint. ' \
        + 'Adjointing will not be available.')
    adjointer = False


class ProblemBase:  # Base class for all problems.

    def __init__(self, options):

        # Store options
        self.options = options

        # Parameters must be defined by subclass
        self.mesh = None

        # time domain and time step
        self.t0 = 0
        self.t = self.t0
        self.T = options['T']  # final time
        self.k = options['k']

        # reset our discretization
        self.Nx = None
        self.Ny = None
        self.Nz = None

        self.solver = None
        self.output_location = ''

    def initial_conditions(self, W, t):
        pass

    def boundary_conditions(self, W, t):
        pass
