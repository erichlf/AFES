__author__ = "Erich L Foster <erichlf@gmail.com>"
__date__ = "2013-08-27"
__license__ = "GNU GPL version 3 or any later version"
#
#   adapted from solverbase.py in nsbench originally developed by
#   Anders Logg <logg@simula.no>
#

from dolfin import *
try:
    from dolfin_adjoint import *

    dolfin.parameters["adjoint"]["record_all"] = True
    adjointer = True
except:
    print "WARNING: Could not import DOLFIN-Adjoint. " \
        + "Adjointing will not be available."
    adjointer = False


from time import time
from os import getpid
from commands import getoutput
import sys

# Common solver parameters
maxiter = default_maxiter = 200
tolerance = default_tolerance = 1e-4


class SolverBase:

    '''
        SolverBase provides a general solver class for the various Solvers. Its
        purpose is take a weak_residual and then solve it using the theta-method
    '''

    def __init__(self, options):

        # Set global DOLFIN parameters
        parameters['form_compiler']['cpp_optimize'] = True
        parameters['allow_extrapolation'] = True
        nonLinearSolver = NewtonSolver()
        prm = nonLinearSolver.parameters
        prm['convergence_criterion'] = 'incremental'
        prm['absolute_tolerance'] = options['absolute_tolerance']
        prm['relative_tolerance'] = options['relative_tolerance']
        prm['report'] = options['monitor_convergence']

        # Set debug level
        set_log_active(options['debug'])

        self.mem = options['check_mem_usage']

        self.saveSolution = options['save_solution']
        if self.saveSolution:
            self.plotSolution = False
        else:
            self.plotSolution = options['plot_solution']
        self.saveFrequency = options['save_frequency']

        # initialize the time stepping method parameters
        try:
            self.theta = options['theta']  # time stepping method
        except:
            self.theta = 0.5

        # Reset some solver variables
        self._time = None
        self._cputime = 0.0
        self._timestep = 0

        # adaptivity options
        self.adaptive = options['adaptive']
        self.adaptRatio = options['adapt_ratio']
        self.maxAdapts = options['max_adaptations']
        self.adaptTOL = options['adaptive_TOL']
        self.onDisk = options['on_disk']

        self.optimize = options['optimize']

        self.steady_state = False

        # set the velocity and pressure element orders
        self.Pu = options['velocity_order']
        self.Pp = options['pressure_order']

        # Reset files for storing solution
        self._ufile = None
        self._pfile = None
        self._uDualfile = None
        self._pDualfile = None
        self.eifile = None
        self.meshfile = None
        self.optfile = None

        # Reset storage for functional values and errors
        self._t = []

        # create plot objects
        self.vizU = None
        self.vizP = None

    def solve(self, problem):
        '''
            This is the general solve class which will determine if adaptivity
            should be used or if a problem is an optimization problem.
        '''
        mesh = problem.mesh
        if not self.steady_state:
            T = problem.T
            t0 = problem.t0

            # adjust time step so that we evenly divide time interval
            k = self.adjust_dt(t0, T, problem.k)
        else:
            T = None
            t0 = None

            k = None

        if self.adaptive:  # solve with adaptivity
            if adjointer:
                mesh, k = self.adaptivity(problem, mesh, T, t0, k)
            else:
                print "WARNING: You have requested adaptivity, but DOLFIN-Adjoint" \
                    + " doesn't appear to be installed."
                print "Solving without adaptivity."

        print 'Solving the primal problem.'
        self.file_naming(problem, n=-1, opt=False)

        # record so that we can evaluate our functional
        if adjointer:
            parameters["adjoint"]["stop_annotating"] = \
                not (self.adaptive or (self.optimize
                                       and 'Optimize' in dir(problem)))

        func = 'functional' in dir(problem)
        W, w, m = self.forward_solve(problem, mesh, t0, T, k, func=func)

        if m is not None:
            print
            print 'The size of the functional is: %0.3G' % m

        # solve the optimization problem
        if(self.optimize and 'Optimize' in dir(problem)):
            if adjointer:
                # give me an end line so that dolfin-adjoint doesn't
                # cover previous prints
                print
                problem.Optimize(self, W, w)

                self.file_naming(problem, n=-1, opt=True)

                parameters["adjoint"]["stop_annotating"] = True
                W, w, m = self.forward_solve(problem, mesh, t0, T, k, func=func)
            else:
                print "WARNING: You have requested Optimization, but" \
                    + " DOLFIN-Adjoint doesn't appear to be installed."
                print "Not running optimization."

        return w

    def adaptivity(self, problem, mesh, T, t0, k):
        COND = 1
        nth = ('st', 'nd', 'rd', 'th')  # numerical descriptors

        # Adaptive loop
        i = 0
        m = 0  # initialize
        while(i <= self.maxAdapts and COND > self.adaptTOL):
            # setup file names
            self.file_naming(problem, n=i, opt=False)
            # save our current mesh
            if not self.plotSolution:
                self.meshfile << mesh

            if i == 0:
                print 'Solving on initial mesh.'
            elif i < len(nth):
                print 'Solving on %d%s adapted mesh.' % (i, nth[i - 1])
            else:
                print 'Solving on %d%s adapted mesh.' % (i, nth[-1])

            # Solve primal and dual problems and compute error indicators
            m_ = m  # save the previous functional value
            W, w, m, ei = self.adaptive_solve(problem, mesh, t0, T, k)
            COND = self.condition(ei, m, m_)
            print 'DOFs=%d functional=%0.5G err_est=%0.5G' \
                % (mesh.num_vertices(), m, COND)

            if i == 0 and self.plotSolution:
                plot(ei, title="Error Indicators.", elevate=0.0)
                plot(mesh, title='Initial mesh', size=((600, 300)))
            elif (i == self.maxAdapts or COND <= self.adaptTOL) \
                    and self.plotSolution:
                plot(ei, title="Error Indicators.", elevate=0.0)
                plot(mesh, title='Final mesh', size=((600, 300)))
                interactive()
            elif not self.plotSolution:
                self.eifile << ei

            # Refine the mesh
            print 'Refining mesh.'
            mesh = self.adaptive_refine(mesh, ei)
            if 'time_step' in dir(problem) and not self.steady_state:
                k = self.adjust_dt(t0, T, problem.time_step(problem.Ubar, mesh))

            adj_reset()  # reset the dolfin-adjoint

            i += 1

        if i > self.maxAdapts and COND > self.adaptTOL:
            s = 'Warning reached max adaptive iterations with' \
                + 'sum(abs(EI))=%0.3G.  Solution may not be accurate.' \
                % COND
            print s

        return mesh, k

    def adaptive_solve(self, problem, mesh, t0, T, k):
        '''
            Adaptive solve applies the error representation to goal-oriented
            adaptivity. This is all done automatically using the weak_residual.
        '''
        print 'Solving the primal problem.'
        parameters["adjoint"]["stop_annotating"] = False

        if not self.steady_state:
            N = int(round((T - t0) / k))

            assert self.onDisk <= 1. or self.onDisk >= 0.
            if self.onDisk > 0:
                adj_checkpointing(strategy='multistage', steps=N,
                                  snaps_on_disk=int(self.onDisk * N),
                                  snaps_in_ram=int((1. - self.onDisk) * N),
                                  verbose=False)

        W, w, m = self.forward_solve(problem, mesh, t0, T, k, func=True)
        adj_html('adjoint.html', 'adjoint')
        parameters["adjoint"]["stop_annotating"] = True
        self._timestep = 0  # reset the time step to zero

        phi = Function(W)

        # Generate error indicators
        Z = FunctionSpace(mesh, "DG", 0)
        z = TestFunction(Z)
        ei = Function(Z, name='Error Indicator')
        LR1 = 0.

        # Generate the dual problem
        if self.steady_state:
            functional = problem.functional(W, w)
        else:
            functional = problem.functional(W, w) * dt
        J = Functional(functional, name='DualArgument')
        timestep = None
        wtape = []
        phi = []

        print
        print 'Solving the dual problem.'
        for (adj, var) in compute_adjoint(J, forget=False):
            # use only the last iteration or the initial condition
            if var.name == 'w' and (timestep != var.timestep
                                    or var.timestep == 0):
                timestep = var.timestep
                if var.timestep == 0 and timestep != var.timestep:
                    iteration = 1
                else:
                    iteration = 0
                # Compute error indicators ei
                wtape.append(DolfinAdjointVariable(w).
                             tape_value(timestep=timestep, iteration=iteration))
                phi.append(adj)
                self.update(problem, None, W, adj, dual=True)

        self._timestep = 0  # reset the time step to zero

        print 'Building error indicators.'

        if not self.steady_state:
            for i in range(0, len(wtape) - 1):
                # the tape is backwards so i+1 is the previous time step
                wtape_theta = self.theta * wtape[i] \
                    + (1. - self.theta) * wtape[i + 1]
                LR1 = self.weak_residual(problem, Constant(k), W, wtape_theta,
                                         wtape[i], wtape[i + 1], z * phi[i],
                                         ei_mode=True)
                ei.vector()[:] += assemble(LR1, annotate=False).array()
        else:
            LR1 = self.weak_residual(problem, W, wtape[0], z * phi[0],
                                     ei_mode=True)
            ei.vector()[:] = assemble(LR1, annotate=False).array()

        return W, w, m, ei

    def condition(self, ei, m, m_):
        '''
            Adaptive stopping criterion for non-Galerkin-orthogonal problems.
            Overload this for Galerkin-orthogonal problems.
            ei - error indicators (non-Galerkin-orthogonal problems)
            m - current functional size (Galerkin-orthogonal problems)
            m_ - previous functional size (Galerkin-orthogonal problems)
        '''
        c = abs(sum(ei.vector()))

        return c

    def forward_solve(self, problem, mesh, t0, T, k, func=False):
        '''
            Here we take the weak_residual and apply boundary conditions and
            then send it to time_stepper for solving.
        '''

        # Define function spaces
        # we do it this way so that it can be overloaded
        W = self.function_space(mesh)

        if not self.steady_state:
            ic = problem.initial_conditions(W)

        # define trial and test function
        wt = TestFunction(W)
        if adjointer:  # only use annotation if DOLFIN-Adjoint was imported
            w = Function(W, name='w')
            if not self.steady_state:
                w_ = Function(ic, name='w_')
        else:
            w = Function(W)
            if not self.steady_state:
                w_ = Function(ic)

        if not self.steady_state:
            theta = self.theta
            w_theta = (1. - theta) * w_ + theta * w

            # weak form of the primal problem
            F = self.weak_residual(problem, Constant(k), W, w_theta, w, w_, wt,
                                   ei_mode=False)

            w, m = self.timeStepper(problem, t0, T, k, W, w, w_, F, func=func)
        else:
            # weak form of the primal problem
            F = self.weak_residual(problem, W, w, wt, ei_mode=False)

            w, m = self.steady_solve(problem, W, w, F, func=func)

        return W, w, m

    # define functions spaces
    def function_space(self, mesh):
        '''
            Sets up a general mixed function space. We assume there are only two
            variables and the first variable is vector valued. To use something
            different the user can overload this in their Solver.
        '''
        V = VectorFunctionSpace(mesh, 'CG', self.Pu)
        Q = FunctionSpace(mesh, 'CG', self.Pp)
        W = MixedFunctionSpace([V, Q])

        return W

    def weak_residual(self, problem, k, W, w, ww, w_, wt, ei_mode=False):

        print "NO WEAK RESIDUAL PROVIDED: You must define a weak_residual for" \
            + " this code to work."
        sys.exit(1)

    # Refine the mesh based on error indicators
    def adaptive_refine(self, mesh, ei):
        '''
            Take a mesh and the associated error indicators and refine
            adapt_ratio% of cells.
        '''
        gamma = abs(ei.vector().array())

        # Mark cells for refinement
        cell_markers = MeshFunction("bool", mesh, mesh.topology().dim())
        gamma_0 = sorted(gamma,
                         reverse=True)[int(len(gamma) * self.adaptRatio) - 1]
        for c in cells(mesh):
            cell_markers[c] = gamma[c.index()] > gamma_0

        # Refine mesh
        mesh = refine(mesh, cell_markers)

        return mesh

    def adjust_dt(self, t0, T, k):
        '''
            Adjust time step so that we evenly divide the time interval, but
            ensure that the new time step is always smaller than the original.
        '''
        d, r = divmod((T - t0), k)
        if r > DOLFIN_EPS:
            k = (T - t0) / (d + 1)

        return k

    def steady_solve(self, problem, W, w, F, func=False):

        self.start_timing()
        bcs = problem.boundary_conditions(W)

        solve(F == 0, w, bcs)

        if func and adjointer:  # annotation only works with DOLFIN-Adjoint
            m = assemble(problem.functional(W, w), annotate=False)
        elif func:
            m = assemble(problem.functional(W, w))
        else:
            m = None

        return w, m

    def timeStepper(self, problem, t, T, k, W, w, w_, F, func=False):
        '''
            Time stepper for solver using theta-method.
        '''
        # Time loop
        self.start_timing()
        if adjointer:
            adj_start_timestep(t)

        bcs = problem.boundary_conditions(W, t)

        # plot and save initial condition
        self.update(problem, t, W, w_)

        if func and adjointer:  # annotation only works with DOLFIN-Adjoint
            m = k * assemble(problem.functional(W, w_), annotate=False)
        elif func:
            m = k * assemble(problem.functional(W, w_))
        else:
            m = None

        while t < T - k / 2.:
            t += k

            if('update' in dir(problem)):
                bcs = problem.update(W, t)

            if 'pre_step' in dir(self):
                self.pre_step(problem, t, k, W, w)

            solve(F == 0, w, bcs)

            if 'post_step' in dir(self):
                self.post_step(problem, t, k, W, w, w_)

            w_.assign(w)

            # Determine the value of our functional
            if func and adjointer:  # annotation only works with DOLFIN-Adjoint
                m += k * assemble(problem.functional(W, w_), annotate=False)
            elif func:
                m += k * assemble(problem.functional(W, w_))

            if adjointer:  # can only use if DOLFIN-Adjoint has been imported
                adj_inc_timestep(t, finished=t>T-k/2.)

            self.update(problem, t, W, w_)

        return w_, m

    def update(self, problem, t, W, w, dual=False):
        '''
            Saves or plots the data at each time step.
        '''
        # Add to accumulated CPU time
        timestep_cputime = time() - self._time
        self._cputime += timestep_cputime

        if t is not None:  # Store time steps
            self._t.append(t)

        if self.saveSolution:  # Save solution
            self.Save(problem, w, dual=dual)
        elif self.plotSolution:  # Plot solution
            self.Plot(problem, W, w)

        # Check memory usage
        if self.mem:
            print 'Memory usage is:', self.getMyMemoryUsage()

        # Print progress
        if not dual:
            s = 'Time step %d finished in %g seconds, ' \
                % (self._timestep, timestep_cputime)
            s += '%g%% done (t = %g, T = %g).' % (100.0 * (t / problem.T), t,
                                                  problem.T)
            sys.stdout.write('\033[K')
            sys.stdout.write(s + '\r')
            sys.stdout.flush()

            # record current time
            self._time = time()

        # Increase time step
        self._timestep += 1

    def Save(self, problem, w, dual=False):
        '''
            Save a variables associated with a time step. Here we assume there
            are two variables where the first variable is vector-valued and the
            second variable is a scalar. If this doesn't fit the particular
            solvers variables the user will need to overload this function.
        '''
        u = w.split()[0]
        p = w.split()[1]

        if self.saveFrequency != 0 \
                and ((self._timestep - 1) % self.saveFrequency == 0
                     or self._timestep == 0 or self._timestep == problem.T):
            if not dual:
                self._ufile << u
                self._pfile << p
            else:
                self._uDualfile << u
                self._pDualfile << p

    def prefix(self, problem):
        '''
            Obtains the beginning of file naming, e.g. Problem Name, Solver
            Name, dimension, etc.
        '''
        # Return file prefix for output files
        p = problem.__module__.split('.')[-1]
        if problem.mesh.topology().dim() > 2:
            p += '3D'
        else:
            p += '2D'
        s = self.__module__.split('.')[-1]

        return problem.output_location + s + p

    def suffix(self, problem):
        '''
            Obtains the run specific data for file naming, e.g. Nx, k, etc.
        '''

        s = 'T' + str(problem.T)
        if problem.Nx is not None:
            s += 'Nx' + str(problem.Nx)
        if problem.Ny is not None:
            s += 'Ny' + str(problem.Ny)
        if problem.mesh.topology().dim() > 2 and problem.Nz is not None:
            s += 'Nz' + str(problem.Nz)

        s += 'K' + str(problem.k)

        return s

    def file_naming(self, problem, n=-1, opt=False):
        '''
            Names our files for saving variables. Here we assume there are
            two variables where the first variable is vector-valued and the
            second variable is a scalar. If this doesn't fit the particular
            solvers variables the user will need to overload this function.
        '''

        s = 'results/' + self.prefix(problem) + self.suffix(problem)

        if n == -1:
            if opt:
                self._ufile = File(s + '_uOpt.pvd', 'compressed')
                self._pfile = File(s + '_pOpt.pvd', 'compressed')
            else:
                self._ufile = File(s + '_u.pvd', 'compressed')
                self._pfile = File(s + '_p.pvd', 'compressed')
            self._uDualfile = File(s + '_uDual.pvd', 'compressed')
            self._pDualfile = File(s + '_pDual.pvd', 'compressed')
            self.meshfile = File(s + '_mesh.xml')
        else:  # adaptive specific files
            if self.eifile is None:  # error indicators
                self.eifile = File(s + '_ei.pvd', 'compressed')
            self._ufile = File(s + '_u%02d.pvd' % n, 'compressed')
            self._pfile = File(s + '_p%02d.pvd' % n, 'compressed')
            self._uDualfile = File(s + '_uDual%02d.pvd' % n, 'compressed')
            self._pDualfile = File(s + '_pDual%02d.pvd' % n, 'compressed')
            self.meshfile = File(s + '_mesh%02d.xml' % n)

    # this is a separate function so that it can be overloaded
    def Plot(self, problem, W, w):
        '''
            Plots our variables associated with a time step. Here we assume
            there are two variables where the first variable is vector-valued
            and the second variable is a scalar. If this doesn't fit the
            particular solvers variables the user will need to overload this
            function.
        '''
        u = w.split()[0]
        p = w.split()[1]

        if self.vizU is None:
            # Plot velocity and pressure
            self.vizU = plot(u, title='Velocity', rescale=True)
            self.vizP = plot(p, title='Pressure', rescale=True, elevate=0.0)
        else:
            self.vizU.plot(u)
            self.vizP.plot(p)

    def getMyMemoryUsage(self):
        '''
            Determines how much memory we are using.
        '''
        mypid = getpid()
        mymemory = getoutput('ps -o rss %s' % mypid).split()[1]
        return mymemory

    def start_timing(self):
        '''
            Start timing, will be paused automatically during update
            and stopped when the end-time is reached.
        '''
        self._time = time()
