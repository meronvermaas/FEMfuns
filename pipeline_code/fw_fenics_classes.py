from dolfin import *
import numpy as np
from time import time
import dispersion_spectra as disp
from collections import Counter
from copy import deepcopy
import os

#parameters["ghost_mode"] = "shared_facet"

def init_funcspace(mesh):

    mesh_funcspace_const = FunctionSpace(mesh, 'DG', 0)
    mesh_funcspace_scalar = FunctionSpace(mesh, 'CG', 2)
    V_ele = FiniteElement("CG", mesh.ufl_cell(), 2)
    mesh_funcspace_mixed = FunctionSpace(mesh, MixedElement([V_ele, V_ele]))

    funcspaces = {
            'const' : mesh_funcspace_const,
            'scalar' : mesh_funcspace_scalar,
            'mixed' : mesh_funcspace_mixed,
    }

    return funcspaces


class FEM_simulation(object):

    def __init__(self, params, mesh_filename = None, mesh_materials_filename = None, mesh_boundaries_filename = None, h5file = None):
        self.params = params
        self.mesh_filename = mesh_filename
        self.mesh_materials_filename = mesh_materials_filename
        self.mesh_boundaries_filename = mesh_boundaries_filename
        self.h5file = h5file
        self.geometry = FEM_geometry(self)
        self.geometry.init_admittivity(self)

        #self.geometry.init_admittivity(self)
        self.e0 = 8.854187817e-12 #F/m
        self.solutions = []

        self.set_empty_params()

        print "done setting up"

    def set_empty_params(self):
        if not hasattr(self.params, 'boundary_markers'):
            self.params.boundary_markers = {}
        if not hasattr(self.params, 'monopole_list'):
            self.params.monopole_list = [None]

    def main(self, solver_method='cg', preconditioner='amg'):
        tic = time()
        print "start solving, this will take some time..."

        self.solver_method = solver_method
        self.preconditioner = preconditioner
        # select all options w.r.t. source, material and electrode
        for monopoles_dict in [self.params.monopole_list[0]]:

            if self.params.material == "resistive":
                self.get_weak_form()
                self.apply_source(monopoles_dict)

                self.phi = Function(self.geometry.funcspace['scalar'])
                iterations = self.get_solution()
                self.solutions.append([self.phi.vector().get_local(), \
                                       [self.params.material,
                                        deepcopy(self.params.volume_markers),
                                        deepcopy(self.params.boundary_markers),
                                        deepcopy(monopoles_dict)
                                       ]
                                      ])

            else:
                if hasattr(self.params, 'frequencies'):
                    for self.frequency in self.params.frequencies:
                        if self.params.material == "dispersive" or not hasattr(self.geometry, 'epsilon'):
                            self.geometry.init_admittivity(self)
                        omega = 2.*np.pi*self.frequency
                        self.get_complex_weak_form(omega)
                        self.apply_complex_source(monopoles_dict)

                        self.phi = Function(self.geometry.funcspace['mixed'])
                        iterations = self.get_solution()
                        self.solutions.append([self.phi.vector().get_local(), \
                                       [self.params.material,
                                        deepcopy(self.params.volume_markers),
                                        deepcopy(self.params.boundary_markers),
                                        deepcopy(monopoles_dict)
                                       ]
                                      ])
                else:
                    print "add [frequencies] attribute to parameters"

            print (time() - tic)/60, "minutes to solve linear system in ", iterations, "iterations"

    def get_solution(self):

        # standard solver. if other needed (no convergence) choose from list_krylov_solver_methods() and list_krylov_solver_preconditioners()

        #if mixed element is used, CG won't converge
        if self.params.material != "resistive" and self.solver_method == 'cg' and self.preconditioner == 'amg':
            solver = KrylovSolver('bicgstab', 'ilu')

        else:
            solver = KrylovSolver(self.solver_method, self.preconditioner)

        prm = solver.parameters
        prm.absolute_tolerance = 1E-8
        #prm.relative_tolerance = 1E-6
        prm.maximum_iterations = 10000
        prm.monitor_convergence = True
        #info(solver.parameters, True)
        #set_log_level(PROGRESS)

        # if mumps solver is wanted:
        #iterations = solve(A, self.phi.vector(), b, 'mumps')
        try:
            iterations = solver.solve(self.A, self.phi.vector(), self.b)
        except RuntimeError:
            print "RuntimeError: Solver did not converge, trying different solver now"
            solver = KrylovSolver('gmres', 'ilu')
            prm = solver.parameters
            prm.absolute_tolerance = 1E-8
            #prm.relative_tolerance = 1E-6
            prm.maximum_iterations = 10000
            prm.monitor_convergence = True
            #info(solver.parameters, True)
            #set_log_level(PROGRESS)
            iterations = solver.solve(self.A, self.phi.vector(), self.b)

        return iterations

    def get_weak_form(self):
        u = TrialFunction(self.geometry.funcspace['scalar'])
        v = TestFunction(self.geometry.funcspace['scalar'])

        ## Define weak form
        self.a = inner(self.geometry.sigma * grad(u), grad(v))*self.geometry.dx
        self.L = Constant(0.)*v*self.geometry.dx

        self.apply_bcs(u, v)

        self.A = assemble(self.a)
        self.b = assemble(self.L)

    def get_complex_weak_form(self, omega):
        (u_r, u_i) = TrialFunctions(self.geometry.funcspace['mixed'])
        (v_r, v_i) = TestFunctions(self.geometry.funcspace['mixed'])

        ## Define weak form
        a_r = inner(self.geometry.sigma*nabla_grad(u_r), nabla_grad(v_r))*self.geometry.dx - \
              inner(self.geometry.sigma*nabla_grad(u_i), nabla_grad(v_i))*self.geometry.dx - \
              inner(self.e0*omega*self.geometry.epsilon*nabla_grad(u_r), nabla_grad(v_i))*self.geometry.dx -\
              inner(self.e0*omega*self.geometry.epsilon*nabla_grad(u_i), nabla_grad(v_r))*self.geometry.dx

        a_i = inner(self.geometry.sigma*nabla_grad(u_r), nabla_grad(v_i))*self.geometry.dx + \
              inner(self.geometry.sigma*nabla_grad(u_i), nabla_grad(v_r))*self.geometry.dx + \
              inner(self.e0*omega*self.geometry.epsilon*nabla_grad(u_r), nabla_grad(v_r))*self.geometry.dx -\
              inner(self.e0*omega*self.geometry.epsilon*nabla_grad(u_i), nabla_grad(v_i))*self.geometry.dx
        self.a = a_r + a_i
        self.L = Constant(0.)*v_r*self.geometry.dx + Constant(0.)*v_i*self.geometry.dx

        self.apply_complex_bcs(u_r, u_i, v_r, v_i)

        self.A = assemble(self.a)
        self.b = assemble(self.L)

    def apply_bcs(self, u, v):
        for bound in self.params.boundary_markers.items():
            if len(bound[1]) == 4:
                if bound[1][3] == 'ext':
                    self.a += Constant(bound[1][2]) * u * v * self.geometry.ds(int(bound[1][0]))
                    self.L += Constant(bound[1][2]) * Constant(bound[1][1]) * v * self.geometry.ds(int(bound[1][0]))
                elif bound[1][3] == 'int':
                    self.a += Constant(bound[1][2]) * u * v * self.geometry.dS(int(bound[1][0]))
                    self.L += Constant(bound[1][2]) * Constant(bound[1][1]) * v * self.geometry.dS(int(bound[1][0]))

    def apply_complex_bcs(self, u_r, u_i, v_r, v_i):
        for bound in self.params.boundary_markers.items():
            if len(bound[1]) == 4:
                if bound[1][3] == 'ext':
                    a_r = Constant(np.real(bound[1][2])) * u_r * v_r * self.geometry.ds(int(bound[1][0])) \
                        - Constant(np.real(bound[1][2])) * u_i * v_i * self.geometry.ds(int(bound[1][0])) \
                        - Constant(np.imag(bound[1][2])) * u_i * v_r * self.geometry.ds(int(bound[1][0])) \
                        - Constant(np.imag(bound[1][2])) * u_r * v_i * self.geometry.ds(int(bound[1][0]))

                    a_i = Constant(np.real(bound[1][2])) * u_i * v_r * self.geometry.ds(int(bound[1][0])) \
                        + Constant(np.real(bound[1][2])) * u_r * v_i * self.geometry.ds(int(bound[1][0])) \
                        + Constant(np.imag(bound[1][2])) * u_i * v_i * self.geometry.ds(int(bound[1][0])) \
                        - Constant(np.imag(bound[1][2])) * u_i * v_i * self.geometry.ds(int(bound[1][0]))

                    self.a += a_r+a_i

                    L_r = - Constant(np.real(bound[1][2])) * Constant(bound[1][1]) * v_r * self.geometry.ds(int(bound[1][0])) \
                          + Constant(np.imag(bound[1][2])) * Constant(bound[1][1]) * v_i * self.geometry.ds(int(bound[1][0]))

                    L_i = - Constant(np.real(bound[1][2])) * Constant(bound[1][1]) * v_i * self.geometry.ds(int(bound[1][0])) \
                          - Constant(np.imag(bound[1][2])) * Constant(bound[1][1]) * v_r * self.geometry.ds(int(bound[1][0]))

                    self.L += L_r+L_i
                elif bound[1][3] == 'int':
                    a_r = Constant(np.real(bound[1][2])) * u_r * v_r * self.geometry.dS(int(bound[1][0])) \
                        - Constant(np.real(bound[1][2])) * u_i * v_i * self.geometry.dS(int(bound[1][0])) \
                        - Constant(np.imag(bound[1][2])) * u_i * v_r * self.geometry.dS(int(bound[1][0])) \
                        - Constant(np.imag(bound[1][2])) * u_r * v_i * self.geometry.dS(int(bound[1][0]))

                    a_i = Constant(np.real(bound[1][2])) * u_i * v_r * self.geometry.dS(int(bound[1][0])) \
                        + Constant(np.real(bound[1][2])) * u_r * v_i * self.geometry.dS(int(bound[1][0])) \
                        + Constant(np.imag(bound[1][2])) * u_i * v_i * self.geometry.dS(int(bound[1][0])) \
                        - Constant(np.imag(bound[1][2])) * u_i * v_i * self.geometry.dS(int(bound[1][0]))

                    self.a += a_r+a_i

                    L_r = - Constant(np.real(bound[1][2])) * Constant(bound[1][1]) * v_r * self.geometry.dS(int(bound[1][0])) \
                          + Constant(np.imag(bound[1][2])) * Constant(bound[1][1]) * v_i * self.geometry.dS(int(bound[1][0]))

                    L_i = - Constant(np.real(bound[1][2])) * Constant(bound[1][1]) * v_i * self.geometry.dS(int(bound[1][0])) \
                          - Constant(np.imag(bound[1][2])) * Constant(bound[1][1]) * v_r * self.geometry.dS(int(bound[1][0]))

                    self.L += L_r+L_i


    def apply_source(self, monopoles_dict):
        if monopoles_dict:
            for monopole in monopoles_dict['monopoles']:
                tmp_point = Point(*monopole[0:3])

                delta = PointSource(self.geometry.funcspace['scalar'], tmp_point, monopole[3])
                delta.apply(self.b)

    def apply_complex_source(self, monopoles_dict):
        if monopoles_dict:
            for monopole in monopoles_dict['monopoles']:
                tmp_point = Point(*monopole[0:3])

                delta0 = PointSource(self.geometry.funcspace['mixed'].sub(0), tmp_point, monopole[3])
                delta0.apply(self.b)
                delta1 = PointSource(self.geometry.funcspace['mixed'].sub(1), tmp_point, monopole[3])
                delta1.apply(self.b)

    def get_poisson_values(self,coords):

        values = []
        self.phi.set_allow_extrapolation(True)

        if len(np.shape(coords)) == 2:
            for idx in range(coords.shape[0]):
                phi_tmp = self.phi(coords[idx,:])
                if phi_tmp.size == 1: #self.phi.function_space().ufl_element().family() == 'Mixed':
                    values.append(phi_tmp)
                else:
                    values.append(phi_tmp[0]+phi_tmp[1]*1j)
        else:
            f = project(Constant(1.), self.geometry.funcspace['scalar'])
            for idx in coords:
                phi_tmp = assemble(self.phi*self.geometry.ds(idx))/assemble(f*self.geometry.ds(idx))
                values.append(phi_tmp)

        return np.asarray(values)

    def get_error_values(self,ana_values,coords):

        fem_values = self.get_poisson_values(coords)
        diff = ana_values - fem_values
        self.RMSE = np.sqrt(np.mean((diff)**2))
        self.RD = np.mean(np.abs(diff)/np.amax(np.abs(ana_values)))

        print "RMSE ", self.RMSE, "RD", self.RD

    def save_pvds(self, allsols = False, filename = None):
        '''
        save either the current (last solved) solution (default)
        or set allsolls = True and save all solutions of simulations
        automatic filename (based on params) which is long (default)
        or manual filename=pvds_dir/string.pvd input
        '''

        if not os.path.exists('pvds_dir'):
            os.makedirs('pvds_dir')

        if allsols and filename:
            print "saving all solutions with one filename does not work. filename is useable for the latest simulation only" 
            return

        if allsols:
            for solution in self.solutions:
                if not filename:
                    filename = 'pvds_dir/' + self.mesh_filename.split('/',1)[1].split('.',1)[0]
                    filename += '_'+solution[1][0] #add material type
                    for idx in range(1,4):
                        if solution[1][idx]:
                            filename += '_'+''.join("%s=%r" % (key,val) for (key,val) in solution[1][idx].iteritems())
                    filename = filename.replace(" ", "")
                    filename = filename.replace("'", "")
                    filename += '.pvd'

                if solution[1][0] == 'resistive':
                    projfunc = Function(self.geometry.funcspace['scalar'])
                else:
                    projfunc = Function(self.geometry.funcspace['mixed'])
                projfunc.vector().set_local(solution[0])
                File(filename) << projfunc
        else:
            if not filename:
                filename = 'pvds_dir/' + self.mesh_filename.split('/',1)[1].split('.',1)[0]
                filename += '_' + self.params.material
                filename += '_' + ''.join("%s=%r" % (key,val) for (key,val) in self.params.volume_markers.iteritems())
                if self.params.boundary_markers:
                    filename += '_' + ''.join("%s=%r" % (key,val) for (key,val) in self.params.boundary_markers.iteritems())
                if self.params.monopole_list[-1]:
                    filename += '_' + ''.join("%s=%r" % (key,val) for (key,val) in self.params.monopole_list.iteritems())
                filename = filename.replace(" ", "")
                filename = filename.replace("'", "")
                filename += '.pvd'
            File(filename) << self.phi

    def find_domain_midpoint(self,domain, marker):
        It_domains = SubsetIterator(domain, marker)
        coords = []
        normals = []
        for el in It_domains:
            coords.append([el.midpoint().x(),el.midpoint().y(),el.midpoint().z()])
            el_facet = Facet(el.mesh(), el.index())
            normals.append(el_facet.normal())

        normal_coords = []
        for ii in range(len(normals)):
            normal_coords.append([normals[ii][0],normals[ii][1],normals[ii][2]])
        normal_coords = np.asarray(normal_coords)
        avgnormal = [np.mean(normal_coords[:,0]),np.mean(normal_coords[:,1]),np.mean(normal_coords[:,2])]

        coords = np.asarray(coords)
        midx = (coords[:,0].max()+coords[:,0].min())/2.
        midy = (coords[:,1].max()+coords[:,1].min())/2.
        midz = (coords[:,2].max()+coords[:,2].min())/2.
        domain_midpoint = np.array([midx,midy,midz])
        return domain_midpoint, coords, avgnormal

    def source_locations(self, domainmarker, locationdomain, locationmarker, distance, inner_distance, nr_points, orientations):
        '''
        Compute points where to apply a pointsource
        domainmarker: domain where points are to be positioned
        locationdomain: cell or facet domain
        locationmarker: region where distance to sources is taken (e.g. electrode)
        distance: distance of selected points to the locationmarker
        inner_distance: how close can individual points be to each other
        nr_points: amount of point locations to produce
        orientations: list containing (a combination of) 'rad', 'tan', or 'mono'

        note that refinement might be needed because of small distance between sink and source

        '''

        # midpoint of locationmarker
        midpoint, coords, self.avgnormal = self.find_domain_midpoint(locationdomain, locationmarker)
        p = Point(midpoint)

        # loop over cells in domain where points are to be selected
        sourcepoints = []
        It_domains = SubsetIterator(self.geometry.domain,domainmarker)
        for c in It_domains:
            tmpcoords = Point(c.midpoint().x(), c.midpoint().y(), c.midpoint().z())

            # break loop if requested nr of points has been added
            if len(sourcepoints) == nr_points:
                break

            # distance from current cell to locationmarker
            threshold = distance/10.
            if p.distance(tmpcoords) > distance or p.distance(tmpcoords) < distance-threshold:
                continue

            # loop over already created points and check distance between each other
            skip = False
            for pnt in sourcepoints:
                if pnt.distance(tmpcoords) < inner_distance:
                    skip = True
                    break
            if skip:
                continue

            sourcepoints.append(tmpcoords)
        if nr_points > len(sourcepoints):
            print "Less locations than requested created, adjust distance or inner_distance if more are needed."

        self.sourcepoints = np.asarray(sourcepoints)
        self.add_monopole_list(orientations)

    def add_monopole_list(self, orientations):
        for pnt in self.sourcepoints:
            for orientation in orientations:
                if orientation == 'rad':
                    source_coord = [(pnt[0]+0.05*self.avgnormal[0]), (pnt[1]+0.05*self.avgnormal[1]), (pnt[2]+0.05*self.avgnormal[2]), 1.]
                    snk_coord = [(pnt[0]-0.05*self.avgnormal[0]), (pnt[1]-0.05*self.avgnormal[1]), (pnt[2]-0.05*self.avgnormal[2]), -1.]
                    self.params.monopole_list.append({'monopoles': [source_coord, snk_coord], 'name':'rad'})
                elif orientation == 'tan':
                    source_coord = [(pnt[0]+0.05*-self.avgnormal[2]), (pnt[1]+0.05*self.avgnormal[1]), (pnt[2]+0.05*self.avgnormal[0]), 1.]
                    snk_coord = [(pnt[0]-0.05*-self.avgnormal[2]), (pnt[1]-0.05*self.avgnormal[1]), (pnt[2]-0.05*self.avgnormal[0]), -1.]
                    self.params.monopole_list.append({'monopoles': [source_coord, snk_coord], 'name':'tan'})
                elif orientation == 'mono':
                    source_coord = [pnt[0],pnt[1],pnt[2], 1.]
                    self.params.monopole_list.append({'monopoles': [source_coord], 'name':'mono'})
                else:
                    print "name types of orientations list of sources can be 'mono', 'rad' or 'tan'"

    def set_dispersive_params(self):
        '''
        Change parameters to run with dispersive material.
        Conductivity value should now be a string of the name in the cole cole table
        Type help(disp.get_gabriel_params) to look for the names
        '''

        self.params.material = 'dispersive'
        if not hasattr(self.params, 'frequencies'):
            self.params.frequencies = [1000] # 1000 Hz as default value

        disp_tissues = [tissue[0] for i, tissue in enumerate(disp.colecole_params)]
        for idx, vol in enumerate(self.params.volume_markers.items()):
            found = False
            for tissue in disp_tissues:
                # check if name is part of colecole name, or some standard exceptions
                if vol[0].lower() in tissue.lower() or \
                   (vol[0].lower() == 'csf' and tissue == 'Cerebro_Spinal_Fluid') or \
                   (vol[0].lower() == 'ventricles' and tissue == 'Cerebro_Spinal_Fluid') or \
                   (vol[0].lower() == 'skull' and tissue == 'Bone_Cortical') or \
                   (vol[0].lower() == 'scalp' and tissue == 'Skin_Dry'):
                    self.params.volume_markers.items()[idx][1][1] = tissue
                    found = True
                    break
            if not found:
                print 'could not automatically set ', vol[0], ' to dispersive setting \n'
                print 'please do this manually, options can be found here: help(disp.get_gabriel_params)'

class FEM_geometry(FEM_simulation):

    def __init__(self, fem_sim): #, params, mesh_filename = None, mesh_materials_filename = None, mesh_boundaries_filename = None, h5file = None): #, mesh_filename=None, mesh_materials_filename=None, mesh_boundaries_filename = None, h5file=None):

        if fem_sim.h5file:
            self.mesh = Mesh()
            hdf = HDF5File(self.mesh.mpi_comm(), fem_sim.h5file, "r")
            hdf.read(self.mesh, "/mesh", False)
            self.domain = MeshFunction("size_t", self.mesh)
            hdf.read(self.domain, "/domains")
            self.boundaries = MeshFunction("size_t", self.mesh)
            hdf.read(self.boundaries, "/boundaries")
        self.init_geometry(fem_sim)

    def init_geometry(self,fem_sim):
        """initialize geometry and materials"""

        if not hasattr(self,'mesh'):
            self.mesh = Mesh(fem_sim.mesh_filename)

        if not hasattr(self,'domain'):
            self.domain = MeshFunction('size_t', self.mesh, fem_sim.mesh_materials_filename)

        if not hasattr(self,'boundaries'):
            if hasattr(fem_sim, 'mesh_boundaries_filename'):
                self.boundaries = MeshFunction('size_t', self.mesh, fem_sim.mesh_boundaries_filename)
            else:
                self.boundaries = compute_boundaries(self.mesh, self.domain)

        self.dx = dx(domain = self.mesh, subdomain_data = self.domain)
        self.dS = dS(domain = self.mesh, subdomain_data = self.boundaries)
        self.ds = ds(domain = self.mesh, subdomain_data = self.boundaries)

        self.funcspace = init_funcspace(self.mesh)


    def init_admittivity(self,fem_sim):
        # make sigma and if needed epsilon function
        self.sigma = Function(self.funcspace['const'])
        self.conductivities = np.empty(len(fem_sim.params.volume_markers),dtype=float)

        if ((fem_sim.params.material == "capacitive") or (fem_sim.params.material == "dispersive")) and not hasattr(self, 'epsilon'):
            self.epsilon = Function(fem_sim.geometry.funcspace['const'])
            self.permittivities = np.empty(len(fem_sim.params.volume_markers),dtype=float)

        for idx, dommark in enumerate(sorted(fem_sim.params.volume_markers.items(), key=lambda x: x[1])):
            if fem_sim.params.material == "resistive":
                self.conductivities[idx] = dommark[1][1] * fem_sim.params.unit
            elif fem_sim.params.material == "capacitive":
                self.conductivities[idx] = np.real(dommark[1][1]) * fem_sim.params.unit
                self.permittivities[idx] = np.imag(dommark[1][1]) * fem_sim.params.unit
            elif fem_sim.params.material == "dispersive": # and hasattr(femsim, 'frequency'):
                try:
                    tmpeps, tmpsig = disp.colecole_gabriel(fem_sim.frequency, dommark[1][1])
                    self.conductivities[idx] = tmpsig * fem_sim.params.unit
                    self.permittivities[idx] = tmpeps * fem_sim.params.unit
                except AttributeError:
                    return
            else:
                raise ValueError("no known tissue type, options are resistive, capacitive and dispersive")

        palette = np.sort(np.unique(self.domain.array()[:]))
        index = np.digitize(self.domain.array().ravel(),palette,right=True)
        self.sigma.vector()[:] = self.conductivities[index]
        if fem_sim.params.material == "capacitive" or fem_sim.params.material == "dispersive":
            self.epsilon.vector()[:] = self.permittivities[index]

    def cells_in_region(self, distance, point):

        refine_markers = MeshFunction('bool',self.mesh,3) #CellFunction('bool', mesh)
        refine_markers.set_all(False)

        for c in cells(mesh):
            if c.midpoint().x() < xmax and c.midpoint().x() > xmin and \
             c.midpoint().y() < ymax and c.midpoint().y() > ymin and \
             c.midpoint().z() < zmax and c.midpoint().z() > zmin:
                refine_markers[c] = True


def make_pulse(time,width,npoints):
    Fs = npoints/time #sampling frequency
    t = np.arange(-time/2.,time/2.,1./Fs)
    x = np.zeros(len(t))
    x[int(len(t)/2-width*Fs/2.)+1:int(len(t)/2+width*Fs/2.)+1] = 1
    x -= np.mean(x)
    X = np.fft.fft(x)
    freqs = np.fft.fftfreq(len(t), d=1./Fs)
    return t,x,freqs,X

def make_cos(time,npoints):
    Fs = npoints/time #sampling frequency
    t = np.arange(-time/2.,time/2.,1./Fs)
    x = np.cos(2.*np.pi*(np.arange(-npoints/2.,npoints/2.)/npoints))
    X = np.fft.fft(x)
    freqs = np.fft.fftfreq(len(t), d=1./Fs)
    return t,x,freqs,X

def make_alphafunction(time,npoints):
    tau = 1e-3 #1ms
    t = np.linspace(0,time,npoints)
    alphafun = t/tau * np.exp(-t/tau)
    alphafun -= np.mean(alphafun)
    Fs = npoints/time
    X = np.fft.fft(alphafun)
    freqs = np.fft.fftfreq(npoints, d=1./Fs)
    return t,alphafun, freqs, X

#timeduration = 1e-2
#npoints = 2000
#width = 1e-3
#t,x,frequencies,X = make_pulse(timeduration,width,npoints)
#t,x,frequencies,X = make_cos(1e-6,150)
#t,x,frequencies,X = make_alphafunction(1e-2,100)

class CubeDomain(SubDomain):
        def __init__(self, ddict, size):
                SubDomain.__init__(self)
                self.xmin, self.xmax = ddict['xmin']-size/2., ddict['xmax']+size/2.
                self.ymin, self.ymax = ddict['ymin']-size/2., ddict['ymax']+size/2.
                self.zmin, self.zmax = ddict['zmin']-size/2., ddict['zmax']+size/2.
        def inside(self, x, on_boundary):
                return ( between(x[0], (self.xmin, self.xmax)) and
                                 between(x[1], (self.ymin, self.ymax)) and
                                 between(x[2], (self.zmin, self.zmax)) )

class FEM_subgeometry:

        def __init__(self, geometry, size, marker):
                self.parent_geometry = geometry
                self.size = size
                self.marker = marker
                self.markers = geometry.markers
                self.submesh_geometry()

        def submesh_geometry(self):

                new_region = CellFunction('size_t', self.parent_geometry.mesh)
                new_region.set_all(0)

                minmax_dict = find_domain_midpoint(self.parent_geometry.domain, self.parent_geometry.markers[self.marker])

                cuberegion = CubeDomain(minmax_dict, self.size)
                cuberegion.mark(new_region,1)
                self.mesh = SubMesh(self.parent_geometry.mesh, new_region, 1)

                ######
                ## transfer domain markers
                self.domain = CellFunction('size_t', self.mesh)
                self.domain.set_all(0)

                cmap = self.mesh.data().array('parent_cell_indices', 3)

                for sub_mesh_cell in cells(self.mesh):
                        parent_cell = cmap[sub_mesh_cell.index()]
                        self.domain.array()[sub_mesh_cell.index()] = self.parent_geometry.domain.array()[parent_cell]

                ######
                ## recompute facet markers
                self.facet_markers = compute_boundaries(self.mesh, self.domain)

                self.dx = dx(subdomain_data = self.domain)
                self.dS = dS(subdomain_data = self.facet_markers)
                self.ds = ds(domain = self.mesh)

                self.funcspace = init_funcspace(self.mesh)

def compute_boundaries(mesh, domains):
    ''' Calculate the interior face boundaries '''

    print "Computing boundaries..."

    tic = time()
    domains_counter = Counter(domains.array())
    boundaries = FacetFunction('size_t', mesh)

    D = mesh.topology().dim()
    mesh.init(D-1,D)

    for idx_domain_value in range(len(domains_counter)-1,-1,-1):
        # get value of the subregion
        domain_value = domains_counter.most_common()[idx_domain_value][0]
        # loop over cells
        for c in cells(mesh):
            # check if cell is in current subdomain_value region
            if domains[c] ==  domain_value :
                for f in facets(c):
                    # check neighbouring cells
                    for d in f.entities(D):
                        # check if facet is at inner boundary and has not been assigned yet
                        if domains[int(d)] != domains[c] and boundaries[f] == 0:
                            ## Needed to have separate boundary between electrode recording surface and outer surface
                            ## custom step now, typically this def won't be used anyway:
                            if (domains[c] == 4 or domains[c] == 5 or domains[c] == 6) and domains[int(d)] == 2:
                                continue
                            boundaries[f] = domain_value
    print "Computing boundaries took", (time()-tic)/60., "minutes"
    return boundaries
