from dolfin import *
import numpy as np
from time import time
import dispersion_spectra as disp
import cantrell as ct
from collections import Counter
from copy import deepcopy
import os

#parameters["ghost_mode"] = "shared_facet"

def init_funcspace(mesh, order):

    mesh_funcspace_const = FunctionSpace(mesh, 'DG', 0)
    mesh_funcspace_scalar = FunctionSpace(mesh, 'CG', order)

    funcspaces = {
            'const' : mesh_funcspace_const,
            'scalar' : mesh_funcspace_scalar
    }

    return funcspaces


class FEM_simulation(object):

    def __init__(self, params, mesh_filename = None, mesh_materials_filename = None, mesh_boundaries_filename = None, h5file = None, order = 1):
        self.params = params
        self.order = order
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

        print( "done setting up")

    def set_empty_params(self):
        if not hasattr(self.params, 'boundary_markers'):
            self.params.boundary_markers = {}
        if not hasattr(self.params, 'monopole_list'):
            self.params.monopole_list = [None]

    def main(self, lagrange_multiplier = False, solver_method='cg', preconditioner='ilu', relative_tolerance=None, absolute_tolerance=1e-8, maximum_iterations=1000):
        print( "start solving, this will take some time...")

        self.solver_method = solver_method
        self.preconditioner = preconditioner
        self.absolute_tolerance = absolute_tolerance
        self.maximum_iterations = maximum_iterations
        self.relative_tolerance = relative_tolerance

        # select all options w.r.t. source, material and electrode
        for monopoles_dict in self.params.monopole_list:
            tic = time()
            if self.params.material == "resistive":
                if not lagrange_multiplier:
                    self.get_weak_form()
                    self.apply_source(monopoles_dict)

                    self.phi = Function(self.geometry.funcspace['scalar'])
                else:
                    if 'lagrange_mp' not in self.geometry.funcspace.keys():
                        P1 = FiniteElement("CG", self.geometry.mesh.ufl_cell(), self.order)
                        R = FiniteElement('Real', self.geometry.mesh.ufl_cell(), 0)
                        mixedelement_list = [P1]
                        for ii in [len(x) for x in self.params.boundary_markers.values()]:
                            if ii == 3:
                                mixedelement_list.append(R)
                        mesh_funcspace_lagrange_mp = FunctionSpace(self.geometry.mesh, MixedElement(mixedelement_list))
                        self.geometry.funcspace['lagrange_mp'] = mesh_funcspace_lagrange_mp

                    self.get_lagrange_mp_weak_form()
                    self.apply_lagrange_mp_source(monopoles_dict)

                    self.phi = Function(self.geometry.funcspace['lagrange_mp'])

                iterations = self.get_solution()
                if lagrange_multiplier:
                    phitmp = self.phi.split(deepcopy=True)
                    self.phi = phitmp[0]
                self.solutions.append([self.phi.vector().get_local(), \
                                       [self.params.material,
                                        deepcopy(self.params.volume_markers),
                                        deepcopy(self.params.boundary_markers),
                                        deepcopy(monopoles_dict)
                                       ]
                                      ])

            else:
                if hasattr(self.params, 'frequencies'):
                    if 'mixed' not in self.geometry.funcspace.keys():
                        P1 = FiniteElement("CG", self.geometry.mesh.ufl_cell(), self.order)
                        mesh_funcspace_mixed = FunctionSpace(self.geometry.mesh, MixedElement([P1,P1]))
                        self.geometry.funcspace['mixed'] = mesh_funcspace_mixed

                    all_phi = []
                    for self.frequency in np.unique(np.abs(self.params.frequencies)):

                        if self.params.material == "dispersive" or not hasattr(self.geometry, 'epsilon'):
                            self.geometry.init_admittivity(self)

                        omega = 2.*np.pi*self.frequency
                        self.get_complex_weak_form(omega)

                        self.apply_complex_source(monopoles_dict)

                        self.phi = Function(self.geometry.funcspace['mixed'])
                        iterations = self.get_solution()

                        if hasattr(self.params, 'x'):
                            all_phi.append(self.phi.vector().get_local())
                            print((len(all_phi)/len(np.unique(np.abs(self.params.frequencies)))*100), ' percent of the frequencies is calculated')
                        else:
                            self.solutions.append([self.phi.vector().get_local(), \
                                                   [self.params.material,
                                                    deepcopy(self.params.volume_markers),
                                                    deepcopy(self.params.boundary_markers),
                                                    deepcopy(monopoles_dict)
                                                   ]
                                                  ])

                    if hasattr(self.params, 'x'):

                        all_phi = np.asarray(all_phi)
                        all_phi.resize([all_phi.shape[0],int(all_phi.shape[1]/2),2])
                        all_phi = np.vectorize(complex)(all_phi[:,:,0], all_phi[:,:,1])
                        #negative frequencies
                        all_phi_conj = np.flip(np.conj(all_phi[1:,:]), 0)
                        if self.params.x.shape[0]%2 != 1:
                            all_phi = all_phi[0:-1,:]
                        all_phi = np.concatenate((all_phi,all_phi_conj), axis=0)*self.params.X[:,None]
                        ifft_all_phi = np.fft.ifftn(all_phi,axes=(0,))

                        self.solutions.append([ifft_all_phi, \
                                               [self.params.material,
                                                deepcopy(self.params.volume_markers),
                                                deepcopy(self.params.boundary_markers),
                                                deepcopy(monopoles_dict),
                                                self.params.t,
                                                self.params.x
                                               ]
                                              ])

                else:
                    print( "add [frequencies] attribute to parameters")

            print( (time() - tic)/60, "minutes to solve linear system in ", iterations, "iterations")

        ## free up some memory
        #self.geometry.funcspace.pop('mixed', None)
        #self.geometry.funcspace.pop('lagrange_mp', None)

    def get_solution(self):

        # standard solver. if other needed (no convergence) choose from list_krylov_solver_methods() and list_krylov_solver_preconditioners()

        #if mixed element is used, CG won't converge
        solver = KrylovSolver(self.solver_method, self.preconditioner)

        prm = solver.parameters
        prm['absolute_tolerance'] = self.absolute_tolerance
        if self.relative_tolerance:
            prm['relative_tolerance'] = self.relative_tolerance
        prm['maximum_iterations'] = self.maximum_iterations
        prm['monitor_convergence'] = True
        #info(solver.parameters, True)
        #set_log_level(PROGRESS)

        # if mumps solver is wanted:
        #iterations = solve(A, self.phi.vector(), b, 'mumps')
        try:
            iterations = solver.solve(self.A, self.phi.vector(), self.b)
        except RuntimeError:
            print( "RuntimeError: Solver did not converge, trying with gmres and ilu and relative tolerance at 1e-1 (this will result in a lower accuracy)")
            solver = KrylovSolver('gmres', 'ilu')
            prm = solver.parameters
            prm['absolute_tolerance'] = 1E-8
            prm['relative_tolerance'] = 1E-1
            prm['maximum_iterations'] = 10000
            prm['monitor_convergence'] = True

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

        self.A.ident_zeros()
        # apply the Dirichlet BC on assmbled matrix
        [bc.apply(self.A, self.b) for bc in self.bcs]

    def get_lagrange_mp_weak_form(self):
        u = TrialFunctions(self.geometry.funcspace['lagrange_mp'])
        v = TestFunctions(self.geometry.funcspace['lagrange_mp'])

        ## Define weak form
        self.a = inner(self.geometry.sigma * grad(u[0]), grad(v[0]))*self.geometry.dx
        self.L = Constant(0.)*v[0]*self.geometry.dx

        self.apply_lagrange_mp_bcs(u, v)

        self.A = assemble(self.a, keep_diagonal=True)
        self.b = assemble(self.L, keep_diagonal=True)

        self.A.ident_zeros()

        # apply the Dirichlet BC on assmbled matrix
        [bc.apply(self.A, self.b) for bc in self.bcs]


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

        self.apply_complex_bcs(u_r, u_i, v_r, v_i, omega)

        self.A = assemble(self.a)
        self.b = assemble(self.L)

    def apply_bcs(self, u, v):

        self.bcs = []
        for bound in self.params.boundary_markers.items():
            if len(bound[1]) == 4:
                if bound[1][3] == 'ext':
                    self.a += Constant(bound[1][2]) * u * v * self.geometry.ds(int(bound[1][0]))
                    self.L += Constant(bound[1][2]) * Constant(bound[1][1]) * v * self.geometry.ds(int(bound[1][0]))
                elif bound[1][3] == 'int':
                    self.a += Constant(bound[1][2]) * avg(u) * avg(v) * self.geometry.dS(int(bound[1][0]))
                    self.L += Constant(bound[1][2]) * Constant(bound[1][1]) * avg(v) * self.geometry.dS(int(bound[1][0]))
            elif len(bound[1]) == 2:
                self.bcs.append(DirichletBC(self.geometry.funcspace['scalar'], Constant(bound[1][1]), self.geometry.boundaries, bound[1][0]))
            else:
                #no BCs were applied, do nothing
                continue

    def apply_lagrange_mp_bcs(self, u, v):

        self.bcs = []
        lagrange_cnt = 0
        for bound in self.params.boundary_markers.items():
            if len(bound[1]) == 4:
                if bound[1][3] == 'ext':
                    self.a += Constant(bound[1][2]) * u[0] * v[0] * self.geometry.ds(int(bound[1][0]))
                    self.L += Constant(bound[1][2]) * Constant(bound[1][1]) * v[0] * self.geometry.ds(int(bound[1][0]))
                elif bound[1][3] == 'int':
                    self.a += Constant(bound[1][2]) * avg(u[0]) * avg(v[0]) * self.geometry.dS(int(bound[1][0]))
                    self.L += Constant(bound[1][2]) * Constant(bound[1][1]) * avg(v[0]) * self.geometry.dS(int(bound[1][0]))
            if len(bound[1]) == 3:
                if bound[1][2] == 'ext':
                    self.a += Constant(bound[1][1]) * u[0] * v[0] * self.geometry.ds(int(bound[1][0])) \
                           - Constant(bound[1][1]) * u[lagrange_cnt+1] * v[0] * self.geometry.ds(int(bound[1][0])) \
                           - Constant(bound[1][1]) * (u[0] - u[lagrange_cnt+1]) * v[lagrange_cnt+1] * self.geometry.ds(int(bound[1][0]))
                elif bound[1][2] == 'int':
                    self.a += Constant(bound[1][1]) * avg(u[0]) * avg(v[0]) * self.geometry.dS(int(bound[1][0])) \
                           - Constant(bound[1][1]) * avg(u[lagrange_cnt+1]) * avg(v[0]) * self.geometry.dS(int(bound[1][0])) \
                           - Constant(bound[1][1]) * (avg(u[0]) - avg(u[lagrange_cnt+1])) * avg(v[lagrange_cnt+1]) * self.geometry.dS(int(bound[1][0]))
                lagrange_cnt += 1
            elif len(bound[1]) == 2:
                self.bcs.append(DirichletBC(self.geometry.funcspace['lagrange_mp'].sub(0), Constant(bound[1][1]), self.geometry.boundaries, bound[1][0]))
            else:
                #no BCs were applied, do nothing
                continue

    def apply_complex_bcs(self, u_r, u_i, v_r, v_i, omega):
        for bound in self.params.boundary_markers.items():
            if bound[1][2] == 'dispersive':
                interface_val = ct.double_layer_impedance(self.frequency)
                interface_val *= self.params.unit**2
            else:
                interface_val = bound[1][2] * self.params.unit**2
                interface_val = complex(interface_val.real, interface_val.imag*np.pi*2*self.frequency)
            if len(bound[1]) == 4:
                if bound[1][3] == 'ext':
                    a_r = Constant(np.real(interface_val)) * u_r * v_r * self.geometry.ds(int(bound[1][0])) \
                        - Constant(np.real(interface_val)) * u_i * v_i * self.geometry.ds(int(bound[1][0])) \
                        - Constant(np.imag(interface_val)) * u_i * v_r * self.geometry.ds(int(bound[1][0])) \
                        - Constant(np.imag(interface_val)) * u_r * v_i * self.geometry.ds(int(bound[1][0]))

                    a_i = Constant(np.real(interface_val)) * u_i * v_r * self.geometry.ds(int(bound[1][0])) \
                        + Constant(np.real(interface_val)) * u_r * v_i * self.geometry.ds(int(bound[1][0])) \
                        + Constant(np.imag(interface_val)) * u_i * v_i * self.geometry.ds(int(bound[1][0])) \
                        - Constant(np.imag(interface_val)) * u_i * v_i * self.geometry.ds(int(bound[1][0]))

                    self.a += a_r +a_i

                    L_r = Constant(np.real(interface_val)) * Constant(bound[1][1]) * v_r * self.geometry.ds(int(bound[1][0])) \
                          - Constant(np.imag(interface_val)) * Constant(bound[1][1]) * v_i * self.geometry.ds(int(bound[1][0]))

                    L_i = Constant(np.real(interface_val)) * Constant(bound[1][1]) * v_i * self.geometry.ds(int(bound[1][0])) \
                          + Constant(np.imag(interface_val)) * Constant(bound[1][1]) * v_r * self.geometry.ds(int(bound[1][0]))

                    self.L += L_r +L_i
                elif bound[1][3] == 'int':
                    a_r = Constant(np.real(interface_val)) * avg(u_r) * avg(v_r) * self.geometry.dS(int(bound[1][0])) \
                        - Constant(np.real(interface_val)) * avg(u_i) * avg(v_i) * self.geometry.dS(int(bound[1][0])) \
                        - Constant(np.imag(interface_val)) * avg(u_i) * avg(v_r) * self.geometry.dS(int(bound[1][0])) \
                        - Constant(np.imag(interface_val)) * avg(u_r) * avg(v_i) * self.geometry.dS(int(bound[1][0]))

                    a_i = Constant(np.real(interface_val)) * avg(u_i) * avg(v_r) * self.geometry.dS(int(bound[1][0])) \
                        + Constant(np.real(interface_val)) * avg(u_r) * avg(v_i) * self.geometry.dS(int(bound[1][0])) \
                        + Constant(np.imag(interface_val)) * avg(u_i) * avg(v_i) * self.geometry.dS(int(bound[1][0])) \
                        - Constant(np.imag(interface_val)) * avg(u_i) * avg(v_i) * self.geometry.dS(int(bound[1][0]))

                    self.a += a_r+a_i

                    L_r = - Constant(np.real(interface_val)) * Constant(bound[1][1]) * avg(v_r) * self.geometry.dS(int(bound[1][0])) \
                          + Constant(np.imag(interface_val)) * Constant(bound[1][1]) * avg(v_i) * self.geometry.dS(int(bound[1][0]))

                    L_i = - Constant(np.real(interface_val)) * Constant(bound[1][1]) * avg(v_i) * self.geometry.dS(int(bound[1][0])) \
                          - Constant(np.imag(interface_val)) * Constant(bound[1][1]) * avg(v_r) * self.geometry.dS(int(bound[1][0]))

                    self.L += L_r+L_i


    def apply_source(self, monopoles_dict):
        if monopoles_dict:
            for monopole in monopoles_dict['monopoles']:
                tmp_point = Point(*monopole[0:3])

                delta = PointSource(self.geometry.funcspace['scalar'], tmp_point, monopole[3])
                delta.apply(self.b)

    def apply_lagrange_mp_source(self, monopoles_dict):
        if monopoles_dict:
            for monopole in monopoles_dict['monopoles']:
                tmp_point = Point(*monopole[0:3])

                delta = PointSource(self.geometry.funcspace['lagrange_mp'].sub(0), tmp_point, monopole[3])
                delta.apply(self.b)

    def apply_complex_source(self, monopoles_dict):
        if monopoles_dict:
            for monopole in monopoles_dict['monopoles']:
                tmp_point = Point(*monopole[0:3])

                delta0 = PointSource(self.geometry.funcspace['mixed'].sub(0), tmp_point, monopole[3])
                delta0.apply(self.b)
                delta1 = PointSource(self.geometry.funcspace['mixed'].sub(1), tmp_point, monopole[3])
                delta1.apply(self.b)

    def get_poisson_values(self, coords, allsols = False, filename=None):

        if allsols or type(allsols)==int:
            values = []
            for idx, solution in enumerate(self.solutions):
                if (type(allsols)==int and idx==allsols) or (allsols == True and type(allsols)!=int):
                    if solution[1][0] == 'resistive':
                        projfunc = Function(self.geometry.funcspace['scalar'])
                    else:
                        if 'mixed' not in self.geometry.funcspace.keys():
                            P1 = FiniteElement("CG", self.geometry.mesh.ufl_cell(), self.order)
                            mesh_funcspace_mixed = FunctionSpace(self.geometry.mesh, MixedElement([P1,P1]))
                            self.geometry.funcspace['mixed'] = mesh_funcspace_mixed
                        projfunc = Function(self.geometry.funcspace['mixed'])
                    projfunc.vector().set_local(solution[0])
                    projfunc.set_allow_extrapolation(True)
                    values.append(self.loop_coords(projfunc, coords))
        else:
            self.phi.set_allow_extrapolation(True)
            values = self.loop_coords(self.phi, coords)
        if not os.path.exists('results'):
            os.makedirs('results')

        if filename:
            np.save(os.path.join('results',filename),values)

        return np.asarray(values)

    def loop_coords(self,phi, coords):
        values = []
        if len(np.shape(coords)) == 2:
            for idx in range(coords.shape[0]):
                try:
                    phi_tmp = phi(coords[idx,:])
                except RuntimeError:
                    phi_tmp = np.float64(0)
                if phi_tmp.size == 1: #self.phi.function_space().ufl_element().family() != 'Mixed':
                    values.append(phi_tmp)
                else:
                    values.append(phi_tmp[0]+phi_tmp[1]*1j)
        else:
            f = project(Constant(1.), self.geometry.funcspace['scalar'], solver_type="cg", preconditioner_type="amg")
            for idx in coords:
                try:
                    phi_tmp = assemble(phi*self.geometry.ds(idx))/assemble(f*self.geometry.ds(idx))
                except ZeroDivisionError:
                    phi_tmp = assemble(phi*self.geometry.dS(idx))/assemble(f*self.geometry.dS(idx))
                if  phi.vector().get_local().shape[0] == self.geometry.funcspace['scalar'].dim():
                    elec_boundary = DirichletBC(self.geometry.funcspace['scalar'], 1, self.geometry.boundaries,idx)
                    elec_func = Function(self.geometry.funcspace['scalar'])
                    elec_boundary.apply(elec_func.vector())
                    elec_values = phi.vector()[elec_func.vector() == 1]
                    tmpstd = np.std(elec_values)
                values.append([phi_tmp, elec_values, elec_values.min(), elec_values.max()])
        return values

    def get_error_values(self,ana_values,coords):

        fem_values = self.get_poisson_values(coords)
        diff = ana_values - fem_values[:,0]
        self.RMSE = np.sqrt(np.mean((diff)**2))
        self.RD = np.mean(np.abs(diff)/np.amax(np.abs(ana_values)))

        print( "RMSE ", self.RMSE, "RD", self.RD)

    def save_pvds(self, allsols = False, filename = 'fem_solution.pvd', subdomain_marker=None):
        '''
        save either the current (last solved) solution (default)
        or set allsolls = True and save all solutions of simulations
        manual filename=pvds_dir/string.pvd input
        automatic filename (based on params) is removed for now, will be added in future versions
        subdomain: if you want to export a region of the solution (e.g. only grey matter) give the marker as an input
        '''

        if not os.path.exists('pvds_dir'):
            os.makedirs('pvds_dir')

        if subdomain_marker:
            submesh = SubMesh(self.geometry.mesh, self.geometry.domain, subdomain_marker)

        if allsols:
            for idx, solution in enumerate(self.solutions):
                if filename[-4:] =='.pvd':
                    tmp_filename = filename[:-4]
                else:
                    tmp_filename = filename

                tmp_filename += str(idx)
                if tmp_filename[-4:]!='.pvd':
                    tmp_filename += '.pvd'

                if solution[1][0] == 'resistive':
                    projfunc = Function(self.geometry.funcspace['scalar'])
                    if subdomain_marker:
                        subV = FunctionSpace(submesh, "CG", self.order)
                else:
                    projfunc = Function(self.geometry.funcspace['mixed'])
                    if subdomain_marker:
                        P1 = FiniteElement("CG", submesh, self.order)
                        subV = FunctionSpace(submesh, MixedElement([P1,P1]))

                projfunc.vector().set_local(solution[0])
                if subdomain_marker:
                    subprojfunc = project(projfunc, subV)
                    File(os.path.join('pvds_dir',tmp_filename)) << subprojfunc
                else:
                    File(os.path.join('pvds_dir',tmp_filename)) << projfunc
        else:
            if filename[-4:]!='.pvd':
                filename += '.pvd'

            if subdomain_marker:
                if self.params.material == 'resistive':
                    subV = FunctionSpace(submesh, "CG", self.order)
                else:
                    P1 = FiniteElement("CG", submesh, self.order)
                    subV = FunctionSpace(submesh, MixedElement([P1,P1]))
                subprojfunc = project(self.phi, subV)
                File(os.path.join('pvds_dir',filename)) << subprojfunc
            else:
                File(os.path.join('pvds_dir',filename)) << self.phi

    def get_all_vert_normals(self, domain, domain_marker, boundary_marker):
        '''
        Get normals for vertices in a domain (domain meshfunction, domain_marker) relative to the closest boundary (boundary_marker)
        '''

        tic = time()
        print("If a large domain/many points are needed, this might take a while.")
        It_domain = SubsetIterator(domain, domain_marker)
        pointlist = []
        for c in It_domain:
            for f in facets(c):
                for v in vertices(f):
                    pointlist.append([v.point().x(),v.point().y(),v.point().z(),v.index()])
        pointlist = np.asarray(pointlist)
        self.unique_index, idx_of_idx = np.unique(pointlist[:,3], return_index=True)
        self.point_coords_list = pointlist[idx_of_idx,0:3]

        #self.dists_indices, self.normals = self.get_closest_boundary(boundary_marker, self.point_coords_list)
        self.get_closest_boundary(boundary_marker, self.point_coords_list)
        print("Ready after ", (time()-tic)/60, " minutes")


    def get_closest_boundary(self, boundary_marker, location_point_coords):
        '''
        input boundary marker to compare point(s) to and list with point(s) of interest
        returns the facet within the given boundary that is closest to the location_point and its normal
        '''
        tic = time()
        It_boundaries = SubsetIterator(self.geometry.boundaries, boundary_marker)

        bound_midpoints = []
        for f in It_boundaries:
            bound_midpoints.append([f.midpoint().x(),f.midpoint().y(),f.midpoint().z(),f.index()])
        bound_midpoints = np.asarray(bound_midpoints)

        distances = np.sqrt(np.sum((location_point_coords[:,None]-bound_midpoints[:,0:3])**2, axis=2))
        closest_args = np.argmin(distances,axis=1)
        self.distances=distances

        self.closest_normals = []
        for idx in range(location_point_coords.shape[0]):
            closest_normal = Facet(self.geometry.mesh, int(bound_midpoints[closest_args[idx], 3])).normal()
            self.closest_normals.append(closest_normal)

        '''
        dist2loc = []
        for f in It_boundaries:
            for location_point in location_point_coords:
                if type(location_point) != dolfin.cpp.mesh.Point:
                    location_point = Point(location_point)
                dist2loc.append([location_point.distance(f.midpoint()), f.index()])

        dist2loc = np.asarray(dist2loc)
        dist2loc = dist2loc.reshape([int(dist2loc.shape[0]/location_point_coords.shape[0]), location_point_coords.shape[0], 2])
        closest_args = np.argmin(dist2loc[:,:,0], axis=0)
        closest_normals = []
        for idx in range(location_point_coords.shape[0]):
            closest_normal = Facet(self.geometry.mesh, int(dist2loc[closest_args[idx], idx,1])).normal()
            closest_normals.append(closest_normal)

        return dist2loc[closest_args,range(location_point_coords.shape[0]),:], closest_normals
        '''

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
            print("Less locations than requested created, adjust distance or inner_distance if more are needed.")

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
                    print("name types of orientations list of sources can be 'mono', 'rad' or 'tan'")

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
                    self.params.volume_markers[vol[0]] = tissue
                    found = True
                    break
            if not found:
                print('could not automatically set ', vol[0], ' to dispersive setting \n')
                print('please do this manually, options can be found here: help(disp.get_gabriel_params)')

    def make_pulse(self,time,width,npoints):
        Fs = npoints/time #sampling frequency
        self.params.t = np.arange(-time/2.,time/2.,1./Fs)
        self.params.x = np.zeros(len(self.params.t))
        self.params.x[int(len(self.params.t)/2-width*Fs/2.)+1:int(len(self.params.t)/2+width*Fs/2.)+1] = 1
        self.params.X = np.fft.fft(self.params.x)
        self.params.frequencies = np.fft.fftfreq(len(self.params.t), d=1./Fs)

    def make_cos(self,time,npoints):
        Fs = npoints/time #sampling frequency
        self.params.t = np.arange(-time/2.,time/2.,1./Fs)
        self.params.x = np.cos(2.*np.pi*(np.arange(-npoints/2.,npoints/2.)/npoints))
        self.params.X = np.fft.fft(self.params.x)
        self.params.frequencies = np.fft.fftfreq(len(self.params.t), d=1./Fs)

    def make_alphafunction(self,time,npoints):
        tau = 1e-3 #1ms
        self.params.t = np.linspace(0,time,npoints)
        self.params.x = self.params.t/tau * np.exp(-self.params.t/tau)
        Fs = npoints/time
        self.params.X = np.fft.fft(self.params.x)
        self.params.frequencies = np.fft.fftfreq(len(self.params.t), d=1./Fs)

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

        self.funcspace = init_funcspace(self.mesh, fem_sim.order)


    def init_admittivity(self,fem_sim):
        # make sigma and if needed epsilon function
        self.sigma = Function(self.funcspace['const'])
        self.conductivities = np.empty(len(fem_sim.params.volume_markers),dtype=float)

        if ((fem_sim.params.material == "capacitive") or (fem_sim.params.material == "dispersive")) and not hasattr(self, 'epsilon'):
            self.epsilon = Function(self.funcspace['const'])
            self.permittivities = np.empty(len(fem_sim.params.volume_markers),dtype=float)

        #anisotropy tensor
        try:
            if list(fem_sim.params.volume_markers.items())[0][1].shape[1] == 9:
                aniso_matrix = []
                for ii in range(9):
                    tmp_DGfunc = Function(self.funcspace['const'])
                    tmp_DGfunc.vector().set_local(list(fem_sim.params.volume_markers.items())[0][1][:,ii])
                    aniso_matrix.append(tmp_DGfunc)
                self.sigma = as_matrix(((aniso_matrix[0], aniso_matrix[1], aniso_matrix[2]), (aniso_matrix[3], aniso_matrix[4], aniso_matrix[5]), (aniso_matrix[6], aniso_matrix[7], aniso_matrix[8])))
                return
        except:
            pass

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

    def cells_in_region(self, distance, point, minimum_inradius):

        refine_markers = MeshFunction('bool',self.mesh,3)
        refine_markers.set_all(False)

        for c in cells(self.mesh):
            if c.midpoint().distance(point) < distance and c.inradius() > minimum_inradius:
                refine_markers[c] = True
        return refine_markers

    def refine_sources(self, pointsources, distance, minimum_inradius, femsim):

        '''
        refine to a minimum cell inradius in a region around points
        pointsources: list with dolfin Point instances around which refinement is done
        distance: float distance from points where minimum inradius is needed
        minimum_inradius: desired inradius
        femsim: relevant (current) FEMsim class
        '''

        print("WARNING: this will overwrite the mesh and functionspaces in the FEM_geometry class. \nSolutions with the original mesh can't be projected to the right (unrefined) functionspace anymore")
        print("This will be fixed in future versions")

        parameters["refinement_algorithm"] = "plaza_with_parent_facets"


        orig_cells = self.mesh.num_cells()
        while True:
            sum_markers = 0
            for pointsource in pointsources:
                refine_markers = self.cells_in_region(distance, pointsource, minimum_inradius)

                refined_mesh = refine(self.mesh, refine_markers)
                self.mesh = refined_mesh


                refined_bounds = adapt(self.boundaries, self.mesh)
                self.boundaries = refined_bounds
                refined_domain = adapt(self.domain, self.mesh)
                self.domain = refined_domain
                for functype in self.funcspace:
                    adapt(self.funcspace[functype], self.mesh)
                    self.funcspace[functype] = self.funcspace[functype].child()


                sum_markers += np.sum(refine_markers.array())

            if sum_markers == 0:
                break

        # refined boundaries not on a marked boundary need to be manually set to zero
        self.boundaries.array()[self.boundaries.array()[:]>1e10] = 0

        print("Number of cells increase from ", orig_cells, " to ", self.mesh.num_cells(), " cells")

        self.init_admittivity(femsim)

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

    print("Computing boundaries...")

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
    print( "Computing boundaries took", (time()-tic)/60., "minutes")
    return boundaries

def get_pots_near_elec(mesh, solution, distpoint, radius, sourcepoint, thresh):
    coords = []
    distances = []
    for c in cells(mesh):
        pnt = Point(c.midpoint().x(), c.midpoint().y(), c.midpoint().z())
        if distpoint.distance(pnt) < radius and distpoint.distance(sourcepoint) > thresh:
            coords.append([c.midpoint().x(), c.midpoint().y(), c.midpoint().z()])
            distances.append(distpoint.distance(pnt))
    coords = np.asarray(coords)
    distances = np.asarray(distances)
    pots = []
    for ii in range(len(coords)):
        pots.append(solution(coords[ii]))
    pots = np.asarray(pots)
    dist_sort = np.argsort(distances)
    return pots[dist_sort], coords[dist_sort], distances[dist_sort]
