% This function creates a tetrahedral two-sphere mesh with cylindrical electrodes lying on top of the inner sphere (i.e., intracranial electrodes)
% which can be used to perform (e.g., ECoG) forward simulations, particularly in FEMfuns (Finite Element Method for useful neuroscience simulations, https://github.com/meronvermaas/FEMfuns).
% It has been created in the IntoTheBrain project and the functions can be used to use FEMfuns forward solutions in FieldTrip.
% The structure of this script is more or less:
% 1. Create two triangulated spheres representing the brain and skull and electrode coordinates on top of the brain surface (FieldTrip workflow)
% 2. Create cylinders at each electrode coordinate. (embedded FEMfuns-FieldTrip workflow)
%   a. make the cylinder surface crossing through the spherical surface it is combined with
%   b. refine the spherical surface where the cylinder crosses it
%   c. cut off the part of the cylindrical surface sticking into the sphere
% 3. Combine all the surfaces into one surface with unique faces and nodes (embedded FEMfuns-FieldTrip workflow)
% 4. Make a tetrahedral mesh struct (embedded FEMfuns-FieldTrip workflow)
% 5. Create a sourcemodel (FieldTrip workflow)
% 6. Create the leadfield (embedded FEMfuns-FieldTrip workflow)
%   a. Export parameters and mesh geometry
%   b. Run FEMfuns
%   b. Import forward solution

addpath('/PATHTO/fieldtrip')

ft_defaults;
% this requires the external iso2mesh toolbox
ft_hastoolbox('iso2mesh', 1);

% Create a spherical volume conductor with two spheres of radius 7 and 10 cm at the origin
csvol.o = [0, 0,0];
csvol.r = [7 10];
cfg = [];
cfg.numvertices = 1000;
csbnd = ft_prepare_mesh(cfg, csvol);

%pick a few electrode positions on top of the sphere
sel = find(csbnd(1).pos(:,3)>0); sel = sel(1:100:end);
elec = [];
elec.elecpos = csbnd(1).pos(sel,:);
for i=1:length(sel)
  elec.label{i} = sprintf('elec%d', i);
end
elec.unit = 'cm';
% update the electrode sets to the latest standards
elec = ft_datatype_sens(elec);

% combine the inner skull surface (i.e., brain sphere) with electrode surfaces and get inner points of the electrodes
dp_elec = 0.5; %height  of the electrode cylinder
rad_elec = 0.2; %radius of the electrode cylinder
[dented_elsurf,elecmarkers] = add_electrodes(csbnd(1), elec.elecpos, rad_elec, dp_elec);
merged_surfs = add_surf(dented_elsurf,csbnd(2));

%create volumetric tetrahedral mesh
[tet_node,tet_elem] = s2m(merged_surfs.pos,merged_surfs.tri, 1, 1, 'tetgen', [point_in_surf(csbnd(1));point_in_surf(csbnd(2));elecmarkers]);
%label the electrode surface where they make contact with the brain
el_faces = label_surf(tet_elem, 3:length(elec.elecpos)+2, 1);

%Construct the FT mesh structure
mesh.unit = 'cm';
mesh.pos = tet_node;
mesh.tet = tet_elem(:,1:4);
mesh.tri = el_faces(:,1:3);
mesh.boundary = el_faces(:,4);
mesh.boundarylabel = elec.label;
mesh.tissue = tet_elem(:,5);
mesh.tissuelabel = [{'brain'}, {'skull'},elec.label(:)'];

% construct a vol to create the FT sourcemodel
vol.pos = mesh.pos;
vol.tet = mesh.tet;
vol.tissue = mesh.tissue;
vol.tissuelabel = mesh.tissuelabel;
vol.unit = mesh.unit;
vol.type = 'simbio';

cfg                 = [];
cfg.resolution      = 5; %in mm
cfg.headmodel       = vol;
cfg.inwardshift     = 1; %shifts dipoles away from surfaces
sourcemodel         = ft_prepare_sourcemodel(cfg);

% conductivities for brain, scalp and metal electrodes are set
conductivities = [0.33 0.01 1e10 1e10 1e10 1e10 1e10];
lf_rec = femfuns_leadfield(mesh,conductivities,sourcemodel,elec);

% Instead dipole sources, a stimulating and ground electrode is set.
% For boundary options, look for example here https://github.com/meronvermaas/FEMfuns/blob/master/separability/parameters_discelecins.py
sourcemodel.inside(:) = false;
elec.label{1} = ['stim_' elec.label{1}];
elec.label{2} = ['ground_' elec.label{2}];
mesh.boundarylabel = elec.label;
elec.params{1} = {1, 100, 'int'}; %stimulating electrode (e.g. 1mV) with surface conductance
elec.params{2} = {0, 100, 'int'};
elec.params{3} = {100, 'int'}; %recording electrodes with surface conductance
elec.params{4} = {100, 'int'};
elec.params{5} = {100, 'int'};
conductivities = [0.33 0.01 0 0 0 0 0];
lf_stim = femfuns_leadfield(mesh,conductivities,sourcemodel,elec);
