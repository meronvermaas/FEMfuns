% This function creates a tetrahedral mesh of the brain, skull and scalp with cylindrical electrodes lying on top of the inner sphere (i.e., intracranial electrodes)
% which can be used to perform (e.g., ECoG) forward simulations, particularly in FEMfuns (Finite Element Method for useful neuroscience simulations, https://github.com/meronvermaas/FEMfuns).
% It has been created in the IntoTheBrain project and the functions can be used to use FEMfuns forward solutions in FieldTrip.
% The structure of this script is more or less:
% 1. Segment the anatomical mri and create surfaces (FieldTrip workflow, using spm12 and iso2mesh)
% 2. Create cylinders at each electrode coordinate. (embedded FEMfuns-FieldTrip workflow, using trident)
%   a. make the cylinder surface crossing through the brain surface it is combined with
%   b. refine the brain surface where the cylinder crosses it
%   c. cut off the part of the cylindrical surface sticking into the sphere
% 3. Combine all the surfaces into one surface with unique faces and nodes (embedded FEMfuns-FieldTrip workflow)
% 4. Make a tetrahedral mesh struct (embedded FEMfuns-FieldTrip workflow, using iso2mesh)
% 5. Create a sourcemodel (FieldTrip workflow)
% 6. Create the leadfield (embedded FEMfuns-FieldTrip workflow)
%   a. Export parameters and mesh geometry
%   b. Run FEMfuns
%   b. Import forward solution

addpath('/PATHTO/fieldtrip')

ft_defaults;
% this requires the external iso2mesh toolbox
ft_hastoolbox('iso2mesh', 1);

% Read in the MRI data
mri = ft_read_mri('Subject01.mri'); % the anatomical mri can be found at FieldTrip: ftp://ftp.fieldtriptoolbox.org/pub/fieldtrip/tutorial/Subject01.zip

% reslice the volume to ensure homogenous voxels
cfg     = [];
cfg.dim = mri.dim;
mri     = ft_volumereslice(cfg,mri);

% segment the anatomical mri into 3 different tissue types 
cfg           = [];
cfg.output    = {'brain','skull','scalp'};
segmentedmri  = ft_volumesegment(cfg, mri);

% create surfaces at the borders of the different tissue-types
cfg        = [];
cfg.tissue={'brain','skull','scalp'};
cfg.numvertices = 3000;
cfg.method = 'iso2mesh';
bnd = ft_prepare_mesh(cfg,segmentedmri);

% get coordinates lying within these surfaces (needed for labeling)
insidepoints = zeros(length(bnd),3);
for ii = 1:length(bnd)
    insidepoints(ii,:) = point_in_surf(bnd(ii));
end

% pick 4 electrodes on top of the brain surface and create the fieldtrip electrode structure
el1 = [51.298 -26.242 106.026]; el2 = [38.926 -37.6445 105.965]; el3 = [31.47 -28.13 113.007]; el4 = [45.06 -13.467 112.423];
[~,I1] = min(abs(sum(bnd(1).pos-el1,2))); [~,I2] = min(abs(sum(bnd(1).pos-el2,2))); [~,I3] = min(abs(sum(bnd(1).pos-el3,2))); [~,I4] = min(abs(sum(bnd(1).pos-el4,2)));
sel = [I1 I2 I3 I4];
elec = [];
elec.elecpos = bnd(1).pos(sel,:);
for i=1:length(sel)
  elec.label{i} = sprintf('elec%d', i);
end
elec.unit = 'mm';
% update the electrode sets to the latest standards
elec = ft_datatype_sens(elec);

% combine the brain surface with electrode surfaces and get inner points of the electrodes
dp_elec = 0.5; %height  of the electrode cylinder
rad_elec = 2; %0.2; %radius of the electrode cylinder
[dented_elsurf,elecmarkers] = add_electrodes(bnd(1), elec.elecpos, rad_elec, dp_elec);

% Add the skull and scalp to the brain again
% TODO this should be done automatically in add_electrodes()
merged_surfs = dented_elsurf;
for ii = 2:length(bnd)
    merged_surfs = add_surf(merged_surfs,bnd(ii));
end

%create volumetric tetrahedral mesh
[tet_node,tet_elem] = s2m(merged_surfs.pos,merged_surfs.tri, 1, 1, 'tetgen', [insidepoints; elecmarkers]);

%label the electrode surface where they make contact with the brain
el_faces = label_surf(tet_elem, length(bnd)+1:length(elec.elecpos)+length(bnd), 1);

%Construct the FT mesh structure
mesh.unit = 'cm';
mesh.pos = tet_node;
mesh.tet = tet_elem(:,1:4);
mesh.tri = el_faces(:,1:3);
mesh.boundary = el_faces(:,4);
mesh.boundarylabel = elec.label;
mesh.tissue = tet_elem(:,5);
mesh.tissuelabel = [{'brain'}, {'skull'},{'scalp'},elec.label(:)'];

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
conductivities = [0.33 0.01 0.3 1e10 1e10 1e10 1e10 1e10];
lf_rec = femfuns_leadfield(mesh,conductivities,sourcemodel,elec);