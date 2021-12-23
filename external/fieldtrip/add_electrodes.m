function [newsurface, elecmarkers] = add_electrodes(surface, elec_coords, rad_elec, depth_elec, elec_normals, elec_shape)

% ADD_ELECTRODES creates electrodes in a volumetric head segmentation
% or triangulated surface and returns a tetrahedral mesh which is the input
% for a volume conduction model using the finite element method in FEMfuns
% Currently only cylinder electrodes (default) can be added.
% NOTE: no checks are done (yet) to see if the electrodes will self-intersect with each other, 

if ~exist('elec_shape', 'var') || isempty(elec_shape)
    elec_shape = 'cylinder';
end

if ~exist('elec_normals', 'var') || isempty(elec_normals)
    % find the orientation of the electrodes based on the nearest triangle
    elec_normals = electrode_normal(surface, elec_coords);
end

%create triangulated electrodess surfaces
for i = 1:size(elec_coords,1)
    if strcmp(elec_shape,'cylinder')
        electrodesurface = cylinder_electrode(elec_coords(i,:), elec_normals(i,:), rad_elec, depth_elec);
    else
        error('Electrode shape %s can not be used.', elec_shape)
    end
    savetri(['elec' num2str(i) '.tri'],electrodesurface.pos,electrodesurface.tri(:,1:3))
end

% find the location of the trident binary
str = which('trident.glnx86');
[p, f, x] = fileparts(str);
trident = fullfile(p, f);

% the following determines the trident executable to use
switch lower(computer)
  case {'maci' 'maci64'}
    % apple computer
    trident = [trident '.maci'];
  case {'glnx86' 'glnxa64'}
    % linux computer
    trident = [trident '.glnx86'];
  case {'win32', 'win64'}
    % windows computer
    trident = [trident '.exe'];
  otherwise
    ft_error('there is no trident executable for your platform');
end

% write out surface where electrodes will be added (currently, intended to for inner skull)
savetri('dent_surf.tri',surface.pos,surface.tri)

% indent surface (inner skull) with electrodes

for i = 1:size(elec_coords,1)
    cmd = sprintf('%s elec%i.tri dent_surf.tri elec%i_out.tri dent_surf.tri', trident, i, i);

    status = system(cmd);
    if status ~= 0
        error('Cylinder for electrode %g could not be pushed in inner skull',i)
    end
end

% combine the surface with the electrodes
[newsurface.pos, newsurface.tri] = loadtri('dent_surf.tri');
elecmarkers = [];
for i = 1:size(elec_coords,1)
    [elec_surf.pos,elec_surf.tri] = loadtri(['elec' num2str(i) '_out.tri']);
    newsurface = add_surf(newsurface,elec_surf,0);

    pointinsurf = point_in_surf(elec_surf);
    elecmarkers = [elecmarkers; pointinsurf];
end

% clean up the tri files
system('rm dent_surf.tri');
for i = 1:size(elec_coords,1)
    system(['rm elec' num2str(i) '.tri']);
    system(['rm elec' num2str(i) '_out.tri']);
end
