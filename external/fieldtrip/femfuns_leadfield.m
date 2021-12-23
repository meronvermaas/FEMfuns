function [grid] = femfuns_leadfield(mesh,conductivities,sourcemodel,elec,frequencies,removetmp)
% FEMFUNS_LEADFIELD computes the forward solution using FEMfuns and returns a FieldTrip leadfield struct.
% Inputs:
% - a FieldTrip mesh struct (with both tetrahedra and optionally boundaries)
% - conductivities array
% - FieldTrip sourcemodel
% - FieldTrip elec struct 
% Which are written as a temporary parameter file that is used by FEMfuns.
% An example solution is written out in a pvds_dir for inspection.
% Allowing for capacitive and dispersive properties is a work in progress.

if ~exist('frequencies', 'var') && (~isreal(conductivities) || ischar(conductivities) || isstring(conductivities))
    error('Capacitive or dispersive materials require frequencies input');
end
    
%assert(~(~exist('frequencies', 'var')&&(~isreal(conductivities) || ischar(conductivities) || isstring(conductivities))), 'Capacitive or dispersive materials require frequency input');

fileID = fopen('params.py','w');
fprintf(fileID,'import numpy as np\n\n');
if isreal(conductivities) && ~(ischar(conductivities) || isstring(conductivities))
    fprintf(fileID,'material = "resistive"\n\n');
elseif ~isreal(conductivities)
    fprintf(fileID,'material = "capacitive"\n\n');
elseif ischar(conductivities) || isstring(conductivities)
    fprintf(fileID,'material = "dispersive"\n\n');
    error('Dispersive tissue has not been fully implemented yet')
end

if strcmp(mesh.unit,'m')
    unit = 1;
elseif strcmp(mesh.unit,'cm')
    unit = 1e-2;
elseif strcmp(mesh.unit,'mm')
    unit = 1e-3;
elseif strcmp(mesh.unit,'um')
    unit = 1e-6;
else
    error("This unit is not supported, m, cm, mm and um can be used")
end
fprintf(fileID,'unit = %1.2e\n\n',unit);

%TODO check if length conductivites is correct / at least as long as number
%of labels
fprintf(fileID,'volume_markers = {\n');
vol_mark = '';
for ii = 1:length(mesh.tissuelabel)
    vol_mark = sprintf('%s''%s'': [%d, %3.2f],\n',vol_mark,mesh.tissuelabel{ii}, ii, conductivities(ii));
end
vol_mark = vol_mark(1:end-2);
fprintf(fileID,'%s\n}\n\n',vol_mark);

% TODO more boundary condition option in should be included here
if isfield(mesh,'boundarylabel')
    fprintf(fileID,'boundary_markers = {\n');
    boundvals = unique(mesh.boundary);
    bnd_mark = '';
    %%%%%%TODO add functionality and check if elec_props is given with surface_conductances and int/ext with
    %%%%%%optional potential for stimulating electrode potential
    if isfield(elec,'params') && length(elec.params) < length(mesh.boundarylabel)
        for ii = length(elec.params)+1:length(mesh.boundarylabel)
            elec.params{ii} = {};
        end
    end
    for ii = 1:length(mesh.boundarylabel)
        if isfield(elec,'params') && ~isempty(elec.params{ii})
            bnd_mark = sprintf('%s''%s'': [%d ',bnd_mark,mesh.boundarylabel{ii}, boundvals(ii));
            for jj = 1:length(elec.params{ii})
                if isa(elec.params{ii}{jj},'double')
                    bnd_mark = sprintf('%s, %f',bnd_mark,elec.params{ii}{jj});
                elseif isa(elec.params{ii}{jj},'char')
                    bnd_mark = sprintf('%s, ''%s''',bnd_mark,elec.params{ii}{jj});
                else
                    error('electrode params should be a double or char parameter')
                end
            end
            bnd_mark = sprintf('%s],\n',bnd_mark);
        else
            bnd_mark = sprintf('%s''%s'': [%d],\n',bnd_mark,mesh.boundarylabel{ii}, boundvals(ii));
        end
    end
    bnd_mark = bnd_mark(1:end-2);
    fprintf(fileID,'%s\n}\n\n',bnd_mark);
else
    %just fill the field
    fprintf(fileID,'boundary_markers = {\n}\n\n');
end

% write dipole sources as 3 bipoles in x, y and z
fprintf(fileID,'monopole_list = [\n');
monopolelist = '';
dippos = sourcemodel.pos(sourcemodel.inside,:);
for ii = 1:size(dippos,1)
    monopole1 = sprintf('%f,' , dippos(ii,:)+[1e-2 0 0]);
    monopole2 = sprintf('%f,' , dippos(ii,:)-[1e-2 0 0]);
    monopolelist = sprintf('%s{''monopoles'': [[%s 1],[%s -1]], ''name'': ''%s''},\n',monopolelist,monopole1,monopole2, ['xpole' num2str(ii)]);
    
    monopole1 = sprintf('%f,' , dippos(ii,:)+[0 1e-2 0]);
    monopole2 = sprintf('%f,' , dippos(ii,:)-[0 1e-2 0]);
    monopolelist = sprintf('%s{''monopoles'': [[%s 1],[%s -1]], ''name'': ''%s''},\n',monopolelist,monopole1,monopole2, ['ypole' num2str(ii)]);
    
    monopole1 = sprintf('%f,' , dippos(ii,:)+[0 0 1e-2]);
    monopole2 = sprintf('%f,' , dippos(ii,:)-[0 0 1e-2]);
    monopolelist = sprintf('%s{''monopoles'': [[%s 1],[%s -1]], ''name'': ''%s''},\n',monopolelist,monopole1,monopole2, ['zpole' num2str(ii)]);
end
monopolelist = monopolelist(1:end-2);
if isempty(monopolelist)
    monopolelist = 'None';
end
fprintf(fileID,'%s\n]\n\n',monopolelist);

ellist = '';
for ii = 1:numel(elec.elecpos)/3
    if isfield(elec,'params') && isempty(elec.params{ii})
        ellist = sprintf('%s[',ellist);
        
        ellist = [ellist sprintf('%f,',elec.elecpos(ii,:))];
        ellist = ellist(1:end-1);
        ellist = sprintf('%s],',ellist);
    end
    
end
ellist = ellist(1:end-1);
fprintf(fileID, 'elecpos = np.array([%s])',ellist);

if exist('frequencies', 'var')
    fprintf(fileID, 'frequencies = np.array([%s])', sprintf('%f, ',frequencies));
    error('This is a work in progress and will be available in the next version (soon)')
end

fclose(fileID);

node = mesh.pos;
elem = [mesh.tet mesh.tissue];
bnd = [mesh.tri mesh.boundary];
meshfilename = "tmp_export.msh";
savemsh_femfuns_bnd(node,elem,bnd,meshfilename)

outmeshfilename = "tmp_export.xml";
cmd = sprintf("sh femfuns_caller.sh %s %s",meshfilename,outmeshfilename);
status = system(cmd);
assert(status == 0, 'Running FEMfuns forward computation did not succeed')

if exist('removetmp','var') && removetmp
    %remove tmp files
    system('rm tmp_export* params.py');
end

lf = load('femfuns_lfmatrix.mat');

grid.dim = sourcemodel.dim;
grid.pos = sourcemodel.pos;
grid.unit = sourcemodel.unit;
grid.inside = sourcemodel.inside;
grid.cfg = sourcemodel.cfg;

if sum(sourcemodel.inside) > 0
    grid.leadfield = cell(1,size(grid.pos,1));
    lf = reshape(lf.lf,3,[],size(elec.elecpos,1));
    inside_idxs = find(grid.inside);
    for ii = 1:length(inside_idxs)
        grid.leadfield{inside_idxs(ii)} = squeeze(lf(:,ii,:))';
    end
else
    grid.leadfield = lf.lf;
end

grid.label = elec.label;

end