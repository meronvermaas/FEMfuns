function savemsh_femfuns_bnd(node,elem,bnd,fname,rname,bname)
%adjust based on correct gmsh format
%https://www.manpagez.com/info/gmsh/gmsh-2.2.6/gmsh_63.php
%
% savemsh_femfuns(node,elem,fname,rname)
%
% save a tetrahedral mesh to GMSH mesh format (version that works for FEniCS)
%
% input:
%      node: input, node list, dimension (nn,3)
%      elem: input, tetrahedral mesh element list, dimension (ne,4) or (ne,5) for multi-region meshes
%      fname: output file name
%      rname: name of the regions, cell-array of strings (optional)
%
% -- this function is adapted from iso2mesh toolbox (http://iso2mesh.sf.net)
%

if (~exist('rname', 'var')) , rname = {}; end
if (~exist('bname', 'var')) , bname = {}; end
if size(elem,2) < 5 , elem(:,5) = 1 ; end
if size(bnd,2) < 4 , bnd(:,4) = 1 ; end
if any(bnd(:,4) == 0), bnd(bnd(:,4)==0,:) = []; end
%if any(bnd(:,4) == 0), bnd(:,4) = bnd(:,4) + 1 ; end
fid = fopen(fname,'wt');
if(fid==-1)
    error('You do not have permission to save mesh files.');
end
nbNodes = size (node,1);
reg = unique (elem(:,5));
bnds = unique(bnd(:,4));
nbRegion = length(reg);
nbBoundary = length(bnds);
nbElementsTet = size (elem, 1);
nbElementsBnd =  size(bnd,1);

% Create the skeleton of the mesh structure
M.Info.version = [];

M.Nodes.nb = 0;
M.Nodes.x = [];
M.Nodes.y = [];
M.Nodes.z = [];

M.Elements.nb = 0;
M.Elements.type = zeros (0, 0, 'uint8');
M.Elements.tableOfNodes = zeros (0, 0, 'uint32');
M.Elements.regionTet = zeros (0, 0, 'uint16');
M.Elements.regionBnd = zeros (0, 0, 'uint16');

M.Boundaries.nb = 0;
M.Boundaries.name = {};
M.Boundaries.dimension = [];

M.Regions.nb = 0;
M.Regions.name = {};
M.Regions.dimension = [];

% Build the table of nodes

M.Nodes.nb = nbNodes;
M.Nodes.x = node(:,1);
M.Nodes.y = node(:,2);
M.Nodes.z = node(:,3);
clear node

% Build the table of elements

M.Elements.nb = nbElementsTet + nbElementsBnd;
M.Elements.type = [uint8(2*ones(nbElementsBnd, 1)); uint8(4*ones(nbElementsTet, 1))];
M.Elements.tableOfNodesTet = uint32(elem(:,1:4));
M.Elements.tableOfNodesBnd = uint32(bnd(:,1:3));
M.Elements.regionTet = uint16(elem(:,5));
M.Elements.regionBnd = uint16(bnd(:,4));
clear elem bnd

% Build the table of regions and boundaries

M.Boundaries.nb = max(bnds);
for k = 1 : nbBoundary
   if length(bname) < k , bname{k} = sprintf('boundary_%d', k) ; end
   M.Boundaries.name{bnds(k)} = sprintf ('%s', bname{k});
   M.Boundaries.dimension(bnds(k)) = 2;
end

M.Regions.nb = max(reg);
for k = 1 : nbRegion
   if length(rname) < k , rname{k} = sprintf('region_%d', k) ; end
   M.Regions.name{reg(k)} = sprintf ('%s', rname{k});
   M.Regions.dimension(reg(k)) = 3;
end

% Writhe the header
fprintf (fid, '$MeshFormat\n2.2 0 8\n$EndMeshFormat\n');

% % % Write the physical names
% % if M.Regions.nb > 0 && M.Boundaries.nb >0
% %    fprintf (fid, '$PhysicalNames\n');
% %    fprintf (fid, '%d\n', M.Regions.nb+M.Boundaries.nb);
% %    
% %    for r = 1 : M.Boundaries.nb
% %       name = M.Boundaries.name{r};
% %       if isempty (name)
% %          name = sprintf ('Boundary_%d', r);
% %       end
% %       fprintf (fid, '%d %d "%s"\n', M.Boundaries.dimension(r), r, name);
% %    end
% %    
% %    for r = 1 : M.Regions.nb
% %       name = M.Regions.name{r};
% %       if isempty (name)
% %          name = sprintf ('Region_%d', r);
% %       end
% %       fprintf (fid, '%d %d "%s"\n', M.Regions.dimension(r), r, name);
% %    end
% %    fprintf (fid, '$EndPhysicalNames\n');
% % end

% Write the nodes
fprintf (fid, '$Nodes\n');
fprintf (fid, '%d\n', size(M.Nodes.x,1));
buffer = [ 1:M.Nodes.nb ; M.Nodes.x' ; M.Nodes.y' ; M.Nodes.z' ];
fprintf (fid, '%d %10.10f %10.10f %10.10f\n', buffer);
fprintf (fid, '$EndNodes\n');

% Write the elements
%
% In order to accelerate the printing, the elements are grouped by (homogeneous) type.
%

fprintf (fid, '$Elements\n');
fprintf (fid, '%d\n', M.Elements.nb);

type = unique (M.Elements.regionBnd);
tmp = 0;
for k = 1 : length(type)
    et = find(M.Elements.regionBnd == type(k));
    buffer = zeros (8 , length(et));
    buffer(1,:) = tmp + 1 : tmp+length(et);
    buffer(2,:) = 2; %3-node triangle element https://www.manpagez.com/info/gmsh/gmsh-2.2.6/gmsh_63.php
    buffer(3,:) = 2;
    buffer(4,:) = M.Elements.regionBnd(et);
    buffer(5,:) = M.Elements.regionBnd(et); %telt door vanaf 1
    for n = 1 : 3
          buffer(5+n,:) = M.Elements.tableOfNodesBnd(et,n);
    end
    tmp = tmp+length(et);
    
    elementFormat = '%d %d %d %d %d %d %d %d \n';
    fprintf (fid, elementFormat, buffer);
end

type = unique (M.Elements.regionTet);
for k = 1 : length(type)
    et = find(M.Elements.regionTet == type(k));
    buffer = zeros (9 , length(et));
    buffer(1,:) = tmp + 1 : tmp+length(et);
    buffer(2,:) = 4; %4-node tetrahedron element https://www.manpagez.com/info/gmsh/gmsh-2.2.6/gmsh_63.php
    buffer(3,:) = 2;
    buffer(4,:) = M.Elements.regionTet(et);
    buffer(5,:) = M.Elements.regionTet(et);
    for n = 1 : 4
          buffer(5+n,:) = M.Elements.tableOfNodesTet(et,n);
    end
    tmp = tmp+length(et);
    
    elementFormat = '%d %d %d %d %d %d %d %d %d \n';
    fprintf (fid, elementFormat, buffer);
end

% % % for h = 1 : blockSize : ceil(length(M.Elements.type)/blockSize)*blockSize
% % %    e = h : min(length(M.Elements.type) , h+blockSize-1);  % = elements being considered
% % %    type = unique (M.Elements.type(e));                    % = types of elements found in this group
% % %    
% % %    %
% % %    % Process each type of element separately
% % %    %
% % %    for k = 1 : length(type)
% % %       if type(k) == 0, continue; end
% % %       et = e(find(M.Elements.type(e) == type(k)));  % = elements of the group of the same type
% % %       
% % %       %
% % %       % Determine the format for printing the elements
% % %       %
% % %       elementFormat = '%d %d %d %d %d %d'; %\n';
% % %       for n = 1 : 4
% % %          elementFormat = [ elementFormat '%d ' ];
% % %       end
% % %       elementFormat = [ elementFormat '\n' ];
% % %       
% % %       %
% % %       % Collect in a buffer all the data of the elements of index (et)
% % %       %
% % %       %buffer = zeros (10 , length(et));
% % %       buffer = zeros (9 , length(et));
% % %       buffer(1,:) = et;
% % %       buffer(2,:) = type(k);
% % %       buffer(3,:) = 3;
% % %       buffer(4,:) = M.Elements.region(et);
% % %       buffer(5,:) = M.Elements.region(et);
% % %       %buffer(6,:) = 0;
% % %       for n = 1 : 4
% % %          %buffer(6+n,:) = M.Elements.tableOfNodes(et,n);
% % %          buffer(5+n,:) = M.Elements.tableOfNodes(et,n);
% % %       end
% % %       
% % %       %
% % %       % Print all the homogeneous elements in the group with a single instruction
% % %       %
% % %       fprintf (fid, elementFormat, buffer);
% % %    end
% % % end
fprintf (fid, '$EndElements\n');


fclose(fid);
