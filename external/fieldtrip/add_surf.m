function newsurf = add_surf(surface1, surface2, marker1, marker2)


% ADD_SURF combine two surfaces and remove duplicate nodes and faces
% surface1 is the base, new nodes and faces of surface2 are added to it
% NOTE: self-intersections or non-connected surfaces are not resolved

% boundary labeling
if size(surface1.tri,2) < 4 && size(surface2.tri,2) < 4
    if (~exist('marker1', 'var') && ~exist('marker2', 'var'))
        marker1 = 1;
        marker2 = 2;
    elseif (~exist('marker1', 'var') && exist('marker2', 'var'))
        marker1 = marker2+1;
    elseif (exist('marker1', 'var') && ~exist('marker2', 'var'))
        marker2 = marker1+1;
    end
    surface1.tri(:,4) = marker1;
    surface2.tri(:,4) = marker2;
elseif size(surface1.tri,2) == 4 && size(surface2.tri,2) < 4
    if (~exist('marker2', 'var')) , marker2 = max(surface1.tri(:,4))+1; end
    surface2.tri(:,4) = marker2;
elseif size(surface1.tri,2) < 4 && size(surface2.tri,2) == 4
    if (~exist('marker1', 'var')) , marker1 = max(surface2.tri(:,4))+1; end
    surface1.tri(:,4) = marker1;
elseif size(surface1.tri,2) == 4 && size(surface2.tri,2) == 4
    surf1lab = unique(surface1.tri(:,4));
    surf2lab = unique(surface2.tri(:,4));
    labidx = ismember(surf1lab,surf2lab);
    for ii = surf1lab(labidx)
        newlab = max(unique([surface1.tri(:,4) surface2.tri(:,4)]))+1;
        surface1.tri(surface1.tri(:,4)==ii,4) = newlab;
    end
end

%add new nodes
flag = ~ismember(surface2.pos,surface1.pos,'rows');
newsurf.pos = [surface1.pos; surface2.pos(flag,:)];
newsurf.tri = surface1.tri;

%loop over faces of the to be added surface
for i = 1:size(surface2.tri,1)
    %find the new face indices
    tmpfc = [];
    for j = 1:3
        [flag,idx] = ismember(surface2.pos(surface2.tri(i,j),:), newsurf.pos, 'rows');
        if length(idx) == 1
            tmpfc = [tmpfc idx];
        else
            error('Only one vertex coordinate index should be found. Remove duplicates first')
        end
    end
    
    tmpfc = [tmpfc surface2.tri(i,4)];

    %check if it is a new face, then add
    flag = ~ismember(sort(tmpfc(:,1:3)),sort(newsurf.tri(:,1:3),2),'rows');
    if sum(flag)==1
        newsurf.tri = [newsurf.tri;tmpfc];
    end
end
