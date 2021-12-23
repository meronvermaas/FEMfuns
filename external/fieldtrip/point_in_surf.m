function point = point_in_surf(surface, depth)

% POINT_IN_SURF returns a point that is directly below the surface.
% the depth of it under the surface is default set to 1e-6 to make sure the
% point is not in another surface.

if ~exist('depth', 'var') || isempty(depth)
    depth = 1e-6;
end

% center of the surface
surfcenter = mean(surface.pos,1);

% calculate the normal of the first triangle
p1 = surface.pos(surface.tri(1,1),:);
p2 = surface.pos(surface.tri(1,2),:);
p3 = surface.pos(surface.tri(1,3),:);
N = cross(p2-p1, p3-p1,2);
N = N ./ sqrt(sum(N.^2,2));

%centroid of the triangle
ctd = (p1+p2+p3)/3;

%make sure the normal points outside wrt the center of the surface
sdot = sign(dot(ctd-surfcenter,N,3));
N = sdot .* N;

%get a point just below the surface
point = ctd - N*depth;
