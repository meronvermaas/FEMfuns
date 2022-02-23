function point = point_in_surf(surface)

% POINT_IN_SURF returns a point that is directly below the surface,
% by taking a face and using the normal the find this coordinate.

%select points of a triangle with fairly similarly sized edges
for ii = 1:size(surface.tri,1)
    p1 = surface.pos(surface.tri(ii,1),1:3);
    p2 = surface.pos(surface.tri(ii,2),1:3);
    p3 = surface.pos(surface.tri(ii,3),1:3);
    if std([norm(p2-p1) norm(p3-p1) norm(p3-p2)]) < 1
        break
    end
end

%TODO raise error if no triangle was selected with similar sized edges

% calculate the normal of the triangle
N = cross(p2-p1, p3-p1,2);
N = N ./ sqrt(sum(N.^2,2));

%centroid of the triangle
ctd = (p1+p2+p3)/3;

%get a point just below the surface
depth = norm(p3-p1)*1e-6; % fraction of the edge size under the triangle surface
point = ctd + N*depth;
