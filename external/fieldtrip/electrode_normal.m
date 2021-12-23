function elnormal = electrode_normal(surface, electrode_coordinates, surfsmooth)


% ELECTRODE_NORMAL returns the electrode normal based on the orientation
% of the nearest triangle of a triangulated surface.
% Smooth the surface with surfsmooth to ensure local curvature in the surface does not affect the electrode orientation

if ~exist('surfsmooth', 'var')
    surfsmooth = 0.5;
end

% Check that all the elements are correctly oriented
surface.tri = meshreorient(surface.pos,surface.tri);

% Keep $surfsmooth$% (default 50%) of the elements after the sampling to ensure sufficiently smooth surface
if surfsmooth ~=1 %when 1, don't resample
    [surface.pos,surface.tri] = meshresample(surface.pos,surface.tri,surfsmooth);
end

% calculate the normals of each triangle taking the cross product of two sides of it
surface.N = surfacenorm(surface.pos,surface.tri);

% find orientation of electrodes taking the normal of the nearest triangle of the surface
surface.ctd = meshcentroid(surface.pos,surface.tri);
[D,I] = pdist2(surface.ctd, electrode_coordinates, 'euclidean', 'Smallest',1);

elnormal = surface.N(I,:);
