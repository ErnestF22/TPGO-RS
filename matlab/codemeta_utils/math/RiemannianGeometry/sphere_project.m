function x=sphere_project(x,radius)
if ~exist('radius','var') || isempty(radius)
    radius=1;
end
x=radius*cnormalize(x/radius);