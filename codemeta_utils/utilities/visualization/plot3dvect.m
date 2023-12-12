%Plot a 3-D vector
%function plot3dvect(o,v,name,style)
%Inputs
%   o      origin of axes (3x1)
%   v      the vector (3x1)
%   name   (optional) show this string at the endpoint of the vector
%   style  (optional) specifies line style (see help plot)

%%AUTORIGHTS%%

function plot3dvect(o,v,name,style)
vect=[o o+v];
if(nargin<4)
    plot3(vect(1,:),vect(2,:),vect(3,:))
else
    plot3(vect(1,:),vect(2,:),vect(3,:), style)
end
if(nargin>2)
    text(vect(1,2),vect(2,2),vect(3,2), [' ' name]);
end
