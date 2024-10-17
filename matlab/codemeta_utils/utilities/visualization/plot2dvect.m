%Plot a 2-D vector
%function plot2dvect(o,v,name,style)
%Inputs
%   o      origin of axes (3x1)
%   v      the vector (3x1)
%   name   (optional) show this string at the endpoint of the vector
%   style  (optional) specifies line style (see help plot)

%%AUTORIGHTS%%

function plot2dvect(o,v,name,style)
vect=[o o+v];
if(nargin<4)
    plot(vect(1,:),vect(2,:))
else
    plot(vect(1,:),vect(2,:), style)
end
if(nargin>2)
    text(vect(1,2),vect(2,2), [' ' name]);
end
