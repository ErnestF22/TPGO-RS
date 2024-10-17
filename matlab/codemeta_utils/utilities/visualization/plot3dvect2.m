%PLOT3DVECT2(O,V,NAME,STYLE)
%
%O      first point of the vector (3x1)
%V      second point of the vector (3x1)
%NAME   (optional) show this string at the endpoint of the vector
%STYLE  (optional) specifies line style (see help plot)
function plot3dvect2(o,v,name,style)
vect=[o v];
if(nargin<4)
    plot3(vect(1,:),vect(2,:),vect(3,:))
else
    plot3(vect(1,:),vect(2,:),vect(3,:), style)
end
if(nargin>2)
    text(vect(1,2),vect(2,2),vect(3,2), [' ' name]);
end
