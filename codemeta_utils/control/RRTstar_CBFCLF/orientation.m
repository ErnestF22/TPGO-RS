% to find the orientation of an ordered triplet (p,q,r)
% function returns the following values:
% 0 : Colinear points
% 1 : Clockwise points
% 2 : Counterclockwise
function flag = orientation(p, q, r)
    val = ((q(2) - p(2)) * (r(1) - q(1))) - ((q(1)- p(1)) * (r(2) - q(2))); 
    if (val > 0)
%         Clockwise orientation 
        flag = 1;
    elseif (val < 0)
%          Counterclockwise orientation 
        flag = 2;
    else
%       Colinear orientation 
        flag = 0;
    end
end