% The main function that returns true if  
% the line segment 'p1-q1' and 'p2-q2' intersect. 
function flag = doIntersect(p1,q1,p2,q2)
%Find the 4 orientations required for
% the general and special cases
o1 = orientation(p1, q1, p2);
o2 = orientation(p1, q1, q2);
o3 = orientation(p2, q2, p1);
o4 = orientation(p2, q2, q1);
  
% General case 
if ((o1 ~= o2) && (o3 ~= o4))
        flag =  true;
  
%Special Cases
%p1 , q1 and p2 are colinear and p2 lies on segment p1q1 
elseif ((o1 == 0) && onSegment(p1, p2, q1))
        flag = true;
  
%p1 , q1 and q2 are colinear and q2 lies on segment p1q1 
elseif ((o2 == 0) && onSegment(p1, q2, q1))
        flag = true;
  
% p2 , q2 and p1 are colinear and p1 lies on segment p2q2 
elseif ((o3 == 0) && onSegment(p2, p1, q2))
        flag = true;
  
% p2 , q2 and q1 are colinear and q1 lies on segment p2q2 
elseif ((o4 == 0) && onSegment(p2, q1, q2))
        flag = true;
  
else
    flag = false;
end
  
end