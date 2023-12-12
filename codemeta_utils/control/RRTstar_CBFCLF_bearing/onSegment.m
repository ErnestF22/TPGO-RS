%Given three colinear points p, q, r, the function checks if  
%point q lies on line segment 'pr' 
function flag = onSegment(p, q, r)
flag = false;
    if (q(1) <= max(p(1), r(1))) && (q(1) >= min(p(1), r(1)))
          if (q(2) <= max(p(2), r(2))) && (q(2) >= min(p(2), r(2)))
              flag = true;
          end
    end
end