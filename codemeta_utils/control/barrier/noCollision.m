function nc = noCollision(n2, n1, o)
    A = [n1(1) n1(2)];
    B = [n2(1) n2(2)];
    
    % Find line passing through points n1 and n2
    if ( abs(n1(1)-n2(1)) < 1e-5)
        m = Inf;
        b = (n1(1)+n2(1))/2;
    elseif ( abs(n1(2)-n2(2)) < 1e-5 )
        m = 0;
        b = (n1(2)+n2(2))/2;
    else
        coeff = polyfit([n1(1) n2(1)], [n1(2) n2(2)],1);
        m = coeff(1);
        b = coeff(2);
    end
    
    nc = 1; % Assume no collision by default
    
    % Check for collision with all obstacles, if any collision return true
    for i = 1:size(o,1)
        o_loop = o(i,:);
        
%         % CIRCULAR OBSTACLES
%         radius = o_loop(3)/2;
%         centerX = o_loop(1) + radius;
%         centerY = o_loop(2) + radius;
%         
%         %Check if A or B is inside the obstacle
%         if ( (n1(1)-centerX)^2 + (n1(2)-centerY)^2 <= radius^2 ||...
%                 (n2(1)-centerX)^2 + (n2(2)-centerY)^2 <= radius^2)
%             nc = 0;
%             return;
%         end
%         
%         % Check if the line segment between A and B intersect the circle
%         [xout, yout] = linecirc(m,b,centerX,centerY,radius);
%         for j = 1:2
%             if ( ~isnan(xout(j)) )
%                 % Check if intersect point is between point A and B
%                 if ( xout(j) >= min(n1(1),n2(1)) &&...
%                         xout(j) <= max(n1(1),n2(1)) &&...
%                         yout(j) >= min(n1(2),n2(2)) &&...
%                         yout(j) <= max(n2(2),n2(2)) )
%                     nc = 0;
%                     return;
%                 end
%             end            
%         end
%         % CIRCULAR OBSTACLES
        
        
        % RECTANGULAR OBSTACLES
        obs = [o_loop(1) o_loop(2) o_loop(1)+o_loop(3) o_loop(2)+o_loop(4)];
    
        C1 = [obs(1),obs(2)];
        D1 = [obs(1),obs(4)];
        C2 = [obs(1),obs(2)];
        D2 = [obs(3),obs(2)];
        C3 = [obs(3),obs(4)];
        D3 = [obs(3),obs(2)];
        C4 = [obs(3),obs(4)];
        D4 = [obs(1),obs(4)];

        % Check if path from n1 to n2 intersects any of the four edges of the
        % obstacle

        ints1 = ccw(A,C1,D1) ~= ccw(B,C1,D1) && ccw(A,B,C1) ~= ccw(A,B,D1); 
        ints2 = ccw(A,C2,D2) ~= ccw(B,C2,D2) && ccw(A,B,C2) ~= ccw(A,B,D2);
        ints3 = ccw(A,C3,D3) ~= ccw(B,C3,D3) && ccw(A,B,C3) ~= ccw(A,B,D3);
        ints4 = ccw(A,C4,D4) ~= ccw(B,C4,D4) && ccw(A,B,C4) ~= ccw(A,B,D4);
        if ints1==0 && ints2==0 && ints3==0 && ints4==0
            nc = 1; % No collision
        else
            nc = 0; % Collision detected
            return;
        end
        % RECTANGULAR OBSTACLES
        
    end    
end