%R=rot(w,theta) computes the 3 by 3 matrix corresponding to a rotation 
%around the unit axis w by amount theta

function R = rot(w,theta)

w=shiftdim(w);

if(size(w,2)==1)
    temp=norm(w);
%     if temp < 1e-12
%       R = eye(3,3);
%     else 
      if nargin==1
        theta=temp;
      end
      if(temp<2*eps)
          theta=0;
          w=[0;0;0];
      else
          w=w/temp;
      end
      R = eye(3,3) - hat3(w)*sin(theta) + hat3(w)^2*(1-cos(theta));
%     end
else
    if (size(w,2)==3)
        temp=norm(vee(w));
%         if temp < 1e-12
%           R = eye(3,3);
%         else 
          if nargin==1
            theta=temp;
          end
          if(temp<2*eps)
              theta=0;
              w=[0;0;0];
          else
              w=w/temp;
          end
          R = eye(3,3) + w*sin(theta) + w^2*(1-cos(theta));
%         end
    else
        error('w neither a vector nor skew-symmetric???')
    end
end