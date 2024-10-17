function retval = eye3d(a,b,c)
%EYE3D Return eye matrix with size axb in each element of the c dimension
retval = repmat(eye(a,b), 1, 1, c);
end

