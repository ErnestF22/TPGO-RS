function bool_is_rot = check_is_rotation(A)
% CHECK_IS_ROTATION Returns a boolean stating whether matrix A is on SO(N)
% manifold i.e., if its determinant is 1 and A * A' = A' * A = eye(N)

if (det(A) < 1-1e-6 || det(A) > 1+1e-6)
    disp("check_is_rotation() -> determinant = 1 check failed")
    fprintf("det = %g\n", det(A))
    bool_is_rot = boolean(0);
    return;
end
if max(abs(A * A' - A' * A)) > 1e-6
    disp("check_is_rotation() -> A*A' == A'*A check failed")
    bool_is_rot = boolean(0);
    return;
end
if max(abs(A * A' - eye(size(A)))) > 1e-6
    disp("check_is_rotation() -> A*A' == I check failed")
    bool_is_rot = boolean(0);
    return;
end
bool_is_rot = boolean(1);

end %file function