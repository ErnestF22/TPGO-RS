function x = solve_quadratic_lambdas(a,b,c)
% Returns 2-by-N array of roots for coefficients a,b,c (scalars or arrays)
% Usage: r = solveQuadratic(1, -3, 2)   % -> [2; 1]
D = b.^2 - 4.*a.*c;                    % discriminant

if D < 0
    x = 1;
    return;
else
    sqrtD = sqrt(D);
    x1 = (-b + sqrtD) / (2.*a);
    x2 = (-b - sqrtD) / (2.*a);
    if x1 < 1 && x2 >=1
        x = x2;
    elseif  x2 < 1 && x1 >=1
        x = x1;
    elseif x1 >= 1 && x2 >= 1
        x = max(x1, x2);
    else
        %which case are we in?
        x = 1;
    end
end

end %file function