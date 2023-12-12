% Bound hat3(w)^2 matrices

close all; clear all; clc;

for i = 1:1e5
    w = randn(3,1); 
    m2 = randn;     
    m3 = rand*5;
    k = rand*10; 

    % Use the negative eigenvalue of [2*m2 m3;m3 0] to bound matrix since
    % eigenvalues of hat3(w)^2 are (0,-1)
    % Note: eig([2*m2 m3;m3 0](x)hat3(w)^2) = eig([2*m2 m3;m3 0])(x)eig(hat3(w)^2))
    eVal = eig(k*norm(w)^2*kron([2*m2 m3;m3 0],hat3(w/norm(w))^2) + ...
        k*norm(w)^2*(m2-abs(m2)-abs(m3))*eye(6));
    if ( m2-abs(m2)-abs(m3) > 1e-5)
        fprintf('scalar greater than 0\n');
    end
    if any(eVal > 1e-5)
        error("Pos eigenvalues");
    end
end

fprintf('Done!\n');