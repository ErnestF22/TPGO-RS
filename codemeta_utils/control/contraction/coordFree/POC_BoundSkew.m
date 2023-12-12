% Bound hat3(w) matrices

clear all; close all; clc;

for i = 1:1e5
    A = rand(2,2);
    B = (A+A')/2+2*eye(2);
    m = [B(1,1);B(1,2);B(2,2)];
    kv = 0.1;
    w = randn(3,1);
    theta = (m(2)/4+m(3)*kv/4)*norm(w);

    AA = theta*[zeros(3) -hat3(w/norm(w));hat3(w/norm(w)) zeros(3)];
    BB = theta*[eye(3) zeros(3);zeros(3) eye(3)];

    eVal = eig(AA-BB);
    if (any(eVal > 1e-5))
        error('AA-BB is not <= 0');
    end
    
    if ( any(m < 0) )
        fprintf('m < 0\n');
        m
        
    end
end

fprintf('Test done!\n')