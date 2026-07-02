function lambdas_out = ssom_manopt_step3(tijs, edges, R, T)


num_edges = size(edges, 1);
lambdas_out = zeros(num_edges, 1);
for ee = 1:num_edges
    
    tij = tijs(:, ee);

    ii = edges(ee, 1);
    jj = edges(ee, 2);
    Ri = R(:,:,ii);
    % Rj = R(:,:,jj);
    Ti = T(:,ii);
    Tj = T(:,jj);

    a = tij' * tij;
    b1 = tij' * Ri' * Ti;
    b2 = -tij' * Ri' * Tj;
    b3 = Ti' * Ri * tij;
    b4 = -Tj' * Ri * tij;
    b = b1 + b2 + b3 + b4;
    c1 = Ti' * Ti;
    c2 = -Ti' * Tj;
    c3 = -Tj' * Ti; % OBS: should be equal to c2
    c4 = Tj' * Tj;
    c = c1 + c2 + c3 + c4;

    x = solve_quadratic_lambdas(a,b,c);

    lambdas_out(ee) = x;
end


end %file function
