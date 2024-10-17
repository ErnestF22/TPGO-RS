function hrt = compute_hrt(x,u,Tijs,edges)
    W = zeros(size(x.R));
    num_edges = size(edges, 1);
    for e = 1:num_edges
        ii = edges(e,1);
        jj = edges(e,2); 
        v_ti = u.T(:,ii);
        v_tj = u.T(:,jj);
        T_ij = Tijs(:, e);
        w_ij = 2 * (v_ti - v_tj)*T_ij';
        W(:,:,ii) = W(:,:,ii) + w_ij;
    end
    hrt = stiefel_tangentProj(x.R, W);
end