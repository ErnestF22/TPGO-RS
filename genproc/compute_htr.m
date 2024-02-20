function htr = compute_htr(x,u,Tijs,edges)
    htr = zeros(size(x.T));
    N = size(x.R,3);
    num_edges = size(edges, 1);
    for e = 1:num_edges
        ii = edges(e,1);
        jj = edges(e,2); 
        v_ri = u.R(:,:,ii);
%         v_rj = u.R(:,:,jj);
        %
        BIJ = zeros(N,1);
        BIJ(ii) = 1;
        BIJ(jj) = -1;
        %
        T_ij = Tijs(:, e);
        w_ij = 2 * BIJ * T_ij' * v_ri';
        htr = htr + w_ij';
    end
end
