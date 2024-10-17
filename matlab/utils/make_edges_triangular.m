function edges = make_edges_triangular(N)
%MAKE_EDGES_TRIANGULAR Returns all correspondences between 2 vectors of
%N elements, i.e. a correspondences 2x(N*(N-1)/2) elements as follows:
%correspondences = [1 2; 1 3; ...; 1 N; 2 3; ...; 2 N; ...; N-1 N]'
    edges = zeros(2, 0.5*N*(N-1));
    full_id = 1;
    for ii = 1:N-1
        for jj = ii+1:N
            new_corr = [ii jj]';
            edges(:, full_id) = new_corr;
            full_id = full_id + 1;
        end
    end
end

