function edges_full = make_edges_full(N)
%MAKE_CORRESPONDENCES_FULL Returns all correspondences between 2 vectors of
%N elements, i.e. a correspondences 2x(N^2-N) elements as follows:
%correspondences = [1 2; 1 3; ...; 1 N; 2 1; 2 3 2 4; ...; 2 N; ...; N 1; N 2; ...; N-1 N]'
    edges_full = zeros(N*N-N, 2);
    full_id = 1;
    for ii = 1:N
        for jj = 1:N
            if ii == jj
                continue;
            end
            new_corr = [ii jj]';
            edges_full(full_id, :) = new_corr;
            full_id = full_id + 1;
        end
    end
end

