function tijs_mat = tijs_vec_2_tijs_mat(tijs_vec, edges, N)
%TIJS_VEC_2_TIJS_MAT: Deloy position-wise the tijs_vec that comes from codemeta's
%testNetwork .gij member

d = size(tijs_vec, 1);

tijs_mat = zeros(d*N,N);
for ee = 1:size(edges, 1)
    ii = edges(ee, 1);
    jj = edges(ee, 2);
    tijs_mat(ii*d-(d-1):ii*d, jj) = tijs_vec(:,ee);
end

end %function