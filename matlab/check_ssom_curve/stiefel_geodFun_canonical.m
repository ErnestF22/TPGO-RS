function geod = stiefel_geodFun_canonical(Y, H)
%generate geodesic emanate from Y in direction H
%from corollary 2.2 in Edelman-Arias-Smith (page 12)

if size(Y) ~= size(H)
    error("sizes don't match")
end

nrs = size(Y, 1);
d = size(Y, 2);
N = size(Y, 3);

if N>1
    error("Currently implemented only for single-Stiefel")
end


for ii = 1:N
    Y_ii = Y(:,:,ii);
    H_ii = H(:,:,ii);
    K_ii = (eye(nrs) - Y_ii * Y_ii') * H_ii;
    A_ii = Y_ii' * H_ii;
    % disp("max(abs(A_ii+A_ii'), [], ""all"")");
    % disp(max(abs(A_ii+A_ii'), [], "all"));
    [Q_ii, R_ii] = qr(K_ii, "econ");
    tmp = [A_ii, -R_ii'; R_ii, zeros(size(R_ii))];
    tmp_t = @(t) expm(t * tmp) * [eye(d); zeros(size(R_ii))];
    get_first_d_rows = @(A, d) A(1:d, :);
    get_second_d_rows = @(A, d) A(d+1:end, :);
    M = @(t) get_first_d_rows(tmp_t(t), d);
    N = @(t) get_second_d_rows(tmp_t(t), d);
    geod = @(t) Y_ii * M(t) + Q_ii * N(t);
end


end %file function
