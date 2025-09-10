function geod = stiefel_geodFun_eucl(Y, H, A, S0)
%generate geodesic emanate from Y in direction H
%from corollary 2.2 in Edelman-Arias-Smith (page 12)

if size(Y) ~= size(A)
    error("sizes don't match")
end

% nrs = size(Y, 1);
d = size(Y, 2);
N = size(Y, 3);

if N>1
    error("Currently implemented only for single-Stiefel")
end


for ii = 1:N
    Y_ii = Y(:,:,ii);
    H_ii = H(:,:,ii);
    A_ii = A(:,:,ii);
    Y0 = Y_ii;
    Ydot0 = H_ii;
    % S_t = @(t) expm(t * A_ii) * S0 * expm(-t * A_ii);
    geod = @(t) [Y0, Ydot0] * expm(t * [A_ii, -S0; eye(size(A_ii)), A]) * ...
        eye(2*d, d) * expm(-t * A_ii);
end


end %file function
