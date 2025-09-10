function test_stiefel_geodFun_euclidean

nrs = 4; 
d = 3; 
N = 1;
Y = make_rand_stiefel_3d_array(nrs, d, N);

M = stiefelfactory(nrs, d, N);
H = M.randvec(Y);
H = H / norm(H);
disp("check_is_tangent_stiefel(Y, H)")
disp(check_is_tangent_stiefel(Y, H))

A = 10 * rand(d, d, N);
A = 0.5*(A - multitransp(A));
A = A / norm(A);

S = 10 * rand(d, d, N);
S = 0.5*(S+multitransp(S));
S = S / norm(S);

geod_euclidean = stiefel_geodFun_eucl(Y, H, A, S);

ts = linspace(0, 1, 100);

for t = ts
    geod_euclidean_t = geod_euclidean(t);
    disp("t")
    disp(t)
    disp("check_is_on_stiefel(geod_euclidean_t)");
    disp(check_is_on_stiefel(geod_euclidean_t));
end

end %file function
