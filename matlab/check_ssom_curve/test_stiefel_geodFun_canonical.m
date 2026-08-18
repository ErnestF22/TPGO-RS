function test_stiefel_geodFun_canonical

nrs = 4; 
d = 3; 
N = 1;
Y = make_rand_stiefel_3d_array(nrs, d, N);

M = stiefelfactory(nrs, d, N);
H = M.randvec(Y);


geod_canonical = stiefel_geodFun_canonical(Y, H);

ts = linspace(-1, 1, 100);

for t = ts
    geod_canonical_t = geod_canonical(t);
    disp("check_is_on_stiefel(geod_canonical_t)");
    disp(check_is_on_stiefel(geod_canonical_t));
end

end %file function
