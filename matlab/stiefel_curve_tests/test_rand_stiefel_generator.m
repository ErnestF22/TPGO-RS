iter_max = 1000;

% array_type = "double";

eps = 1e-6;

nrs = 4;
d = 3;
N = 5;

broken = boolean(0);
for ii = 1:iter_max
    M_stiefel = make_rand_stiefel_3d_array(nrs, d, N);
%     disp(M_stiefel)
    if max( ...
        eye3d(d,d,N) - multiprod(multitransp(M_stiefel), M_stiefel), ...
            [], "all") > eps
        disp("broken");
        disp(M_stiefel);
        broken = boolean(1);
        break;
    end
    if broken
        break;
    end
end


