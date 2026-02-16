function is_tg_bool = check_is_tangent_stiefel(X, Y, thr)

if ~exist('thr','var')
  thr = 1e-6;
end

x_is_on_stiefel_bool = check_is_on_stiefel(X, thr); %first, check whether X is on Stiefel

if ~x_is_on_stiefel_bool
    error("X is NOT on Stiefel")
end



if (size(X) ~= size(Y))
    error("Bad sizes")
end

N = size(X, 3);

is_tg_bool = true;
for ii = 1:N
    X_ii = X(:,:,ii);
    Y_ii = Y(:,:,ii);
    if ~is_equal_floats(X_ii' * Y_ii, -Y_ii' * X_ii)
        is_tg_bool = false;
        break;
    end
end

end %file function
