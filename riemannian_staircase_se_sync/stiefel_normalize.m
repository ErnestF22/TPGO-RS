function v_out = stiefel_normalize(x_in, v_in)
%STIEFEL_NORMALIZE Normalize tangent vector v_in (tg. to x_in)
%

if size(v_in, 3) < 2
    v_out = v_in ./ stiefel_metric(x_in, v_in, v_in, "canonical");
else
    v_out = zeros(size(v_in));
    for ii = 1:size(v_in, 3)
        v_in_ii = v_in(:,:,ii);
        v_out(:,:,ii) = v_in_ii ./ ...
            sum(stiefel_metric(x_in, v_in_ii, v_in_ii, "canonical"));
    end
end


end

