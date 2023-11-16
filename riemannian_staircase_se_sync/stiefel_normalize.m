function v_out = stiefel_normalize(v_in)
%STIEFEL_NORMALIZE Summary of this function goes here
%   Detailed explanation goes here
if size(v_in, 3) < 2
    v_out = v_in ./ stiefel_metric([], v_in, v_in, "euclidean");
else
    v_out = zeros(size(v_in));
    for ii = 1:size(v_in, 3)
        v_in_ii = v_in(:,:,ii);
        v_out(:,:,ii) = v_in_ii ./ stiefel_metric([], v_in_ii, v_in_ii, "euclidean");
    end
end


end

