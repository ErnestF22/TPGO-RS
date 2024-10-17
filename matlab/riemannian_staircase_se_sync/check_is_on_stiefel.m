function is_on_stiefel_bool = check_is_on_stiefel(x,thr)
%CHECK_IS_ON_STIEFEL Check if x is on (3D) Stiefel manifold
% i.e., x'(:,:,ii) * x(:,:,ii) = eye(size(x,2)) for all ii = 1:size(x,3)

if ~exist('thr','var')
  thr = 1e-6;
end

val = max(abs( multiprod(multitransp(x), x) - ...
    eye3d(size(x,2), size(x,2), size(x,3)) ), [], 'all');
if val < thr
    is_on_stiefel_bool = boolean(1);
else
    is_on_stiefel_bool = boolean(0);
end

end %file function

