function rxe = retraction_stiefel(x,e)
%RETRACTION_STIEFEL Compute and outputs rxe = retraction on Stiefel manifold
% x is on the Stiefel manifolds and e is on the tangent space to the
% manifold at x

%TODO: add check size(x) == size(e)

N = size(x,3);
p = size(x,2);

if N > 1
    rxe = zeros(size(x));
    i_p = eye3d(p,p,N);
    for ii = 1:N
        x_ii = x(:,:,ii);
        e_ii = e(:,:,ii);
        snd_term_ii = inv(sqrtm((i_p(ii) + e_ii' * e_ii)));
        rxe(:,:,N) = (x_ii + e_ii)*snd_term_ii;
    end
else
%     rxe = zeros(size(x));
    i_p = eye(p,p);
    snd_term = inv(sqrtm((i_p + e' * e)));
    rxe = (x + e)*snd_term;
end



end %file function
