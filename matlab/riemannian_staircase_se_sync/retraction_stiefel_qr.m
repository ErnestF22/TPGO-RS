function rxe = retraction_stiefel_qr(x,e)
%RETRACTION_STIEFEL Compute and outputs rxe = retraction on Stiefel manifold
% x is on the Stiefel manifolds and e is on the tangent space to the
% manifold at x

%TODO: add check size(x) == size(e)

N = size(x,3);
p = size(x,2);

if N > 1
    rxe = zeros(size(x));
    for ii = 1:N
        x_ii = x(:,:,ii);
        e_ii = e(:,:,ii);
        [Q,R] = qr(x_ii + e_ii);
        rxe(:,:,ii) = R;
    end
else
%     rxe = zeros(size(x));
    [Q,R] = qr(x + e);
    rxe = R;
end



end %file function
