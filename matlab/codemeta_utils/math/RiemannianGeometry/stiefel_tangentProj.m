%Projects H onto the tangent space of the Stiefel manifold at Y
function H_out=stiefel_tangentProj(Y,H)
sym=@(A) (A+A')/2;

N = size(H,3);
if N > 1
    H_out = zeros(size(H));
    for ii = 1:N
        Hii = H(:,:,ii);
        Yii = Y(:,:,ii);
        H_out(:,:,ii) = Hii-Yii*sym(Yii'*Hii);
    end
else
    H_out=H-Y*sym(Y'*H);
end

end %function
