function A=rot_eucl2RiemGrad(R,A)
for(iN=1:size(A,3))
    A(:,:,iN)=A(:,:,iN)-R*A(:,:,iN)'*R;
end

