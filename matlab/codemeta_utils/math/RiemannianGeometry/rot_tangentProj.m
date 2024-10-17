%Project a matrix A on the tangent space of SO(n) at R
%If A contains multiple vectors (3-D matrix) but R is single (2-D matrix),
%project all the vectors in the tangent space at R
%If also R contains multiple rotations, project A(:,:,iN) on the tangent
%space of R(:,:,iN)
function A=rot_tangentProj(R,A)
NR=size(R,3);
NA=size(A,3);

if NR>1 && NR~=NA
    error('If R contains multiple rotations, A must contain an equal number of matrices')
end

if NR==1
    for iN=1:NA
        A(:,:,iN)=0.5*(A(:,:,iN)-R*A(:,:,iN)'*R);
    end
else
    for iN=1:NA
        A(:,:,iN)=rot_tangentProj(R(:,:,iN),A(:,:,iN));
    end
end

        