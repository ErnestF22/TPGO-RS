%Compute the Lie bracket of two tangent vectors in SO(3)
%function C=rot_bracket(R,A,B)
function C=rot_bracket(R,A,B)
%If R,A,B are vector fields, return a function handle with parameter R
if ( isa(R, 'function_handle') )
    C = @(R) rot_bracket(R,A(R),B(R));
else
    % Compute Lie bracket with the given numerical vectors
    NA=size(A,3);
    NB=size(B,3);
    C=zeros([size(A,1) size(A,2) NA NB]);
    for iNA=1:NA
        for iNB=1:NB
            Ap=R'*A(:,:,iNA);
            Bp=R'*B(:,:,iNB);
            C(:,:,iNA,iNB)=R*(Ap*Bp-Bp*Ap);
        end
    end
end

