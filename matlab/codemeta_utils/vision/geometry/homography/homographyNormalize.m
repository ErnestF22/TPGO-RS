%Normalize scale (and possibly sign) of an homography
%function H=homographyNormalize(H,x1,x2)
%Normalizes second singular value of H to plus or minus one and, if x1,x2
%are provided, chooses also the correct sign using depts
function H=homographyNormalize(H,x1,x2)
flagNormalizeSign=false;
if exist('x1','var')
    flagNormalizeSign=true;
    if ~exist('x2','var') || isempty(x2)
        [x1,x2]=homographySplitData(x1);
    end
end

N=size(H,3);
if flagNormalizeSign && size(x2,3)~=N
    error('size(x2,3) and size(H,3) must be the same.')
end
if N>1
    for iN=1:N
        if flagNormalizeSign
            H(:,:,iN)=homographyNormalize(H(:,:,iN),x1,x2(:,:,iN));
        else
            H(:,:,iN)=homographyNormalize(H(:,:,iN));
        end
    end
else
    s=svd(H);
    H=H/s(2);

    if flagNormalizeSign
        lambdap=homographyTriangulateDepths(H,x1,x2);
        lambdan=homographyTriangulateDepths(-H,x1,x2);

        %choose solution with the most coherent signs
        slambdap=sum(prod(sign(lambdap)));
        slambdan=sum(prod(sign(lambdan)));
        if slambdan<slambdap
            H=-H;
        end
    end
end
