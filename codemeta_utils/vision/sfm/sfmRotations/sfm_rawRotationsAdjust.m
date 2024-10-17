%Adjust a vector of rotations with a linear transformation
%function [RiK,K]=sfm_rawRotationsAdjust(Ri,varargin)
%Finds the matrix K such that RiVec*K best approximates a stack of
%rotations
function [RiK,K]=sfm_rawRotationsAdjust(Ri,varargin)
method='spd';
flagLeft=false;
flagStacked=size(Ri,3)==1;

%optional parameters
ivarargin=1;
while ivarargin<=length(varargin)
    switch(lower(varargin{ivarargin}))
        case 'method'
            ivarargin=ivarargin+1;
            method=lower(varargin{ivarargin});
        case 'left'
            flagLeft=true;
        case 'right'
            flagLeft=false;
        otherwise
            disp(varargin{ivarargin})
            error('Argument not valid!')
    end
    ivarargin=ivarargin+1;
end

if flagStacked
    Ri=matUnstack(Ri);
end
if flagLeft
    Ri=multitransp(Ri);
end
NR=size(Ri,3);
idxA=reshape(1:6*NR,6,NR);
b=repmat([1;0;0;1;0;1],NR,1);
switch method
    case 'linear'
        A=zeros(6*NR,6);
        for iR=1:NR
            RCur=Ri(:,:,iR);
            Ai=[    RCur(1,1)^2,           2*RCur(1,1)*RCur(1,2),           2*RCur(1,1)*RCur(1,3),    RCur(1,2)^2,           2*RCur(1,2)*RCur(1,3),    RCur(1,3)^2; ...
                RCur(1,1)*RCur(2,1), RCur(1,1)*RCur(2,2) + RCur(1,2)*RCur(2,1), RCur(1,1)*RCur(2,3) + RCur(1,3)*RCur(2,1), RCur(1,2)*RCur(2,2), RCur(1,2)*RCur(2,3) + RCur(1,3)*RCur(2,2), RCur(1,3)*RCur(2,3); ...
                RCur(1,1)*RCur(3,1), RCur(1,1)*RCur(3,2) + RCur(1,2)*RCur(3,1), RCur(1,1)*RCur(3,3) + RCur(1,3)*RCur(3,1), RCur(1,2)*RCur(3,2), RCur(1,2)*RCur(3,3) + RCur(1,3)*RCur(3,2), RCur(1,3)*RCur(3,3); ...
                   RCur(2,1)^2,           2*RCur(2,1)*RCur(2,2),           2*RCur(2,1)*RCur(2,3),    RCur(2,2)^2,           2*RCur(2,2)*RCur(2,3),    RCur(2,3)^2; ...
                RCur(2,1)*RCur(3,1), RCur(2,1)*RCur(3,2) + RCur(2,2)*RCur(3,1), RCur(2,1)*RCur(3,3) + RCur(2,3)*RCur(3,1), RCur(2,2)*RCur(3,2), RCur(2,2)*RCur(3,3) + RCur(2,3)*RCur(3,2), RCur(2,3)*RCur(3,3); ...
                   RCur(3,1)^2,           2*RCur(3,1)*RCur(3,2),           2*RCur(3,1)*RCur(3,3),    RCur(3,2)^2,           2*RCur(3,2)*RCur(3,3),    RCur(3,3)^2];

            A(idxA(:,iR),:)=Ai;
        end
        g=A\b;
        G=[g(1) g(2) g(3); g(2) g(4) g(5); g(3) g(5) g(6)];
        G=matProjectPSD(G);
    case 'spd'
        A=zeros(6*NR,9);
        idxSel=[1 2 3 5 6 9];
        for iR=1:NR
            Ai=kron(Ri(:,:,iR),Ri(:,:,iR));
            A(idxA(:,iR),:)=Ai(idxSel,:);
        end
        cvx_quiet(true)
        cvx_begin
            variable G(3,3)
            minimize norm(A*G(:)-b)
            subject to
                G==semidefinite(3)
        cvx_end
end

K=chol(G)';
RiK=multiprod(Ri,K);

if flagLeft
    RiK=permute(RiK,[2 1 3]);
end
if flagStacked
    RiK=matStack(RiK);
end
