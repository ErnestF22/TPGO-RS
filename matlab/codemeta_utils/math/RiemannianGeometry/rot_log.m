%function A=rot_log(R1,R2)
%The inverse of R2=rot_exp(A)
%
%   See also rot_exp
%

%%AUTORIGHTS%%

function A=rot_log(R1,R2,varargin)
method='asym';

%these flags control tricks to improve the accuracy
flagRotMult=true;           %use numerically stable but slow (3x slower) rotation multiplication
flagRotSkewMult=false;      %like flagRotMult, but for multiplication of rotation and skew sym matrix
flagEnforceSkew=true;       %enforce tangent vector to be skew symmetric
flagAverageSVDBasis=true;   %average U and V in SVD decomposition of the tangent vector

%optional parameters
ivarargin=1;
while(ivarargin<=length(varargin))
    switch(lower(varargin{ivarargin}))
        case 'method'
            ivarargin=ivarargin+1;
            method=varargin{ivarargin};
        case 'flagrotmult'
            ivarargin=ivarargin+1;
            flagRotMult=varargin{ivarargin};
        otherwise
            error(['Argument ' varargin{ivarargin} ' not valid!'])
    end
    ivarargin=ivarargin+1;
end

N1=size(R1,3);
N2=size(R2,3);
if N1>1 || N2>1
    if N1==1
        R1=repmat(R1,[1 1 N2]);
        N1=N2;
    end
    if N2==1
        R2=repmat(R2,[1 1 N1]);
        N2=N1;
    end
    %recursive call if R contains multiple rotations
    if N1==N2
        A=zeros([size(R1,1) size(R1,2) N1]);
        for iN2=1:N2
            A(:,:,iN2)=rot_log(R1(:,:,iN2),R2(:,:,iN2));
        end
    else
        error('R1 should either contain one rotation or the same number of rotations as R2')
    end
else
    if(~isempty(R1))
        if flagRotMult
            R1t=R1';                %workaround for a known bug in Matlab R2009b
            R=(R1\R2+R1t/R2')/2;
        else
            R=R1'*R2;
        end
    else
        R=R2;
    end

    n=size(R,1);

    switch method
        case 'sym'
            %transform R so that it is the direct product of 2x2 rotations
%             [U,S,V]=svd((R+R')/2);
%             R=V'*R*V;
% 
%             if(mod(n,2)==0)
%                 is=1;
%             else
%                 is=2;
%             end
% 
%             %find angles of the 2x2 rotations and start making A as a direct product of
%             %2x2 skew-symmetric matrices
%             A=zeros(n);
%             while is<n
%                 stheta=R(is+1,is)-R(is,is+1);
%                 ctheta=R(is,is)+R(is+1,is+1);
% 
%                 theta=atan2(stheta,ctheta);
% 
%                 A(is,is+1)=-theta;
%                 A(is+1,is)=theta;
% 
%                 is=is+2;
%             end
% 
%             %undo the similarity transformation
%             A=V*A*V';

            [V,S]=eig((R+R')/2);
            if flagRotMult
                Rt=R';
                Vt=V';
                RV=(Rt\V+R/V')/2;
                R=(V\RV+Vt/RV')/2;
            else
                R=V'*R*V;
            end

            %find angles of the 2x2 rotations and start making A as a direct product of
            %2x2 skew-symmetric matrices
            A=zeros(n);
            for is=1:2:n-1
                stheta=R(is+1,is)-R(is,is+1);
                ctheta=R(is,is)+R(is+1,is+1);

                theta=atan2(stheta,ctheta);

                A(is,is+1)=-theta;
                A(is+1,is)=theta;
            end

            %undo the similarity transformation
            if flagRotSkewMult
                Vt=V';
                A=Vt\A/V;
            else
                A=V*A*V';
            end
                
        case 'asym'
            [U,S,V]=svd((R-R')/2);
            
            if flagAverageSVDBasis
                for is=1:2:n-1
                    s=sign(V(:,is)'*U(:,is+1));
                    U(:,is)=(U(:,is)-s*V(:,is+1))/2;
                    U(:,is+1)=(U(:,is+1)+s*V(:,is))/2;
                end
                if mod(n,2)==1
                    s=sign(V(:,end)'*U(:,end));
                    U(:,end)=(U(:,end)+s*V(:,end))/2;
                end
            end
            
            if flagRotMult
                Rt=R';
                Ut=U';
                RV=(Rt\U+R/U')/2;
                R=(U\RV+Ut/RV')/2;
            else
                R=U'*R*U;
            end

            is=1;
            
            %find angles of the 2x2 rotations and start making A as a direct product of
            %2x2 skew-symmetric matrices
            A=zeros(n);
            while is<n
                stheta=R(is+1,is)-R(is,is+1);
                ctheta=R(is,is)+R(is+1,is+1);

                theta=atan2(stheta,ctheta);

                A(is,is+1)=-theta;
                A(is+1,is)=theta;

                is=is+2;
            end

            %undo the similarity transformation
            if flagRotSkewMult
                Ut=U';
                A=Ut\A/U;
            else
                A=U*A*U';
            end
    end
    if flagEnforceSkew
        %make sure A is skew
        A=(A-A')/2;
    end

    if(~isempty(R1))
        %move to the tangent space at R1
        A=R1*A;
    end
end
