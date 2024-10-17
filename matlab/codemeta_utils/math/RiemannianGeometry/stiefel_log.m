function [H]=stiefel_log(Y1,Y2,varargin)
%recursive call if Y2 contains multiple points
N2=size(Y2,3);
if(N2>1)
    H=zeros([size(Y1) N2]);
    for(iN2=1:N2)
        H(:,:,iN2)=stiefel_log(Y1,Y2(:,:,iN2),varargin{:});
    end
else

    method='adapt';
    %method='SO(n)';

    [~,p]=size(Y1);

    %optional parameters
    ivarargin=1;
    while(ivarargin<=length(varargin))
        switch(lower(varargin{ivarargin}))
            case 'method'
                ivarargin=ivarargin+1;
                method=varargin{ivarargin};
            otherwise
                disp(varargin{ivarargin})
                error('Argument not valid!')
        end
        ivarargin=ivarargin+1;
    end

    if strcmp(method,'adapt')
        if(norm(Y1-Y2,'fro')<1e-8)
            method='so(n)';
        else
            method='so(2p)';
        end
    end
    
    switch lower(method)
        %use geodesics in SO(n)
        case 'so(n)'
            %make Y1 and Y2 square matrices in SO(n)
            YY1=orthCompleteBasis(Y1);
            YY2=orthCompleteBasis(Y2);
            %move YY2 inside the equivalence class in order to minimize
            %distance to YY1
            [YY2]=stiefel_log_findOrthGrad(YY1,YY2,p);
            %compute the log map
            H=YY1*rot_log([],YY1'*YY2);   %projection of H on the vertical space should be zero here

        %use geodesics in SO(2p)
        %Basically the inverse of the algorithm for stiefel_exp
        case {'so(2p)','so(2*p)'}
            %make n by 2p orthogonal matrix with the same span as Y1 and Y2 and
            %first columns equal to Y1
            %Y1Q=[Y1 orth((eye(n)-Y1*Y1')*Y2)];
            Y1Q=orthCompleteBasis(Y1, Y2);
            %find what linear combination of the columns of Y1Q generates Y2
            MN=reshape(kron(eye(p),Y1Q)\Y2(:),size(Y1Q,2),size(Y2,2));
            %complete this to be a matrix in SO(2p)
            MN=orthCompleteBasis(MN);
            %now we solve: find T in SO(2p) s.t. projVertSpace(logm(T*MNbigest))=0
            T=stiefel_log_findOrthGrad(MN',eye(size(MN)),size(Y1,2));
            %compute log map from MN
            H=Y1Q*rot_log([],MN*T);
    end
    H=H(:,1:p);
end
