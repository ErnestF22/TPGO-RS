function [R,T,lambda]=epipolarEToRT(E,x1,x2,varargin)
flagProgressBar=false;

NE=size(E,3);
NP=size(x1,2);

R=zeros(3,3,NE);
T=zeros(3,NE);
lambda=zeros(2,NP,NE);

ivarargin=1;
while(ivarargin<=length(varargin))
    switch(lower(varargin{ivarargin}))
        case 'progressbar'
            flagProgressBar=true;
        otherwise
            disp(varargin{ivarargin})
            error('Argument not valid!')
    end
    ivarargin=ivarargin+1;
end

transfPairs = repmat(struct('T',[],'R',[],'lambda',[]),[4 1]);
if flagProgressBar
    w=getTextWaitBar(NE);
end
for iE=1:NE
    [U,S,V] = svd(E(:,:,iE));

    %U,V must be in SO(3)!!!
    U = det(U)*U;
    V = det(V)*V;

    Rzp  = [0 -1 0; 1 0 0; 0 0 1];
    Rzm  = Rzp';

    transfPairs(1).R = U*Rzp'*V';  THat{1} = U*Rzp*S*U';
    transfPairs(2).R = U*Rzp'*V';  THat{2} = U*Rzm*S*U';
    transfPairs(3).R = U*Rzm'*V';  THat{3} = U*Rzm*S*U';
    transfPairs(4).R = U*Rzm'*V';  THat{4} = U*Rzp*S*U';
    nFrontPoints=zeros(1,4);
    for i=1:4
        transfPairs(i).T = cnormalize(vee3(THat{i}));
        transfPairs(i).lambda=epipolarTriangulateDepths(transfPairs(i).R,transfPairs(i).T,x1,x2);
        nFrontPoints(i) = sum(transfPairs(i).lambda(1,:)>0 & transfPairs(i).lambda(2,:)>0);
    end

    [~,idxMax]=max(nFrontPoints);

    R(:,:,iE)=transfPairs(idxMax).R;
    T(:,iE)=transfPairs(idxMax).T;
    lambda(:,:,iE)=transfPairs(idxMax).lambda;
    if flagProgressBar
        w(iE)
    end
end