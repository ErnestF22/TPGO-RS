function [e, grade, hessOpe]=epipolarCostFromE_SampsonSq(E,x1,x2,varargin)
flagSymmetricHess=false;

%optional parameters
ivarargin=1;
while(ivarargin<=length(varargin))
    switch(lower(varargin{ivarargin}))
        case 'symmetrichess'
            flagSymmetricHess=true;
        case 'flagsymmetrichess'
            ivarargin=ivarargin+1;
            flagSymmetricHess=varargin{ivarargin};
    end
    ivarargin=ivarargin+1;
end

x1=homogeneous(x1,3);
x2=homogeneous(x2,3);

flagComputeGrad=nargout>1;
flagComputeHess=nargout>2;    

%Residuals (arguments of the squares)
numeratorRes=sum(x1.*(E*x2));
denominator1Res=E(1:2,:)*x2;
denominator2Res=E(:,1:2)'*x1;

%Numerator and denominator of the Sampson cost
numerator=numeratorRes.^2;
denominator1=sum(denominator1Res.^2);
denominator2=sum(denominator2Res.^2);
denominator=denominator1+denominator2;

%Sampson cost
e=sum(numerator./denominator);

if flagComputeGrad
    Nx=size(x1,2);
    z=zeros(3,1,Nx);
    
    %reshape variables into a form more amenable for vectorized computations
    x1=permute(x1,[1 3 2]);
    x2=permute(x2,[1 3 2]);
    numeratorRes=shiftdim(numeratorRes,-1);
    denominator1Res=permute(denominator1Res,[1 3 2]);
    denominator2Res=permute(denominator2Res,[1 3 2]);
    numerator=shiftdim(numerator,-1);
    denominator=shiftdim(denominator,-1);
    
    gradNumeratorRes=multiprod(x1,multitransp(x2));

    %Gradient of the numerator and the denominator
    gradNumerator=2*multiprod(numeratorRes,gradNumeratorRes);
    %we can compute the gradient of the cost without going through all the
    %gradients of the individual residuals
    gradDenominator=2*(multitransp([multiprod(x2,multitransp(denominator1Res)),z])...
        +[multiprod(x1,multitransp(denominator2Res)),z]);

    %Gradient of the Sampson cost
    gradeNumerator=multiprod(denominator,gradNumerator)-multiprod(numerator,gradDenominator);
    gradeDenominator=denominator.^2;
    grade=sum(multidiv(gradeNumerator,gradeDenominator),3);

    if flagComputeHess
        %in order to compute the hessian operator for the denominator, we
        %will need the gradients of the individual residuals
        graddenominator11Res=multitransp([x2 z z]);
        graddenominator12Res=multitransp([z x2 z]);
        graddenominator21Res=[x1 z z];
        graddenominator22Res=[z x1 z];

        dgradNumerator=@(dE) 2*multiprod(multitraceprod(dE',gradNumeratorRes),gradNumeratorRes); %symmetric

        dgradDenominator=@(dE) 2*(...
            +multiprod(multitraceprod(dE',graddenominator11Res),graddenominator11Res)...
            +multiprod(multitraceprod(dE',graddenominator12Res),graddenominator12Res)...
            +multiprod(multitraceprod(dE',graddenominator21Res),graddenominator21Res)...
            +multiprod(multitraceprod(dE',graddenominator22Res),graddenominator22Res)...
            ); %symmetric
        
        if ~flagSymmetricHess
            %second derivative
            dgradeNumerator=@(dE) mprtrpr(dE',gradDenominator,gradNumerator)+multiprod(denominator,dgradNumerator(dE))...
                -mprtrpr(dE',gradNumerator,gradDenominator)-multiprod(numerator,dgradDenominator(dE));
            %dgradeNumerator(dE): matrix
            %gradeDenominator: scalar
            %gradeNumerator: matrix
            %dgradeDenominator(dE): scalar
            hessOpe=@(dE) sum(multidiv(...
                (multiprod(dgradeNumerator(dE),gradeDenominator)...
                -2*multiprod(denominator,mprtrpr(dE',gradDenominator,gradeNumerator))),...
                gradeDenominator.^2),3);
        else
            %symmetric hessian
            dgradeNumerator=@(dE) mprtrprSym(dE',gradDenominator,gradNumerator)+multiprod(denominator,dgradNumerator(dE))...
                -mprtrprSym(dE',gradNumerator,gradDenominator)-multiprod(numerator,dgradDenominator(dE));

            hessOpe=@(dE) sum(multidiv(...
                (multiprod(dgradeNumerator(dE),gradeDenominator)...
                -2*multiprod(denominator,mprtrprSym(dE',gradDenominator,gradeNumerator))),...
                gradeDenominator.^2),3);
        end            

    end
end


function C=multitraceprod(A,B)
C=shiftdim(multitrace(multiprod(A,B)),-2);

function D=mprtrpr(A,B,C)
D=multiprod(multitraceprod(A,B),C);

function D=mprtrprSym(A,B,C)
D=(multiprod(multitraceprod(A,B),C)+multiprod(multitraceprod(A,C),B))/2;
