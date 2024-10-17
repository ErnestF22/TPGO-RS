function fixedRankLSQP_test
resetRands(0)
sz=[20,10];
k=3;
nbConstraints=30;
nbInequalities=7;
flagCheckConstraints=true;
methodInequalities='none';
methodConstraints='random';
methodPerturbation='empty';

%"true" reference solution
XReference=randn(sz(1),k)*randn(k,sz(2));

%linear constraints
switch methodConstraints
    case 'random'
        A=randn(nbConstraints,prod(sz));
        b=A*vec(XReference);
    case 'projective'
        A=randn(nbConstraints,prod(sz))*orthComplementProjector(vec(XReference));
        b=zeros(nbConstraints,1);
        if ~strcmp(methodInequalities,'random')
            error('With projective LS constraints, inequalities are needed')
        end
    otherwise
        error('methodConstraints not recognized')
end


switch methodInequalities
    case 'random'
        %random inequalities, one active at the true solution
        C=randn(nbInequalities,prod(sz));
        d=C*vec(XReference);
        idxActive=randi(1,nbInequalities);
        idxChange=[1:idxActive-1 (idxActive+1):nbInequalities]';
        d(idxChange)=d(idxChange)-rand(nbInequalities-1,1);
    case 'null'
        %all inequalities are trivially true for any solution
        C=zeros(nbInequalities,prod(sz));
        d=-ones(nbInequalities,1);
    case 'none'
        %no inequalities are passed
        C=[];
        d=[];
    otherwise
        error('methodInequalities not recognized')
end

%add perturbation to initial condition
switch methodPerturbation
    case 'none'
        XInitial=XReference;
    case 'subspace'
        XInitial=XReference+10*reshape(A'*randn(nbConstraints,1),sz(1),sz(2));
    case 'full'
        XInitial=XReference+10*randn(size(XReference));
    case 'empty'
        %Ask for automatic initialization
        XInitial=[];
    otherwise
        error('methodPerturbation not recognized')
end        


if flagCheckConstraints
    disp('Linear residuals at reference solution should be zero')
    disp((A*vec(XReference)-b)')
    if ~isempty(C)
        disp('Linear inequalities at reference solution should be zero or positive')
        disp((C*vec(XReference)-d)')
    end
end

[XEstimated,U,output]=fixedRankLSQP(A,b,C,d,sz,k,...
    'referenceSolution',XReference,'initialSolution',XInitial,...
    'maxIt',400);

figure(1)
fixedRankLSQP_plotErrors(output)
svd(U.FixedRank)
%keyboard
