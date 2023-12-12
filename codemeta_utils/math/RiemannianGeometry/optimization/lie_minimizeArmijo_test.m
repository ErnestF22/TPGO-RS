clear variables
randn('state',1)
testNumber=1;
d=6;

%init
x0=orth(randn(d));
if det(x0)<0
    [U,S,V]=svd(x0);
    x0=U*diag([ones(1,d-1) -1])*V';
end

switch testNumber
    case 1
        %simple linear functional
        lf=rot_funs();
        A=randn(d);
        f=@(x) 0.5*trace(A'*x);
        gradf=@(x) 2*rot_tangentProj(x,A);
        
        xclose=rot_proj(-A);
        
    case 2
        %minimize the trace of a symmetric matrix with given evals
        lf=rot_funs();
        L=diag([1; -1; zeros(d-2,1)]);
        f=@(x) 0.5*sum(diag(x*L*x').^2);
        euclGradf=@(x) 2*diag(diag((x*L*x')))*x*L;
        gradf=@(x) 2*rot_tangentProj(x,euclGradf(x));
    
    case 3
        %compute d/3-dimensional dominant subspace
        x0=orth(randn(d,ceil(d/3)));

        lf=grassman_funs();
        A=randn(d);
        A=A'*A;
        B=randn(size(x0));
        f=@(x) 0.5*trace(x'*A*x);
        euclGradf=@(x) A*x;
        gradf=@(x) 2*grassman_tangentProj(x,euclGradf(x));
    case 4
        %quadratic form + linear
        x0=orth(randn(d,ceil(d/1.5)));

        lf=grassman_funs();
        A=randn(d);
        A=A'*A;
        B=randn(size(x0));
        f=@(x) 0.5*trace(x'*A*x)+trace(rot_proj(-x'*B)*B'*x);
        euclGradf=@(x) A*x + B*rot_proj(-x'*B);
        gradf=@(x) grassman_tangentProj(x,euclGradf(x));
    case 5
        %like 3 but on Stiefel instead of Grassman
        x0=orth(randn(d,ceil(d/3)));

        lf=stiefel_funs();
        A=randn(d-ceil(d/3),d);
        A=A'*A;
        B=randn(size(x0));
        f=@(x) 0.5*trace(x'*A*x);
        euclGradf=@(x) A*x;
        gradf=@(x) stiefel_tangentProj(x,stiefel_eucl2RiemGrad(x,euclGradf(x)));
        
        s=svd(A);
        disp(['Expected minimum of the function ' num2str(0.5*sum(s(end-size(x0,2)+1:end)))])
    case 6
        %like 4 but on Stiefel instead of Grassman
        x0=orth(randn(d,ceil(d/3)));

        lf=stiefel_funs();
        A=randn(d);
        A=A'*A;
        B=randn(size(x0));
        f=@(x) 0.5*trace(x'*A*x)+trace(B'*x);
        euclGradf=@(x) A*x+B;
        gradf=@(x) stiefel_tangentProj(x,stiefel_eucl2RiemGrad(x,euclGradf(x)));
        
        [x0,R]=qr(A\B);
        x0=x0(:,1:size(B,2));
    case 7
        %data from stiefel_PCA

        lf=stiefel_funs();
        s=load('lie_minimizeArmijo_stiefel_PCA_data');
        x0=s.x0;
        A=s.A;
        B=zeros(size(x0));%s.B;
        
        f=@(x) -0.5*trace(x'*A*x)+trace(B'*x);
        euclGradf=@(x) -A*x+B;
        gradf=@(x) stiefel_tangentProj(x,stiefel_eucl2RiemGrad(x,euclGradf(x)));

%        [x0,R]=qr(-A\B);
%        x0=x0(:,1:size(B,2));
                
end

disp(datestr(now))
par={'exp' 'qr' 'polar'};
for ip=1:length(par)
    itersMax=3000;
    [x,errors]=lie_minimizeArmijo(lf,f,gradf,x0,'itersMax',itersMax,'method','conjgrad','armijoOpts',{'alpha',1,'beta',0.8,'sigma',1e-2},'retractionsOpts',{par{ip}});
    fevals(:,ip)=[errors.feval Inf(1,itersMax-length(errors.feval))];
    ts(:,ip)=[errors.t Inf(1,itersMax-length(errors.t))];
    if isnumeric(par{ip})
        legendtext{ip}=['par = ' num2str(par{ip})];
    else
        legendtext{ip}=['par = ' par{ip}];
    end
        
    disp(errors.exitInfo)
    disp(['Iterations performed: ' num2str(errors.itNum)])
    disp(['Last cost function: ' num2str(f(x))])
    disp(['Last magnitude of gradient: ' num2str(sqrt(lf.metric(x,gradf(x),gradf(x))))])
end

figure(1)
%semilogy(ts,fevals)
plot(ts,fevals)
legend(legendtext)

