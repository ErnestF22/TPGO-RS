function POCmanoptTrifocalRotations

%%
dataset=1;
switch dataset
    case 1
        data=trifocal_dataset(6);
        RTruth=cat(3,data.R(:,:,1)*data.R(:,:,2)',data.R(:,:,1)*data.R(:,:,3)');
        n0=data.l2d(:,:,1);
        n1=data.l2d(:,:,2);
        n2=data.l2d(:,:,3);
    case 2
        Nn=2;
        n0=squeeze(sphere_randn([],[],Nn));
        n1=squeeze(sphere_randn([],[],Nn));
        n2=squeeze(sphere_randn([],[],Nn));
end

% n=3;
% 
% AB=randn(n);
% ABt=AB.';

m=rotationsfactory(3,2);
problem.M=m;
%%

% Define the problem cost function and its derivatives.
% problem.cost = @(X) -X(:).'*ABt(:);
% egrad = @(X) -ABt;
% problem.grad = @(X) m.egrad2rgrad(X, egrad(X));
% 
% ehess = @(X, S) m.zerovec(X);
% problem.hess = @(X, S) m.ehess2rhess(X, egrad(X), ehess(X, S), S);
problem.cost=@(R) 100*cost(R,n0,n1,n2);
problem.grad=@(R) 100*grad(R,n0,n1,n2)/2;
problem.hess=@(R,S) 100*hess(R,S,n0,n1,n2)/2;
%%

%checkgradient(problem)
%checkhessian(problem)

ROpt=trustregions(problem);
disp(reshape(RTruth,3,[]))
disp(problem.cost(RTruth))
disp(reshape(ROpt,3,[]))
disp(problem.cost(ROpt))


function f=cost(R,n0,n1,n2)
f=trifocalRotationCost(R(:,:,1),R(:,:,2),n0,n1,n2);
f=sum(f);

function gf=grad(R,n0,n1,n2)
[~,gf]=trifocalRotationCost(R(:,:,1),R(:,:,2),n0,n1,n2);
gf=rot_hat(eye(3),reshape(sum(gf,2),3,2));

function HfS=hess(R,S,n0,n1,n2)
[~,~,Dgf]=trifocalRotationCost(R(:,:,1),R(:,:,2),n0,n1,n2);
Dgf=sum(Dgf,3);
SVec=rot_vee(eye(3),S);
Hf=(Dgf+Dgf')/2;
HfS=rot_hat(eye(3),reshape(Hf*SVec(:),3,2));

