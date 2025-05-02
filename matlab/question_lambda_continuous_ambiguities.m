function question_lambda_continuous_ambiguities

%% first, we want to generate a geodesic curve on Stiefel

n = 4;
k = 3;
A = rand(k,k);
A_skew = 0.5 * (A - A');

B = rand(1, k);
Z = zeros(n-k);
X = [A_skew, -B'; B, Z];


Delta = [eye(k); zeros(1,k)];
Theta = make_rand_stiefel_3d_array(n, n, 1);

S = Theta * Delta;

ts = 0:0.01:1;
geod = zeros(n,k,length(ts));

idx = 1;

for t = ts
    etx = expm(t * X);
    newm = Theta * etx * Delta; %new matrix belonging to the geodesic
    geod(:,:,idx) = newm;
    % disp("is_equal_floats(newm' * newm, eye(k,k))")
    % disp(is_equal_floats(newm' * newm, eye(k,k)))
    idx = idx + 1;
end

disp("Check if S is actually the starting point")
% disp("is_equal_floats(S, geod(:,:,1))")
disp(is_equal_floats(S, geod(:,:,1)))

%% gradient

geod_grad = zeros(n,k,length(ts));
for t = ts
    etx = expm(t * X);
    newm = Theta * X * etx * Delta; %new matrix belonging to the geodesic
    geod_grad(:,:,idx) = newm;
    % disp("is_equal_floats(newm' * newm, eye(k,k))")
    % disp(is_equal_floats(newm' * newm, eye(k,k)))
    idx = idx + 1;
end

%% funCheckDer() for Stiefel geodesic

fx = @(t) Theta * expm(t * X) * Delta;
v0Stief = stiefel_randTangentNormVector(S);
dfx = @(t) Theta * X * expm(t * X) * Delta;
% funCheckDer(fx, dfx)

% function H=stiefel_eucl2RiemGrad(Y,H)
%Projects H onto the tangent space of the Stiefel manifold at Y

%% Check \dot{R_i} = R_i A + R_i^{\bot} B
% Hp) v0Stief = \dot{R_i}
% v0Stief = S A + R_i^{\bot} B
% A skew-symm
% B arbitrary
% R_i^{\bot} basis for the orthogonal complement of R_i



%% (37)
lambdas = rand(2,1);
Lambda = diag(lambdas);

[Lambdat,dLambdat,~,vLambda0,ddLambdat] = real_randGeodFun(lambdas); % missing output parameter would just be equal to Lambda

%% matrix-form

resetRands();
B = rand;
A = skew3(rand(3,1));
E = eye(3);



%% keyboard()
keyboard()

end %file function
