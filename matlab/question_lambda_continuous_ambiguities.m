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

%% keyboard()
% keyboard()

end %file function
