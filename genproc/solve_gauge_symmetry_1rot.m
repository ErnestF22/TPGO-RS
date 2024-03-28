nrs = 4; %change to nrs = 5 to test it
d = 3;
N = 1;


% R_i = make_rand_stiefel_3d_array(nrs, d, N);
% T_i = 10 * rand(nrs, N);

% POCRotateToMinimizeLastEntries() ...

R = zeros(nrs, d, 1);
R(:,:,1) =[
  -0.812728367573879   0.364762014232424   0.418766102871191;
   0.079417601237271  -0.664534191278519   0.502715344803537;
   0.478086669628382   0.364716187242740   0.738515371539764;
   0.323417039560469   0.540670957568491  -0.162810562263695];

% T = zeros(nrs,N);
T = [5.52989200497386	-3.23868362853368	-7.57164584572133	-0.347137118025300	6.40077586933894;
-4.40050299467541	-2.74607885107208	1.60648817904566	4.71944333560070	-0.729181940525049;
-5.73635735107686	-7.22243705954438	0.0545123672074450	6.10779249344889	4.55609031254877;
-0.423968727476036	-0.177903011520358	0.614132793616925	-0.874411109442508	0.862364178906416];
T = T(:,1);

%%

x = [matStackH(R)];
[Q_svd,~,~]=svd(x*x');
disp('x=')
disp(x)
disp('Q''*x=')
disp(Q_svd'*x)

% transf_out = RT2G(matUnstackH(x), T(1:d, 1:N));

%% by hand (Manopt)
problem.M = rotationsfactory(nrs);
P_zeros = zeros(nrs-d, d);
P = [P_zeros, eye(nrs-d)];

problem.cost = @(Q) norm(P*Q*x);
problem = manoptAD(problem);

% Q_start = [
%     0.2577    0.1364    0.9560    0.0314
%     0.9365   -0.2768   -0.2116   -0.0406
%    -0.2375   -0.9480    0.2015   -0.0655
%     0.0144   -0.0779   -0.0255    0.9965];


% [Q_out, Q_cost, ~, ~] = trustregions(problem);

disp("Q_out")
disp(Q_out)

disp("Q_out * x")
disp(Q_out * x)

% disp("multidet(matUnstackH(Q_out * x, d))");
% Xtrue = matUnstackH(Q_out * x, d);
% disp(multidet(Xtrue(1:d, 1:d, :)));

%% POC

x = R;
Q_gauge = POCRotateToMinimizeLastEntries(x);
disp('x=')
disp(x)
disp('Q_gauge*x=')
disp(Q_gauge*x)


