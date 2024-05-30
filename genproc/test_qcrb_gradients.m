function test_qcrb_gradients
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here


load('Qbnn_data/Qbnn_R.mat', 'R')
load('Qbnn_data/Qbnn_T.mat', 'T')
load('Qbnn_data/Qbnn_edges.mat', 'edges')
load('Qbnn_data/x_rs.mat', 'x_rs')
load('Qbnn_data/Tijs.mat', 'Tijs')

node_deg = 2;

R10=stiefel_randn(4,4,1);
R20=stiefel_randn(2,2,1);

[R1,dR1,~,~,~,ddR1]=rot_geodFun(R10, []);
[R2,dR2,~,~,~,ddR2]=rot_geodFun(R20, []);


node_deg = 2;
i = 1;
j1 = 2;
j2 = 3;
R_i = R(:,:,1);
idx_i_j1 = find(ismember(edges, [i, j1], "rows"));
idx_i_j2 = find(ismember(edges, [i, j2], "rows"));

T_i_j1_j2 = [Tijs(:,idx_i_j1), Tijs(:,idx_i_j2)];
T_i_j1_j2_tilde = - [T(:,i) - T(:,j1), T(:,i) - T(:,j2)]; % !! -

Qa = [-1.64302828263923e-16	0.106320705792361	-0.889261908456345	-0.444869830049638;
    0.993981726033715	-0.108890962053455	-0.0116398298358148	-0.00275700117918123;
    2.81044311504080e-17	0.0248892039991507	-0.444885487710104	0.895241548605309;
    0.109546010018785	0.988038052621039	0.105615696536622	0.0250160529834866];
Qcd_i = eye(4);
Ri = [0.0864141652363848	0.0420429561776974	-0.993397284823312;
    0.256372166990434	0.964519148326362	0.0626995624731963;
    -0.962094412651845	0.260563503917528	-0.0701748073529320;
    0.0343546966691522	-0.00647014480909629	0.0656208487016165];

xf.qc = @(t) R1(t);
xf.rb = @(t) R2(t);
f=@(t) mycost_qcrb(xf(t), node_deg, Qa, Qcd_i, Ri, T_i_j1_j2, T_i_j1_j2_tilde);
egradf=@(t) myegrad_qcrb(xf(t), deg_i, Qa_i, Qcd_i, Qa_i_1, Qa_i_2, Ri, T_i_j1_j2, T_i_j1_j2_tilde);
df=@(t) sum(stiefel_metric(c,egradf(t),dc(t)));
funCheckDer(f,df)



end %file function

