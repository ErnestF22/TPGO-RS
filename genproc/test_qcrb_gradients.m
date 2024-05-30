function test_qcrb_gradients
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here


load('Qbnn_data/Qbnn_R.mat', 'R')
load('Qbnn_data/Qbnn_T.mat', 'T')
load('Qbnn_data/Qbnn_edges.mat', 'edges')
load('Qbnn_data/x_rs.mat', 'x_rs')
load('Qbnn_data/Tijs.mat', 'Tijs')

node_deg = 2;

qc0=make_rand_stiefel_3d_array(4,4,1);
if det(qc0) < 0
    qc0(:,1) = - qc0(:,1);
end
rb0=make_rand_stiefel_3d_array(2,2,1);
if det(rb0) < 0
    rb0(:,1) = - rb0(:,1);
end

disp("check_is_rotation(qc0)")
disp(check_is_rotation(qc0))
disp("check_is_rotation(rb0)")
disp(check_is_rotation(rb0))


% [Rt,dRt,~,~,~,ddRt]=qcrb_geodFun(R0);
[qcRt, dqcRt, ~, ~, ~, ddqcRt] = rot_geodFun(qc0, []);
[rbRt, drbRt, ~, ~, ~, ddrbRt] = rot_geodFun(rb0, []);


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

Qa_i_1 = [0.00282804903025583, 0.0399962827864656];
Qa_i_2 = [0.0626139221207881, -0.997232067403872];

% test qc
f=@(t) mycost_qc(qcRt(t), rb0, node_deg, Qa, Qcd_i, Ri);
egradf=@(t) myegrad_qc(qcRt(t), rb0, node_deg, ...
    Qa, Qcd_i, Qa_i_1, Qa_i_2, Ri);
df=@(t) stiefel_metric([],egradf(t),dqcRt(t));
funCheckDer(f,df)

% test rb
f=@(t) mycost_rb(qc0, rbRt(t), node_deg, Qa, Qcd_i, Ri);
egradf=@(t) myegrad_rb(qc0, rbRt(t), node_deg, ...
    Qa, Qcd_i, Qa_i_1, Qa_i_2, Ri);
df=@(t) stiefel_metric([],egradf(t),drbRt(t));
funCheckDer(f,df)



end %file function

function bool_is_rot = check_is_rotation(A)
    if (det(A) < 1-1e-6 || det(A) > 1+1e-6)
        disp("check_is_rotation() -> determinant = 1 check failed")
        fprintf("det = %g\n", det(A))
        bool_is_rot = boolean(0);
        return;
    end
    if max(abs(A * A' - A' * A)) > 1e-6
        disp("check_is_rotation() -> A*A' == A'*A check failed")
        bool_is_rot = boolean(0);
        return;
    end
    if max(abs(A * A' - eye(size(A)))) > 1e-6
        disp("check_is_rotation() -> A*A' == I check failed")
        bool_is_rot = boolean(0);
        return;
    end
    bool_is_rot = boolean(1);
end


function c_out = mycost_qc(x, Rb, node_deg, Qa, Qcd_i, Ri, Ti, Ti_tilde)
Qcdd_i = x;
rb_i = Rb;
p = size(Ri,1);
d = size(Ri,2);
% A
P = [zeros(p-d,d), eye(p-d)];
A = P * Qcdd_i' * Qcd_i * Ri;
% B
%     Qb_i = eye(p);
%     Qb_i(p-1:end, p-1:end) = rb_i;
Qb_i = blkdiag(eye(p-node_deg), rb_i); %node_deg = 2
B = Qcdd_i' - Qa' * Qb_i * Qa;

% sum of norms
c_out = norm(A(:))^2 + norm(B(:))^2;
end


function grad_out = myegrad_qc(x, Rb, node_deg, ...
    Qa, Qcd_i, Qa_i_1, Qa_i_2, Ri)
Qcdd_i = x;
rb_i = Rb;
p = size(Ri,1);
d = size(Ri,2);
% A
P = [zeros(p-d,d), eye(p-d)];
% Qcd_i = eye(p);
A = P * Qcdd_i' * Qcd_i * Ri;
% B
%     Qb_i = eye(p);
%     Qb_i(p-1:end, p-1:end) = rb_i;
Qb_i = blkdiag(eye(p-node_deg), rb_i); %node_deg = 2
B = Qcdd_i' - Qa' * Qb_i * Qa;
% grad out (struct)
% Qa1 = Qa(1:2, :);
Qa2 = Qa(3:end, :);
grad_out = 2 * Qcd_i * Ri * A' * P + B;
%     grad_out.rb = - Q2 * B * Q2';
% grad_out.rb = - 2 * Qa2 * B * Qa2';

end 

function c_out = mycost_rb(qc, x, node_deg, Qa, Qcd_i, Ri, Ti, Ti_tilde)
Qcdd_i = qc;
rb_i = x;
p = size(Ri,1);
d = size(Ri,2);
% A
P = [zeros(p-d,d), eye(p-d)];
A = P * Qcdd_i' * Qcd_i * Ri;
% B
%     Qb_i = eye(p);
%     Qb_i(p-1:end, p-1:end) = rb_i;
Qb_i = blkdiag(eye(p-node_deg), rb_i); %node_deg = 2
B = Qcdd_i' - Qa' * Qb_i * Qa;

% sum of norms
c_out = norm(A(:))^2 + norm(B(:))^2;
end

function grad_out = myegrad_rb(qc, x, node_deg, ...
    Qa, Qcd_i, Qa_i_1, Qa_i_2, Ri, T_i_j1_j2, T_i_j1_j2_tilde)
Qcdd_i = qc;
rb_i = x;
p = size(Ri,1);
d = size(Ri,2);
% A
P = [zeros(p-d,d), eye(p-d)];
% Qcd_i = eye(p);
A = P * Qcdd_i' * Qcd_i * Ri;
% B
%     Qb_i = eye(p);
%     Qb_i(p-1:end, p-1:end) = rb_i;
Qb_i = blkdiag(eye(p-node_deg), rb_i); %node_deg = 2
B = Qcdd_i' - Qa' * Qb_i * Qa;
% grad out (struct)
% Qa1 = Qa(1:2, :);
Qa2 = Qa(3:end, :);
% grad_out.qc = 2 * Qcd_i * Ri * A' * P + B;
%     grad_out.rb = - Q2 * B * Q2';
grad_out = - 2 * Qa2 * B * Qa2';


end 
