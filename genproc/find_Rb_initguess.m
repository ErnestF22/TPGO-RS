function Rb = find_Rb_initguess(node_deg,p,d,Qa_i,R_i)
%linsyst + procrustes initguess for Rb

% Q_1 \in \real{(p-d)\times 2} as the first 2 columns of P Q_{a,i}'
% Q_2 \in \real{(p-d)\times (p-2)}$ as the last $p-2$ columns of P Q_{a,i}'
% R_1 \in \real{2\times 3} as the first 2 rows of Q_{a,i} \hat{Q}_{c,i}' \tilde{\tilde{R}}_i
% R_2 \in\real{(p-2)\times 3} as the last $p-2$ rows of Q_{a,i} \hat{Q}_{c,i}' \tilde{\tilde{R}}_i

P = [zeros(p-d,d), eye(p-d)];
Qcd_i = eye(p); %??
QQ = P * Qa_i';
RR = Qa_i * Qcd_i' * R_i;
Q1 = QQ(:, 1:2);
Q2 = QQ(:, 3:end);
R1 = RR(1:2,:);
R2 = RR(3:end,:);

A_ls = R2'; %ls = lin. syst.
b_ls = - R1' * Q1';
Ab_transp =  A_ls \ b_ls;
Ab = Ab_transp';
disp('Ab')
disp(Ab)


%Procrustes to find Rb
a_procr = Q2;
b_procr = Ab;
Rb = procrustes_kabsch(a_procr, b_procr, node_deg);




end %file function
