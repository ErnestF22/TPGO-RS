% Code to reduce rank of solution to SDP using "Low-Rank Semidefinite
% Programming: Theory and Applications" by Alex Lemon, etc al.

% Factorizie current X = VV' (since X >=0 can always do)
% Since we only have one equality constraint for our QCQP, we should always
% be able to reduce to a rank 1 matrix (IE: trace(A*X')==1==X(1,1) )
X_approx = X;
A = zeros(6); A(1,1)=1; % constraint X(1,1)=1
C=[0 1/2 0 1/2 0 0;1/2 0 0 0 0 0;0 0 0 0 0 0;1/2 0 0 0 0 0;0 0 0 0 0 0;0 0 0 0 0 0]; % This is the trace of M_nn
syms d_var1 d_var2
for i=1:100
    % Diagonalize X
    [W,XD]=eig(X_approx);
    V=W*sqrt(XD)*W';
    
    % Choose random diagonal as delta matrix
%     Delta = diag(randi([-1,1],6,1));
    Delta = diag(randn(6,1));
    % Replace a random eigenvalue with d_var to solve
    Delta = sym(Delta);
    idx_replace_1 = randi([1,6],1);
    idx_replace_2 = randi([1,6],1);
    while idx_replace_2 == idx_replace_1
        idx_replace_2 = randi([1,6],1);
    end
    Delta(idx_replace_1,idx_replace_1)=d_var1;
    Delta(idx_replace_2,idx_replace_2)=d_var2;
    % Solve such that trace(A*(V*Delta*V')')=0 && trace( (V'*C*V)'*Delta ) == 0
    d_var_num = solve(...
        [trace(A*(V*Delta*V')')==0; trace( (V'*C*V)'*Delta)==0],...
        [d_var1;d_var2]);
    Delta(idx_replace_1,idx_replace_1)=d_var_num.d_var1;
    Delta(idx_replace_2,idx_replace_2)=d_var_num.d_var2;
    Delta = double(Delta); % trace(A*(V*Delta*V')') should equal 0 now
    
    % Find alpha
    [max_eVal, idx_eVal] = max(abs(diag(Delta)));
    alpha = -1/Delta(idx_eVal,idx_eVal);
    
    % Update X
    X_approx = V*(eye(6) + alpha*Delta)*V';
    
    if rank(X_approx) == 1
        break;
    end
end

% Due to numeric approx, result maybe have complex parts that are 0
X_approx = real(X_approx);

% Define new contraction metric
m1 = X_approx(1,2); m2 = X_approx(1,3); m3 = X_approx(1,4); 
m5 = X_approx(1,5); m6 = X(1,6); m4 = 1;
M_nn_new = [m1 m2 m6;m2 m3 m5;m6 m5 m4]