% POC test code for using SDP
cvx_begin sdp
    variable X(3,3) semidefinite
    Q = [0 1 2;1 0 3;2 3 0];
    
    minimize ( trace(Q*X) )
    subject to
%         trace(C*X) <= 0
        X(1,1) == 1
        X(2,2) == 1
        X(3,3) == 1
cvx_end
fprintf('Result from SDP\n');
trace(Q*X)
X

%% Solve its dual
cvx_begin sdp
    variable L(3,3) diagonal
    Q = [0 1 2;1 0 3;2 3 0];
    
    maximize ( trace(L) )
    subject to
        Q >= L
cvx_end
fprintf('Result from dual\n');
trace(L)
L