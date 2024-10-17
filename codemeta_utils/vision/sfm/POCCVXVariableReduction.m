function POCCVXVariableReduction
cvx_quiet(false)
cvx_begin
    variable x(2)
    minimize (norm(x))
    subject to
        x(1)==1;
cvx_end
