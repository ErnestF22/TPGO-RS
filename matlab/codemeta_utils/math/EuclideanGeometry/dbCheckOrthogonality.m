%checks if x has orthonormal columns. If not, go into the debugger
function dbCheckOrthogonality(x)
p=x'*x;
%average distance per entry w.r.t. p being an identity
avgErr=sum(sum(abs(p-eye(size(p)))))/numel(p);
if avgErr>1e-10
    disp(['Average, per-entry error on inner products giving the identity: ' num2str(avgErr)])
    keyboard
end

    