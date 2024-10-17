function l=bearingCluster_scalesSolutionFromL(L,b,idxFix,scaleFix)
if exist('idxFix','var') && ~isempty(idxFix)
    Lij=L(idxFix,:);
    N=null(Lij);

    l=scaleFix*L*Lij'/norm(Lij')^2+L*N*b;
else
    l=L*b;
end

