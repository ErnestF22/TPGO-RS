function [Ri,Ti]=lowRankLocalization_solution_extractProjection(W,WInfo)
[U,S,V]=svd(W,'econ');
[Ri,K]=sfm_rawRotationsAdjust(U(:,1:3));
Ri=rot_proj(matUnstack(Ri));
Ti=K\(S(1:3,1:3)*V(:,1:3)');
if exist("WInfo","var")
    if WInfo.flagRotationAugmented
        Ti=Ti(:,1:end-WInfo.dimAmbient);
    end
    if isfield(WInfo,'iNodeFix')
        RiFix=Ri(:,:,WInfo.iNodeFix);
        Ri=multiprod(Ri,RiFix');
        Ti=RiFix*Ti;
    end
end

%The left factor in W is defined with transposed rotations
Ri=multitransp(Ri);
