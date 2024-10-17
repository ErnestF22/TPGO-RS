%Combine rotations in a cycle
%function RCombined=sfm_rawRotationsCycleCombine_IndexSign(R,idxc,signc)
%The cycle is given with the indeces of the edges (in order) and 
function RCombined=sfm_rawRotationsCycleCombine_IndexSign(R,idxc,signc)
RCombined=eye(3);
for ic=1:length(idxc)
    if signc(ic)>0
        RCombined=RCombined*R(:,:,idxc(ic));
    else
        RCombined=RCombined*R(:,:,idxc(ic))';
    end
end
