function sfm_rawRotationsAdjust_test
NR=3;
K=randn(3);

for flagLeft=[false true]
    RiTruth=rot_randn(eye(3),[],NR);
    Ri=multiprod(RiTruth,K);
    if flagLeft
        Ri=multitransp(Ri);
        RiAdjusted=sfm_rawRotationAdjust(Ri,'left');
    else
        RiAdjusted=sfm_rawRotationAdjust(Ri,'right');
    end
    for iR=1:NR
        disp(RiAdjusted(:,:,iR)*RiAdjusted(:,:,iR)')
    end
end
