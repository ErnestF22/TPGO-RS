function locSegSegmentationPair_test
%   1   2   3   4   5
mA=[1   1   2   2   0];
mB=[0   0   2   1   1];
mD=[1   1   2   0   0];
mI=[0   0   0   0   0];

mDTest=locSegSegmentationPairDiff(mA,mB);
check('mD',mA,mB,mD,mDTest);
mITest=locSegSegmentationPairIntersection(mA,mB);
check('mI',mA,mB,mI,mITest);

function check(testString,mA,mB,mR,mRTest)
if all(mRTest==mR)
    fprintf('\t%s OK\n',testString)
else
    fprintf('\t%s errors:\n',testString)
    disp([mR;mRTest]);
end
