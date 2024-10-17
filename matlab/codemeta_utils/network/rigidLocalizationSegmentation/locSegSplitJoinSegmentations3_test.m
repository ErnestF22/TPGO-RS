function locSegSplitJoinSegmentations3_test
NTests=10;
for iTest=1:NTests
    switch iTest
    %           1   2   3   4   5
        case 1
            mA=[1   1   2   2   0];
            mB=[0   0   2   1   1];
            mR=[1   1   2   0   0];
        case 2
            mA=[0   2   2   1   1];
            mB=[0   0   2   1   1];
            mR=[0   0   2   1   1];
        case 3
            mA=[0   0   2   1   1];
            mB=[1   1   2   2   0];
            mR=[0   0   2   1   1];
        case 4
            mA=[0   0   2   1   1];
            mB=[0   2   2   1   1];
            mR=[0   0   2   1   1];
        case 5
            mA=[1   1   2   2   0];
            mB=[0   0   2   2   1];
            mR=[1   1   2   2   0];
        case 6
            mA=[0   2   2   1   1];
            mB=[0   0   2   2   1];
            mR=[0   0   2   2   1;
                0   2   2   1   0];
        case 7
            mA=[0   0   2   2   1];
            mB=[0   2   2   1   1];
            mR=[0   0   2   2   1];
        case 8
            mA=[0   0   2   2   1];
            mB=[1   1   2   2   0];
            mR=[0   0   2   2   1];
        case 9
            mA=[1   2   2   0];
            mB=[2   1   0   2];
            mR=[1   0   2   0;
                1   1   0   0];
        case 10
            mA=[1   1   1   1   2   0];
            mB=[0   2   1   1   1   1];
            mR=[1   2   0   0   0   0;  %mOut
                0   2   1   1   2   0]; %mIn  
    %     case 
    %         mA=[];
    %         mB=[];
    %         mR=[];
        otherwise
            continue
    end
    fprintf('Test %d\n',iTest)
    check(mA,mB,mR);
end
function check(mA,mB,mR)
mRTest=locSegSplitJoinSegmentations3(mA,mB);
mR=sortrows(mR);
mRTest=sortrows(mRTest);
if numel(mRTest)==numel(mR)
    if all(mRTest==mR)
        fprintf('\tOK\n')
    else
        fprintf('\tErrors: different assignments\n')
        disp([mR;mRTest]);
    end
else
    fprintf('\tError: different dimensions')
    disp([size(mR,1) size(mRTest,1)])
end
