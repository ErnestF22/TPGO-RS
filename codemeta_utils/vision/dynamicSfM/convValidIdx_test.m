function convValidIdx_test
for ny=9:11
    for nh=5:7
        display(nh)
        y=randn(ny,1);
        h=randn(nh,1);
        yhValid=conv(y,h,'valid');
        disp('Errors:')
        for shape={'valid','same','full'}
            yh=conv(y,h,shape{1});
            idxValid=convValidIdx(ny,nh,shape{1});
            disp(norm(yh(idxValid)-yhValid))
        end
    end
end
