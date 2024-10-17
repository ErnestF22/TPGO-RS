function f=medianFieldRegularizedEnergy(uMedian,uReg,idxList,kList,lambda)
f=0;
Nu=numel(uMedian);
for iu=1:Nu
    f=f+medianRegularizedEnergy(uMedian(iu),uMedian(idxList(iu,1:kList(iu))),lambda(2),uReg(iu),lambda(1));
end
