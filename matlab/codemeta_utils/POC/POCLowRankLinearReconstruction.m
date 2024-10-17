function POCLowRankLinearReconstruction
resetRands()
n=7;
k=3;
d=[k*n,n];
nbElements=prod(d);
nbMeasurements=max(round(nbElements*0.8),3*d(2));
A=randn(d(1),k)*randn(k,d(2));
fprintf('# measurements: %d, theoretical min # measurements: %d, # entries: %d\n',nbMeasurements,k^2*n+k*n,prod(d))

EDiag=zeros([size(A),d(2)]);

idxE=reshape(1:k*d(2),k,d(2));
idxDiag=idxE+(0:d(2)-1)*d(1);
for id=1:d(2)
    for ik=1:k
        Ei=zeros(d);
        Ei(idxDiag(ik,id))=1;
        EDiag(:,:,idxE(ik,id))=Ei;
    end
end
nbMeasurementsRandom=nbMeasurements-k*d(2);
Omega=randsample(nbElements,nbMeasurementsRandom);
EMeasurements=zeros([d nbMeasurementsRandom]);
EMeasurements(Omega+(0:nbMeasurementsRandom-1)'*nbElements)=1;

E=cat(3,EDiag,EMeasurements);

cvx_begin quiet
    variable AEstimated(d(1),d(2))
    f=0.000001*norm_nuc(AEstimated);
    for iMeasurement=1:nbMeasurements
        f=f+abs(sum(sum((A-AEstimated).*E(:,:,iMeasurement),1),2));
    end
    minimize f
cvx_end
disp(sum(E,3))
disp(full([AEstimated-A]))
svd(AEstimated)