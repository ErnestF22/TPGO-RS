function d=grassman_metric(~,H1,H2)
NH1=size(H1,3);
NH2=size(H2,3);
d=zeros(NH1,NH2);
for iH1=1:NH1
    for iH2=1:NH2
        d(iH1,iH2)=trace(H1(:,:,iH1)'*H2(:,:,iH2));
    end
end
