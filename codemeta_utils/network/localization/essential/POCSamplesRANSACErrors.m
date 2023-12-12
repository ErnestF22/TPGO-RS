idxPair=5;
idx=7;
a=cellfun(@(x) x.errR(:,idx), errors(idxPair,:),'UniformOutput',false);
a=cat(2,a{:});
plot(mean(a,2));
legend(errors{idxPair,1}.legendTextR{idx})
