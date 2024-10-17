function errorsAggregated=multiCalibration_testRealData_aggregateErrors(errors)
fieldNames=fields(errors{1});
for iField=1:length(fieldNames)
    fn=fieldNames{iField};
    e=cellfun(@(x) x.(fn),errors,'UniformOutput',false);
    errorsAggregated.(fn)=cat(1,e{:});
end

