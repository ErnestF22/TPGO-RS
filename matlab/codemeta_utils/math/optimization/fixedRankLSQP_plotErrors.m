function fixedRankLSQP_plotErrors(output)
recognizedFields.linearResiduals={@semilogy,'LS cost'};
recognizedFields.fixedRankResiduals={@semilogy,'Residual rank constraint'};
recognizedFields.inequalitiesResiduals={@plott,'Residuals inequality constraints'};
recognizedFields.solutionSvals={@semilogyt,'Svals of X'};
recognizedFields.solutionReferenceResiduals={@plot,'Distance from reference solution'};

%list fields and filter out those that cannot be recognized
outputFields=fieldnames(output);
outputFields=outputFields(ismember(outputFields,fieldnames(recognizedFields)));
nbFieldsDisplayed=length(outputFields);

%iterate over fields to be displayed
for iField=1:nbFieldsDisplayed
    subplot(nbFieldsDisplayed,1,iField)
    of=outputFields{iField};
    af=recognizedFields.(of);
    %call action for that field
    af{1}(output.(of))
    title(af{2})
end

function plott(x)
plot(x')

function semilogyt(x)
semilogy(x')
