function genAD_myfun
x=adigatorCreateDerivInput([3 1],'x');
k=adigatorCreateAuxInput([3 1]);
adigator('myfun',{x,k},'myderiv')
