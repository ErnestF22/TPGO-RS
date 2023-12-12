function localization_essential_test
resetRands();
%t_node=testNetworkBuildTestNetwork('methodInit','noisytruth',...
%    'varNoisyTruth',[0.5 0.5]);
maxIt=10000;

optsCommon={'maxit',maxIt,'progressBar'};
t_node=testNetworkBuildTestNetwork('methodInit','rand');
disp('Gradient descent')
[t_nodeOptGrad,outputGrad]=localization_essential_gradient(t_node,optsCommon{:},...
    'optsLieMinimize',{'getErrors'});
disp('Coordinate descent')
[t_nodeOptCoord,outputCoord]=localization_essential_coordinateDescent(t_node,optsCommon{:});

save([mfilename '_data'])

figure(1)
subplot(1,3,1)
testNetworkDisplayErrors(t_node,'n')
title('Init')
subplot(1,3,2)
testNetworkDisplayErrors(t_nodeOptGrad,'n')
title('Gradient descent')
subplot(1,3,3)
testNetworkDisplayErrors(t_nodeOptCoord,'n')
title('Coordinate descent')

figure(2)
semilogy(outputGrad.outputLieMinimize.t,outputGrad.outputLieMinimize.cost,'.-');
hold on
semilogy(outputCoord.tFinal,outputCoord.cFinal,'r*')
hold off
