function admmMed_auxiliaryCost_test
muijk=[0 1; 0 1];
lambdaijk=-[-1;1];
rho=1;

f=@(x) admmMed_auxiliaryCost(muijk,lambdaijk,rho,'zEval',x);

funImage(linspace(-1,2),linspace(-1,2),f,'method','surf');
[z,c]=admmMed_auxiliary(muijk,lambdaijk,rho);
hold on
plotPoints([z;c],'o')
hold off
view(3)