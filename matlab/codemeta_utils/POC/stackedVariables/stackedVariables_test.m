function stackedVariables_test
variables=stackedVariables_set([],{'x',3});
variables=stackedVariables_set(variables,{'theta',1,5});
variables=stackedVariables_set(variables,{'c1',1});
stackedVariables_disp(variables)

disp('Positions for theta')
disp(stackedVariables_positionsGet(variables,'theta'))
disp('Positions for x and theta(2:3)')
disp(stackedVariables_positionsGet(variables,{{'x'},{'theta',2:3}}))
disp('Their stacked dimension')
disp(stackedVariables_stackSize(variables,{{'x'},{'theta',2:3}}))

v=randn(stackedVariables_stackSize(variables),1);
vStruct=stackedVariables_unstack(variables,v);
v2=stackedVariables_stack(variables,vStruct);

assert(~any(v-v2))
