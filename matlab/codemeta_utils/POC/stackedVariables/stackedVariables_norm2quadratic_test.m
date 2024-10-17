function stackedVariables_norm2quadratic_test
variables=stackedVariables_set([],{{'x',3},{'theta',1,5},{'c1',1}});

v=randn(stackedVariables_stackSize(variables),1);
varList={{'x'},{'theta',1}};
vSubStack=stackedVariables_substack(variables,varList,v);
for type={'point','model'}
    switch type{1}
        case 'point'
            A=randn(3,size(vSubStack,1));
            b=randn(3,1);
            [Q,p,c]=stackedVariables_norm2quadratic(variables,varList,A,b,type{1});
            disp([norm(A*vSubStack-b)^2 (v'*Q*v+2*p'*v+c)])
        case 'model'
            A=randn(size(vSubStack,1),1);
            b=randn(1,1);
            [Q,p,c]=stackedVariables_norm2quadratic(variables,varList,A,b,type{1});
            disp([norm(vSubStack'*A-b)^2 (v'*Q*v+2*p'*v+c)])
    end
end

           