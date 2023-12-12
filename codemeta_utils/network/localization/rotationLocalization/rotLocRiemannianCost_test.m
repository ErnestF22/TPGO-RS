function rotLocRiemannianCost_test
flagSeparateGradient=true;
funs=consensus_rot3_almostGlobal_functions('type','squared');
t_node=testNetworkBuildTestNetwork('N',7);
E=t_node.E;
[R,v]=rot_randGeodFun(t_node.Ritruth);
RRel=t_node.Rijtruth;

funCheckDer(@(t) costAndDer(R(t),v(t)))

    function [c,dc]=costAndDer(R,v)
    if flagSeparateGradient
        c=rotLocRiemannianCost(E,R,RRel,funs);
        gradc=rotLocRiemannianCost_grad(E,R,RRel,funs);
    else
        [c,gradc]=rotLocRiemannianCost(E,R,RRel,funs);
    end        
    dc=sum(rot_metric(R,gradc,v));
    end
end


