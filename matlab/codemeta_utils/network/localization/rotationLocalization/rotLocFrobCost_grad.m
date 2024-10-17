function gradc=rotLocFrobCost_grad(E,R,RRel)
gradc=rotLocBaseCost_grad(E,R,RRel,@rotLocFrobCostPair_grad);
