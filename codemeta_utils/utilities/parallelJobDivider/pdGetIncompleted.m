function idx=pdGetIncompleted(name)
load([name '_pdIndicators'], 'served', 'completed');

idx=find(and(served,~completed));
