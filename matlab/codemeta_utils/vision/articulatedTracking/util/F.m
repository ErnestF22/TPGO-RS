function [F1,F2] = F(x,v)



 F1 =   ExpMapSphere(x,v);
 F2 =  ParallelTranspot(x,F1,v); 


end

