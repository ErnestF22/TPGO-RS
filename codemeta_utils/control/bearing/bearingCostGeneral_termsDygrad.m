function Hi=bearingCostGeneral_termsDygrad(Yei,Ygi,dfci,ddfci)

I=eye(size(Yei,1));
PYei=I-Yei*Yei';
Hi=-(dfci*I+ddfci*PYei*Ygi*Yei');