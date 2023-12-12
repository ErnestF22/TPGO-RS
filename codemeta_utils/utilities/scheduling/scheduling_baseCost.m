function f=scheduling_baseCost(data)
f=diag([8,2,8,1,40])*data.Tron.flag*diag([ones(1,6),[2 2],ones(1,18-6-4),[2 2]]);
f(end,9:18)=2*f(end,9:18);
f(f==0)=999;