function optionsExtract_test
optsIn={'type','x','invert','length',5,'flagDummy',true,'RT',eye(3),zeros(2,1),'nameConfuse','RT','dummy'};
optsOut=optionsExtract(optsIn,{'type',1,'invert',0,'RT',2},{'nameConfuse',1});

disp(optsOut)

