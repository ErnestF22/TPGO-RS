%Test numerical accuracy of log_se3 combined with exp_se3

%%AUTORIGHTS%%
function log_se3_test
for(it=1:1000)
    v=10*randn(3,1);
    w=randn(3,1);
    w=w/norm(w)*pi*rand;
    
    [w1,v1]=log_se3(exp_se3(w,v));
    err(it)=norm([w;v]-[w1;v1]);
    
    if(err(it)>1e-6)
        error('Check the code!')
    end
    
end
plot(err)
title('All errors should be small')
xlabel('Trial')
ylabel('Norm of (w,v)-log(exp(w,v))')

