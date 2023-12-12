function POCsubwindowEntropyMax
n=10; %# of pixel levels
x=randi(20,n,1); %random distribution of pixels
cnt=sum(x); %total pixel count
cntSW=ceil(0.75*cnt); %subwindow pixel count
%optimize entropy subject to constraints
cvx_begin
    variable xSW(n)
    maximize(sum(entr(xSW/cntSW)))
    subject to 
        sum(xSW)==cntSW
        xSW<=x
        %xSW>=0     %sometimes this makes CVX fail
cvx_end
%display results
disp([cnt sum(x) cntSW sum(xSW)])
disp([x';xSW'])
plot(1:n,x,1:n,xSW)

