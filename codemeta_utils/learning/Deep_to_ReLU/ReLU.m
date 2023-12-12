% element-wise ReLU function
function y = ReLU(x)
y = zeros(1,length(x));
for i=1:length(x)
    if x(i)>=0
        y(i)=x(i);
    end
end