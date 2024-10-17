function pos = GetPosInChain(x,k,lens)

pos = x(1:3);
for i=2:k
    pos = pos + lens(i-1)*x(6*i-5:6*i-3);
end



