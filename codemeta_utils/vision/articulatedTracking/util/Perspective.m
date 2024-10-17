function h = Perspective(x)

if numel(x)~=3
    error('Should be a 3-vector')
end

if x(3) <= 0
     error('Depth should be positive')
end


h = [x(1)/x(3) ;...
     x(2)/x(3)] ;

end

