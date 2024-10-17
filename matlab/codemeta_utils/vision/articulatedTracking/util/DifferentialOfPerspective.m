function Dh = DifferentialOfPerspective(x)

if numel(x)~=3
    error('Should be a 3-vector')
end

if x(3) <= 0
    error('Depth should be positive')
end


Dh = [1/x(3) 0 -x(1)/x(3)^2;...
      0 1/x(3) -x(2)/x(3)^2] ;
 

 
end

