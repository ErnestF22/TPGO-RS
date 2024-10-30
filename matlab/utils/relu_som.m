function x_out = relu_som(x_in)
%RELU_SOM Apply ReLU to x_in array
% This function is needed as the native MATLAB one only works with dlarrays

x_out = zeros(size(x_in));
x_out(x_in > 0) = x_in(x_in > 0);

end

