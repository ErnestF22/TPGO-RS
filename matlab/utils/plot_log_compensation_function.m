function plot_log_compensation_function

x = linspace(-10,0.999999999, 10000);

N = length(x);
y = zeros(size(x));

a = 1.0;

for ii = 1:N
    x_ii = x(ii);
    y(ii) = -log(1-a*x_ii) - a*x_ii + x_ii^2;
end

figure(1)
plot(x,y)

end
