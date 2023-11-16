A = [[4 0 1]; [2 3 2]; [1 0 4]];
fun_han = @(x) A*x;
v_start = rand(d,1);
stief_normalize_han = @(x) x / sqrt(x' * x);

[lambda_max, v_max] = pim_function(fun_han, v_start, stief_normalize_han, thresh);

disp("v_max")
disp(v_max)


disp("lambda_max")
disp(lambda_max)