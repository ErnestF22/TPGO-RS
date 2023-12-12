% cd '/Library/gurobi912/mac64/matlab'
% gurobi_setup
addpath ('/Library/gurobi912/mac64/matlab')
gurobi_setup
network = ['nnet/ACASXU_run2a_' num2str(i) '_' num2str(j) '_batch_2000.nnet'];
Ab_set = load_network(network);
Ab_set(2:6) = [];
input = [1500;-0.06;3.1;980;960];
layer = size(Ab_set,2);
z_size = 0;
for i=1:layer
    z_size = z_size+size(Ab_set(i).b,1);
end
result = mip(z_size,input,Ab_set);

function result = mip(z_size,input,Ab_set)
model = {};
input_size = size(input,1);
layer = size(Ab_set,2);
x = model.addVars(input_size, vtype='C'); % variables [y1,y2,...,yn]
y = model.addVars(z_size, vtype=GRB.CONTINUOUS); % variables [y1,y2,...,yn]
b_on = model.addVars(z_size, vtype=GRB.BINARY); % Relu activate
b_off = model.addVars(z_size, vtype=GRB.BINARY); % Relu inactivate
model.setObjective(norm(x-input), GRB.MINIMIZE); % objective function

start = 1;
xi = x;
for i=1:layer
    A = Ab_set(i).A;
    b = Ab_set(i).b;
    y_size = size(b,1);
    yi = y(start:start+y_size-1);
    b_offi = b_off(start:start+y_size-1);
    b_oni = b_on(start:start+y_size-1);
    model.addConstr(yi >= 0);
    model.addConstr(A@xi + b - yi - M*b_offi <= 0);
    model.addConstr(A@xi + b - yi + M*b_offi >= 0);
    model.addConstr(yi - M*b_oni <= 0);
    model.addConstr(xi - M*b_oni <= 0);
    xi = y(start:start+y_size-1);
    start = start+y_size;
end


gurobi_write(model, 'mip.lp');

result = gurobi(model, params);

end