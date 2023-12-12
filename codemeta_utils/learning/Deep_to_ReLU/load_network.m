function Ab_set = load_network(network)
% Load txt file
Wb = fileread(network);
% The first line gives the four values:
%     Number of hidden layers in the networks
%     Number of inputs to the networks
%     Number of outputs from the networks
%     Maximum size of any hidden layer
L1 = textscan(Wb,'%d %d %d %d', 1, 'HeaderLines', 3, 'Delimiter',',');
[L,di,do,dm] = L1{:};
% The second line gives the sizes of each layer, including the input and output layers
L2 = cell2mat(textscan(Wb,repmat('%d',[1,L+1]), 1, 'HeaderLines', 4, 'Delimiter',','));
% The eighth line gives weights and biases
S = regexprep(Wb, '[-+]\d\d\d;', 'e$&');
values = cell2mat(textscan(S, '%.10f%.10f%.10f', 'Delimiter',',\t','HeaderLines', 10,'CollectOutput',1));
Wbmat = reshape(values',1,[]);
Wbmat(find(isnan(Wbmat)))=[];
% Assign weight and bias for each layer(in Ab_set), e.g. Ab_set(1).A = A of layer 1
for i=1:L
    dim = L2(i)*L2(i+1);
    Ab_set(i).A = reshape(Wbmat(1:dim),L2(i),L2(i+1))';
    Wbmat(1:dim) = [];
    Ab_set(i).b = reshape(Wbmat(1:L2(i+1)),L2(i+1),1);
    Wbmat(1:L2(i+1)) = [];
end
end