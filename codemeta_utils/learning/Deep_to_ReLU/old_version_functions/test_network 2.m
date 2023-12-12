function [num_region, size_ACT_inregion] = test_network(network,x)
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

%% find the neighborhood region
[y,z] = forward(x,L,Ab_set);
[A,b] = get_cumulative(Ab_set,L,z);
num_region = 0;
size_ACT_inregion = 0;
if sum(sum([A,b],2)~= 0)% check if A=0 or b=0
    [As,basic,d] = getdualAmatrix(z,Ab_set);
    d = 2*d;
    [basic,~,Ar,P]=dual_simplex_xpn(As,basic,d);
    if ~isempty(basic)% if there exists basic
%         disp('Region found');
        [ACT, V]=find_vertices_query_xpn(basic,Ar,d);
        num_region = num_region+1;
        size_ACT_inregion = size(ACT,1);
        Z(1).z = z;
        Z(1).index_set = unique(ACT)';
        z_check_list = z;
    end
else
%     disp('Initial region infeasible');
    Z = [];
end
% while ~isempty(Z)
%     z = Z(1).z;
%     index_set = Z(1).index_set;
%     Z(1) = [];
%     for idx = index_set
%         if size(Z,2)>=20
%             break
%         end
%         idx = idx - d;
%         z_new = z;
%         z_new(idx,:)=~z_new(idx,:);
%         if sum(sum(abs(z_check_list-z_new),1)~=0)==0 % if checked
% %         if ismember(z_new',z_check_list','rows') % if checked
%             continue
%         else
%             z_check_list = [z_check_list,z_new];
%             [A,b] = get_cumulative(Ab_set,L,z_new);
%             if ~sum(sum([A,b],2)~= 0)% check if A=0 or b=0
%                 continue
%             end
%             [As,basic,d] = getdualAmatrix(z_new,Ab_set);
%             d = 2*d;
%             [basic,~,Ar,P]=dual_simplex_xpn(As,basic,d);
%             if isempty(basic)% if there exists basic
%                 continue
%             else
%                 [ACT, V]=find_vertices_query_xpn(basic,Ar,d);
%                 Z(end+1).z = z_new;
%                 Z(end).index_set = unique(ACT)';
%                 num_region = num_region+1;
%                 size_ACT_inregion = [size_ACT_inregion,size(ACT,1)];
%             end
%         end
%     end
%     if size(Z,2)>=20
%         break
%     end
% end
end