function h_input = get_h_input(x,v,A,act_set)
% calculate distance from input point

% h_input = sum(abs(v-x')); % L-1 norm
h_input = sum((v-x').^2); % L-2 norm
% h_input = max(abs(v-x')); % infinity norm

% load('A_initial.mat','A_initial');
% d = size(A_initial,2)-1;
% h_input_set = [];
% for act = act_set
%     act = act-d*2;
%     A1 = A_initial(act,:);
%     b1 = A1(end);
%     A1(end) = [];
%     h_input_set = [h_input_set; abs(A1*x-b1)/norm(A1)];
% end
% h_input = min(h_input_set);