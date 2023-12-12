clear 
close all
addpath('old_functions/')
t1 = 250;
t2 = 250;
for d = 2:6
    for j = 1:t1
        nConsAdded = randi([d,d+7]); % number of constraints
        limL = -15; limU = 15;
        Acons = limL + (limU-limL) * rand([nConsAdded,d]);
        Acons(:,end+1) = limU*rand([nConsAdded,1]);
        % A = [Acube;Acons];
        A = Acons;
        % load counterexample case(if stored)
        % load('At.mat');
        % load('Atset.mat');
        % load('A_v21.mat');
        [As,basic] = prep_table(A,d);
        [basic,result,Ar,~]= dual_simplex_xpn_old(As,basic,2*d);
        
        vertices = [];
        acts = [];
        if ~isempty(basic)
            [acts, vertices] = find_vertices_query_xpn_old(basic,Ar,2*d);
%             vertices = unique(round(vertices,4),'rows');
            acts = unique(acts,'rows');
        end
        
        c = 0;
        UnsuccessfulCase = [];
        Atset = {}; % store counterexample case
        for k = 1:t2
            AT = A;
            AT = transfer_PDM(AT,d);
            Atset{k} = AT; % store counterexample case
            %         AT = Atset{k}; % load counterexample case
            
            [As,basic] = prep_table(AT,d);
            [basic,Ar,~]= dual_simplex_xpn(As,basic);
            verticesT = [];
            actT = [];
            if ~isempty(basic)
                [actT, verticesT] = find_vertices_query_xpn(basic,Ar);
%                 uniqueverticesT = unique(round(verticesT,4),'rows');
                actT = unique(actT,'rows');
            end
            % if size(uniqueverticesT,1)==size(vertices,1)
            if size(actT,1)==size(acts,1)
                c=c+1;
            else
                sizeactT = size(actT,1);
                UnsuccessfulCase = [UnsuccessfulCase;k,sizeactT,isempty(basic)]; % UnsuccessfulCase = [k, number of vertices find, if the basic is empty from dual_simplex]
            end
        end
        if ~isempty(UnsuccessfulCase)
            disp(['j=',num2str(j), ', Unsuccess!']);
            save(['FailMat',num2str(j),'.mat'],'A','Atset');
            save(['Case',num2str(j),'.mat'],'UnsuccessfulCase');
%         else
%             disp(['j=',num2str(j), ', Success!']);
%             continue
        end
    end
    if c == t2
         disp(['d=',num2str(d), ', Success!']);
    else
        disp(['d=',num2str(d), ', Unsuccess!']);
    end
    disp(['d=',num2str(d), ', Completed!']);
end