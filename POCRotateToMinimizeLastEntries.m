function Q_transp =  POCRotateToMinimizeLastEntries(x)

% x = [matStackH(R)];
[Q,~,~]=svd(x*x');
Q_transp = Q';
% disp('x=')
% disp(x)
% disp('Q_transp''*x=')
% disp(Q_transp'*x)

end %file function
