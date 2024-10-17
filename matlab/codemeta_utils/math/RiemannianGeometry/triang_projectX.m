function x=triang_project(X)
x=[projectFromRT(eye(3),zeros(3,1),X) projectFromRT(eye(3),[0;0;-1],X,'references')];
