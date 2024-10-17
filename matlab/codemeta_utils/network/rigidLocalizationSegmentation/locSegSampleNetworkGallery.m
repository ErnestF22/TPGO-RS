function [E,x,R,constraintList,constraintParametersList]=locSegSampleNetworkGallery(sampleNum,varargin)
constraintSampleNum=3;
dimData=3;

%optional parameters
ivarargin=1;
while(ivarargin<=length(varargin))
    switch(lower(varargin{ivarargin}))
        case 'constraintnum'
            ivarargin=ivarargin+1;
            constraintSampleNum=varargin{ivarargin};
        case 'dimdata'
            ivarargin=ivarargin+1;
            dimData=varargin{ivarargin};
        otherwise
            disp(varargin{ivarargin})
            error('Argument not valid!')
    end
    ivarargin=ivarargin+1;
end

fprintf('Graph from gallery: %d\n',sampleNum)
fprintf('Data dimension: %d\n',dimData)

switch sampleNum
    case 1 % two separate components (one line, one triangle)
        E=[
            1,2;
            3,4;
            3,5;
            4,5
            ];
    case 2 % two triangles joined at node 3
        E=[
            1,2;
            1,3;
            2,3;
            3,4;
            3,5;
            4,5
            ];
    case 3 % one 2D rigid component
        E=[
            1,2;
            1,3;
            2,3;
            3,4;
            3,5;
            4,5;
            1,4
            ];
    case 4 % two 3D rigid component
        E=[
            1,2;
            1,3;
            2,3;
            3,4;
            3,5;
            4,5;
            1,5;
            2,5
            ];
    case 5 %loop with 4 elements
        E=[
            1,2;
            2,3;
            3,4;
            4,1
            ];
    case 6 %chain with 2 elements
        E=[
            1,2
            ];
       
    case 7 %chain with 3 elements
        E=[
            1,2;
            2,3;
            ];
    case 8 %loop with 3 elements
        E=[
            1,2;
            2,3;
            3,1
            ];
    case 9 % three triangles joined at node 3
        E=[
            1,2;
            1,3;
            2,3;
            3,4;
            3,5;
            4,5;
            3,6;
            6,7;
            7,3
            ];
    case 10 % full-connected with four nodes
        E=[
            1,2;
            2,3;
            3,4;
            4,1;
            2,4;
            3,1
            ];
    case 11 % two 3D rigid component
        E=[
            1,2;
            1,3;
            2,3;
            3,4;
            3,5;
            4,5;
            1,5;
            2,5;
            4,2
            ];
    case 12
        E=[
            1,2;
            2,3;
            3,4;
            4,1;
            2,4
            ];
    case 13 % a central triangle (134) plus three others connected to central one through a different edge
        E=[
            1,2;
            1,3;
            2,3;
            3,4;
            3,5;
            4,5;
            1,4;
            1,6;
            4,6
            ];
    case 14 % three triangles: one joined only through 3, the other two joined through 3 and 5
        E=[
            1,2;
            1,3;
            2,3;
            3,4;
            3,5;
            3,6;
            4,5;
            5,6
            ];
    case 15 % two squares connected through an edge
        E=[
            1,2;
            1,3;
            2,4;
            3,4;
            3,5;
            4,6;
            5,6
            ];
    case 16 % two triangles connected through the edge (2,3)
        E=[
            1,2;
            1,3;
            2,3;
            3,4;
            2,4
            ];
%         x=[
%             1,0,0;
%             0,0,0;
%             0,1,0;
%             0,0,1
%             ]';
        R=repmat(eye(3),[1 1 4]);
    case 17 % a central triangle (134) plus one other and two thethraedrons connected to central one through a different edge
        E=[
            1,2;
            1,3;
            2,3;
            3,4;
            3,5;
            4,5;
            1,4;
            1,6;
            4,6;
            7,1;
            7,2;
            7,3;
            8,3;
            8,4;
            8,5
            ];
    case 18 % two triangles joined at node 3, addition of relative rotation constraints
        E=[
            1,2;
            1,3;
            2,3;
            3,4;
            3,5;
            4,5
            ];
    case 19 % chain of three triangles joined at node 3 and node 5
        E=[
            1,2;
            1,3;
            2,3;
            3,4;
            3,5;
            4,5;
            5,6;
            6,7;
            7,5
            ];
        x=[
            0,0,0;
            0,1,0;
            0,0,1;
            1,1,1;
            2,0,1;
            2,0,0;
            2,1,0
            ]';
        R=repmat(eye(3),[1 1 7]);
    case 20 %four triangles: two joined by an edge, the other joined by the two endpoints of that edge
        E=[
            1,2;
            2,3;
            3,1;
            3,4;
            3,5;
            3,6;
            4,6;
            5,6;
            6,7;
            6,8;
            7,8
            ];
    case 21 %same as 20, but one of the central triangles is a tethraedron (3569)
        E=[
            1,2;
            2,3;
            3,1;
            3,4;
            3,5;
            3,6;
            4,6;
            5,6;
            6,7;
            6,8;
            7,8;
            3,9;
            5,9;
            6,9
            ];
    case 22 %two disconnected components, relative translations constraints
        E=[
            1,2;
            2,3;
            4,5;
            5,6;
            ];
    case 23 %six triangles all connected through node 1
        E=[
            1,2;
            1,3;
            2,3;
            ...
            1,4;
            1,5;
            4,5;
            ...
            1,6;
            1,7;
            6,7;
            ...
            1,8;
            1,9;
            8,9;
            ...
            1,10;
            1,11;
            10,11;
            ...
            1,12;
            1,13;
            12,13;
            ];
    case 24 %star with 6 nodes connected to node 1 
        E=[
            1,2;
            1,3;
            1,4;
            1,5;
            1,6
            ];
    case 25 %two disconnected edges, relative translations constraints
        E=[
            1,2;
            3,4;
            ];
    case 26 %two disconnected components, one chain and one triangle, similar to 22
        E=[
            1,2;
            2,3;
            4,5;
            5,6;
            4,6
            ];
    case 27 %two disconnected triangles, similar to 22
        E=[
            1,2;
            2,3;
            1,3;
            4,5;
            5,6;
            4,6
            ];
    case 28 %two triangles (123,456) joined  by an edge (34)
        E=[
            1,2;
            2,3;
            1,3;
            3,4;
            4,5;
            5,6;
            4,6
            ];
    case 29 %two squares with one edge in common (34), the first one has the diagonal (24)
        E=[
            1,2;
            2,3;
            3,4;
            4,1;
            2,4;
            3,5;
            5,6;
            6,4
            ];
        x=[
            0.1 1   1.1 0   1   0;
            0   0   1   1   2.1 2;
            0   0   0   0   0   0
            ];
        
    case 30 %illustration 1
        E=[
            1,2;
            2,3;
            3,1;
            3,4;
            4,1;
            4,5;
            5,6;
            7,8;
            8,9;
            9,7;
            9,10;
            8,10;
            10,11;
            11,12;
            12,8
            ];
        x=[
   -2.4309   -3.7212   -2.4309   -1.2788    0.2419    0.3111    0.9562    0.9332    2.6382    2.7074    4.5046    2.5000;
   -1.0673   -2.1784   -3.4064   -1.8860   -0.8041   -3.4942    0.0146    2.3246    0.8041    2.6754    3.9912    4.1959;
    0         0         0         0         0         0         0         0         0         0         0         0     
        ];
    case 31 %illustration 1
        E=[
            1,2;
            2,3;
            3,1;
            3,4;
            4,1;
            4,5;
            5,6;
            7,8;
            8,9;
            9,7;
            9,10;
            8,10;
            10,11;
            11,12;
            12,8;
            12,10;
            12,9
            ];
        x=[
   -2.4309   -3.7212   -2.4309   -1.2788    0.2419    0.3111    0.9562    0.9332    2.6382    2.7074    4.5046    2.5000;
   -1.0673   -2.1784   -3.4064   -1.8860   -0.8041   -3.4942    0.0146    2.3246    0.8041    2.6754    3.9912    4.1959;
    0         0         0.5       0         0.5      -0.5      -0.5       0         0         0.5       1         1.5     
        ];            
    
    case 32
        x=[
   -2.4309   -3.7212   -2.4309   -1.2788    0.2419    0.3111    0.9562    0.9332    2.6382    2.7074    4.5046    2.5000   -2.1544    1.7857    2.3848    0.5415
   -1.0673   -2.1784   -3.4064   -1.8860   -0.8041   -3.4942    0.0146    2.3246    0.8041    2.6754    3.9912    4.1959   -1.7105    2.7632   -2.9386   -4.0789
          ];
    case 40 % one square
        E=[
            1,2;
            2,3;
            3,4;
            4,1;
            ];
        x=[
            0.1 1   1.1 0;
            0.1 0   1   1;
            0   0   0   0;
            ];
    case 41 % two triangles joined at an edge plus two edges in chain
        E=[
            1,2;
            2,3;
            3,1;
            3,4;
            4,1;
            4,5;
            5,6;
            ];
        x=[
   -2.4309   -3.7212   -2.4309   -1.2788    0.2419    0.3111;
   -1.0673   -2.1784   -3.4064   -1.8860   -0.8041   -3.4942;
    0         0         0.5       0         0.5      -0.5      
        ];            
        
    case 42 %two disconnected components, one v and one edge
        E=[
            1,2;
            3,4;
            4,5;
            ];
    otherwise
        error('sampleNum not valid')
end
E=[E;fliplr(E)];
N=max(E(:));

A=zeros(N);
A(sub2ind([N N],E(:,1),E(:,2)))=1;

if ~exist('x','var')
    x=randn(dimData,N);
    x=x-mean(x,2)*ones(1,N);
    [U,S,V]=svd(x,'econ');
    x=graphDrawing(A,dimData)+0.2*randn(size(x));
else
    if dimData==2
        x=x(1:2,:);
    end
end
if ~exist('R','var')
    R=rot_randn(eye(dimData),[],N);
end

switch sampleNum
    case 18
        [constraintList,constraintParametersList,E]=locSegSampleConstraintsGallery(3,E,x,R);
        [constraintList,constraintParametersList,E]=locSegSampleConstraintsGallery(4,E([1:3 7:9],:),x,R,constraintList,constraintParametersList);
%     case {22,25}
%         [constraintList,constraintParametersList,E]=locSegSampleConstraintsGallery(5,E,x,R);
end

if ~exist('constraintList','var')
    fprintf('Constraint from gallery: %d\n',constraintSampleNum);
    [constraintList,constraintParametersList,E]=locSegSampleConstraintsGallery(constraintSampleNum,E,x,R);
end

