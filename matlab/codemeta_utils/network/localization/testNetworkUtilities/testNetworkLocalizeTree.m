function t_node=testNetworkLocalizeTree(t_node,ETree,varargin)

methodAbsolutePoses='references';

%optional parameters
ivarargin=1;
while(ivarargin<=length(varargin))
    switch(lower(varargin{ivarargin}))
        case 'references'
            methodAbsolutePoses='reference';
        case 'poses'
            methodAbsolutePoses='pose';
        case 'methodabsoluteposes'
            ivarargin=ivarargin+1;
            methodAbsolutePoses=lower(varargin{ivarargin});
        otherwise
                error(['Argument ' varargin{ivarargin} ' not valid!'])
    end
    ivarargin=ivarargin+1;
end

structType=testNetworkDetectStructType(t_node);


switch structType
    case 'single'
        t_node.gi(:,:,ETree(1,1))=eye(4);
    case 'array'
        t_node(ETree(1,1)).gi=eye(4);
end

for iEdgeTree=1:size(ETree,1)
    e1=ETree(iEdgeTree,1);
    e2=ETree(iEdgeTree,2);

    lij=1;  %default value if lambdaij members are not present
    switch structType
        case 'single'
            iEdge=find(t_node.E(:,1)==e1 & t_node.E(:,2)==e2,1,'first');
            if isempty(iEdge)
                %try to look for reversed edge
                iEdge=find(t_node.E(:,1)==e2 & t_node.E(:,2)==e1,1,'first');
                if isempty(iEdge)
                    error('Invalid edge in the tree!')
                end
            end
            gij=t_node.gij(:,:,iEdge);
            gi=t_node.gi(:,:,e1);
            if isfield(t_node,'lambdaij')
                lij=t_node.lambdaij(iEdge);
            end
        case 'array'
            gij=t_node(e1).gij(:,:,e2);
            gi=t_node(e1).gi;
            if isfield(t_node(e1),'lambdaij')
                lij=t_node(e1).lambdaij(e2);
            end
    end
    
    gij(1:3,4)=lij*gij(1:3,4);
    switch methodAbsolutePoses
        case 'reference'
            gj=gi*gij;
        case 'pose'
            gj=invg(gij)*gi;
    end

    switch structType
        case 'single'
            t_node.gi(:,:,e2)=gj;
        case 'array'
            t_node(e2).gi=gj;
    end
    
end