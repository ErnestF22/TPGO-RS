function bearingCluster_tikz(x,E,varargin)
NEdges=size(E,1);
membership=ones(1,NEdges);
fid=1;
flagNormalizeScale=true;
flagLabel=false;
formatNameNode='v';
formatNode='vertex';
formatEdge='group';
scale=5;
flagOpenFile=false;
d=size(x,1);

%optional parameters
ivarargin=1;
while(ivarargin<=length(varargin))
    switch(lower(varargin{ivarargin}))
        case 'membership'
            ivarargin=ivarargin+1;
            membership=varargin{ivarargin};
        case 'fileid'
            ivarargin=ivarargin+1;
            fid=varargin{ivarargin};
        case 'flagnormalizescale'
            ivarargin=ivarargin+1;
            flagNormalizeScale=varargin{ivarargin};
        case 'flaglabel'
            ivarargin=ivarargin+1;
            flagLabel=varargin{ivarargin};
        case 'formatnamenode'
            ivarargin=ivarargin+1;
            formatNameNode=varargin{ivarargin};
        case 'formatnode'
            ivarargin=ivarargin+1;
            formatNode=varargin{ivarargin};
        case 'formatedge'
            ivarargin=ivarargin+1;
            formatEdge=varargin{ivarargin};
        case 'scale'
            ivarargin=ivarargin+1;
            scale=varargin{ivarargin};
        otherwise
            error(['Argument ' varargin{ivarargin} ' not valid!'])
    end
    ivarargin=ivarargin+1;
end

if ischar(fid)
    fileName=fid;
    fid=fopen(fileName,'wt');
    flagOpenFile=true;
end

if flagNormalizeScale
    x(1,:)=x(1,:)-min(x(1,:));
    x(2,:)=x(2,:)-min(x(2,:));
    x=x/max(x(:))*scale;
end

NPoints=size(x,2);
for iPoint=1:NPoints
    if flagLabel
        labelString=[',label=45:' num2str(iPoint)];
    else
        labelString='';
    end
    formatCoord=['\\coordinate[' formatNode '%s] (' formatNameNode '%d) at '];
    switch d
        case 2
            fprintf(fid,[formatCoord '(%.2f,%.2f);\n'],labelString,iPoint,x(:,iPoint));
        case 3
            fprintf(fid,[formatCoord '(%.2f,%.2f,%.2f);\n'],labelString,iPoint,x(:,iPoint));
    end
end

c=unique(membership);
for ic=1:length(c)
    idxGroup=find(membership==c(ic));
    for iEdge=1:length(idxGroup)
        fprintf(fid,['\\draw[' formatEdge '%d] '],ic);
        fprintf(fid,['(' formatNameNode '%d) -- (' formatNameNode '%d) '],E(idxGroup(iEdge),:));
        fprintf(fid,';\n');
    end
end

if flagOpenFile
    fclose(fid);
    disp(['Output written to ' fileName])
end
