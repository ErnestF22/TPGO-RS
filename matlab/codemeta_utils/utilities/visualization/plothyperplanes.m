function h=plothyperplanes(b,varargin)

[K,n]=size(b);
styles=char('b','r','g','c','m','y','k');

flaghold=ishold();

for(iplane=1:n)
    h(iplane,1)=draw3dPlane(b(:,iplane),'style',styles(iplane,:),varargin{:});
    hold on
end

if(~flaghold)
    hold off
end
