%function flagGroundTruthPresent=testNetworkCheckGroundTruth(t_node)
%Return true if t_node contains a field called 'gijtruth'.
%If t_node contains a field 'gitruth' but not a field 'gijtruth', automatically
%calls testNetworkAddRelativeGroundTruth and return true.

%%AUTORIGHTS%%

function flagGroundTruthPresent=testNetworkHasGroundTruth(t_node)
if isfield(t_node(1),'gijtruth')
    flagGroundTruthPresent=true;
else
    if isfield(t_node(1),'gitruth')
        t_node=testNetworkAddRelativeGroundTruth(t_node);
        flagGroundTruthPresent=true;
    else
        flagGroundTruthPresent=false;
    end
end
        
