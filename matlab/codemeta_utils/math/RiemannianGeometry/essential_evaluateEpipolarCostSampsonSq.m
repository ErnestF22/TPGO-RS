function [c,gradcQ,hessOpcQ]=essential_evaluateEpipolarCostSampsonSq(Q,x,varargin)
    flagSymmetricHess=false;

    %optional parameters
    ivarargin=1;
    while(ivarargin<=length(varargin))
        switch(lower(varargin{ivarargin}))
            case 'symmetrichess'
                flagSymmetricHess=true;
            case 'flagsymmetrichess'
                ivarargin=ivarargin+1;
                flagSymmetricHess=varargin{ivarargin};
        end
        ivarargin=ivarargin+1;
    end

    flagComputeGrad=nargout>1;
    flagComputeHess=nargout>2;
    

    E=essential_toE(Q);
    asym=@(A) (A-A')/2;
    if ~flagComputeGrad
        c=epipolarCostFromE_SampsonSq(E,x(:,:,1),x(:,:,2));
    else
        if ~flagComputeHess
            [c,gradc]=epipolarCostFromE_SampsonSq(E,x(:,:,1),x(:,:,2));
        else
            [c,gradc,hessOpc]=epipolarCostFromE_SampsonSq(E,x(:,:,1),x(:,:,2));
        end    
        gradcQ=-[asym(gradc*E'); asym(gradc'*E)];
    end
    
    function h=hOp(hatv0)
        hatv01=essential_getR1(hatv0);
        hatv02=essential_getR2(hatv0);
        symm=@(A) (A+A')/2;
        dE= hatv01'*E+E*hatv02;
        h=-[asym(hessOpc(dE)*E'); asym(hessOpc(dE)'*E)];
        if ~flagSymmetricHess
            h=h-[asym(gradc*dE'); asym(gradc'*dE)];
        else
            h=h-[asym(symm(gradc*E')*hatv01+gradc*hatv02'*E');asym(gradc'*hatv01'*E+symm(gradc'*E)*hatv02)];
        end            
    end

    if flagComputeHess
        hessOpcQ=@hOp;
    end
end
