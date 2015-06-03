function [Jac,f]=RCBLP_Jacobian(params,draws,mkt,method)
[C,IA,IC] = unique(mkt.prodid);
J = max(IA);

nK = size(mkt.x2,2);
ns = size(draws.v,1);

[pjt, pijt] = rc_share_safe(mkt.delta,params,draws,mkt);

%one shot for deltas -- this is NOT normalized by 1/pjt
Jac_d=-pijt*diag(draws.w)*pijt'+diag(pjt);

% this is ok
Jac_theta = zeros(J,nK);
for kk=1:size(params.betamask,1),
    Jac_theta(:,kk)= JacSigma(params.betamask(kk,1), params.betamask(kk,2));
end

if(strcmp(method,'sterr')),
    % rescale by 1./pjt to prevent underflow
    Jac = (diag(1./pjt)*Jac_d)\(diag(1./pjt)*Jac_theta);
elseif(strcmp(method,'mpec')),
    Jac = bsxfun(@rdivide,[Jac_theta Jac_d],pjt);
end
f=-log(mkt.sjt)+log(pjt);

    function [JacVec]=JacSigma(xindex,vindex)
        vk=params.dbeta(:,vindex);
        xk=mkt.x2(:,xindex);
        JacVec=(bsxfun(@times,pijt,vk').*bsxfun(@minus,xk,xk'*pijt))*draws.w;
    end
end
