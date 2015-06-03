function delta=solveNewton(deltanew,params,draws,mkt)

[C,IA,IC] = unique(mkt.prodid);
J = max(IA);

g= @(x)(rc_newton(x));
ops=optimset('Display','off','Jacobian','on','DerivativeCheck','off','FinDiffType','central','maxIter',1000,'TolFun',1e-14);
[delta,~,~,~] = fsolve(g,deltanew,ops);

    function [f,Jac]=rcnl_newton(del)
        [pjt, pijt,~,within] = rc_share_safe(del,params,draws,mkt);
        f=-log(mkt.sjt)+log(pjt);
        if nargout > 1
            Jac=bsxfun(@rdivide,-pijt*diag(draws.w)*pijt',pjt)+speye(J);
        end
    end

end