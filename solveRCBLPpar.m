function [fval2,thetaval2,betaval2,SE]=solveRCBLPpar(dtable,draws,theta0,bounds)
% these can be modified
ops=optimset('Display','iter' ,'Algorithm','Interior-Point','GradObj','on','DerivativeCheck','off');
% This can be 'fixed-point' or 'newton': I recommend fixed-point.
method = 'fixed-point';

% These are bounds on parameters that are passed to the function
lb= bounds.lb;
ub= bounds.ub;

ns = length(draws.v);
n=length(dtable.x1);
X=[dtable.price dtable.x1];
Z=[dtable.x1 dtable.z];
K=length(theta0);

f = @(x)evalSingle(x);
theta0 = [ones(K,1)];
W=(Z'*Z)\eye(size(Z,2));

% This bit just checks that the derivatives match the numeric ones
deps = 1e-5;
[gstar]=myDerivCheck(f,theta0)
[fval,gval] = f(theta0);
[gstar gval]

tic
% first step
[fval,thetaval,betaval]=get_results(dtable,theta0);

% Lazy way to display key parameters (fix this later)
thetaval
betaval(1)
save first-step.mat

% update weight matrix and produce second-step
[W]=getCovariance(thetaval,dtable,draws);
[fval2,thetaval2,betaval2]=get_results(dtable,theta0);
[S,SE]=getCovariance(thetaval2,dtable,draws);

    %For a fixed weighting matrix (W) this minimizes the GMM objective
    %this works much better if knitromatlab is installed
    function [fval,that,beta]=get_results(tableA,x0)
        % function handle f is mapped to evalsingle below for a (X,Z,W) 
        [that]=fmincon(f,x0,[],[],[],[],lb,ub,[],ops);
        %[that]=knitromatlab(f,x0,[],[],[],[],lb,ub,[],[],ops);
        
        % After optimization recover the linear parameters and objective 
        thetahat =extract_params(that,draws);
        delta=solveAllShares(tableA,draws,thetahat,method);
        dtable.delta=delta;
        [beta,resid]=ivregression(delta,X,Z,W);
        fval=(resid'*Z)*W*(resid'*Z)';
    end
    
    % evaluate the GMM objective and its gradient once
    function [fval,g]=evalSingle(theta)
        % this extracts the parameters for your specification
        p = extract_params(theta,draws);
        [delta,Jac]=solveAllShares(dtable,draws,p,method);
        dtable.delta=delta;
        [beta,resid]=ivregression(delta,X,Z,W);
        fval=(resid'*Z)*W*(resid'*Z)';
        g=-2*(Jac'*Z)*W*(resid'*Z)';
        %disp(['Price Coeff: ' num2str(beta(1)) ' and f: ' num2str(fval)])
    end

    function [S,SE]=getCovariance(theta,dtable,draws)
        % this extracts the parameters for your specification
        params1=extract_params(theta,draws);
        [dstar,Jac]=solveAllShares(dtable,draws,params1,method);
        [bhat,uhat]=ivregression(dstar,X,Z,W);
        g=bsxfun(@times,uhat,Z);
        S=inv(g'*g);
        
        G=Z'*[Jac X];
        bread=inv(G'*W*G);
        sand=(G'*W*S*W*G);
        
        % sqrt n comes from lack of 1/n term in S and W
        SE=full(sqrt(diag(bread'*sand*bread)./n));
    end

    function [gradstar]=myDerivCheck(func,theta0)
        gradstar = zeros(size(theta0));
        for k=1:length(theta0)
           thetastarA = theta0;
           thetastarB = theta0;

           thetastarA(k) = thetastarA(k) + deps;
           thetastarB(k) = thetastarB(k) - deps;
           gradstar(k) = (func(thetastarA) - func(thetastarB))./(2*deps);
        end
    end
end
