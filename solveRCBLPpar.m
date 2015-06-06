function [fval2,thetaval2,betaval2,SE]=solveRCBLPpar(dtable,draws,theta0,extract_fun)
% these can be modified
ops=optimset('Display','iter' ,'Algorithm','Interior-Point','GradObj','on','DerivativeCheck','off');
% This can be 'fixed-point' or 'newton': I recommend fixed-point.
method = 'fixed-point';

% DO NOT modify below this line
% 

% This is the parameter extraction function
get_params = str2func(extract_fun);

% These are bounds on parameters that are set in the parameter extraction function
init_param=get_params(theta0,draws);
lb= init_param.lb;
ub= init_param.ub;

ns = length(draws.v);
n=length(dtable.x1);
X=[dtable.price dtable.x1];
Z=[dtable.x1 dtable.z];
K=length(theta0);

f = @(x)evalSingle(x);
theta0 = [ones(K,1)];
W=inv(Z'*Z);

% This bit just checks that the derivatives match the numeric ones
myDerivCheck(f,theta0,1e-5)

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
        if exist('knitromatlab'),
            [that]=knitromatlab(f,x0,[],[],[],[],lb,ub,[],[],ops);
        else,
            [that]=fmincon(f,x0,[],[],[],[],lb,ub,[],ops);
        end
        
        % After optimization recover the linear parameters and objective 
        thetahat =get_params(that,draws);
        delta=solveAllShares(tableA,draws,thetahat,method);
        dtable.delta=delta;
        [beta,resid]=ivregression(delta,X,Z,W);
        fval=(resid'*Z)*W*(resid'*Z)';
    end
    
    % evaluate the GMM objective and its gradient once
    function [fval,g]=evalSingle(theta)
        % this extracts the parameters for your specification
        p = get_params(theta,draws);
        [delta,Jac]=solveAllShares(dtable,draws,p,method);
        dtable.delta=delta;
        [beta,resid]=ivregression(delta,X,Z,W);
        fval=(resid'*Z)*W*(resid'*Z)';
        g=-2*(Jac'*Z)*W*(resid'*Z)';
        %disp(['Price Coeff: ' num2str(beta(1)) ' and f: ' num2str(fval)])
    end

    function [S,SE]=getCovariance(theta,dtable,draws)
        % this extracts the parameters for your specification
        params1=get_params(theta,draws);
        [dstar,Jac]=solveAllShares(dtable,draws,params1,method);
        [bhat,uhat]=ivregression(dstar,X,Z,W);

        g=bsxfun(@times,uhat,Z);
        gstar=bsxfun(@minus,g,mean(g));
        S=inv(gstar'*gstar);
        G=Z'*[Jac X];

        SE=full(sqrt(diag(inv(G'*S*G))));
    end
end
