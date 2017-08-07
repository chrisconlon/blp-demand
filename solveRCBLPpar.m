function [results2]=solveRCBLPpar(dtable,draws,theta0,extract_fun)
% these can be modified
ops=optimset('Display','iter' ,'Algorithm','Interior-Point','GradObj','on','GradCon','off','DerivativeCheck','off', 'TolCon',1e-4);
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
%myDerivCheck(f,theta0,1e-5)

tic
% first step
[results1]=get_results(dtable,theta0);
print_results(results1);

% update weight matrix and produce second-step
[W,~]=getCovariance(results1.theta,dtable,draws);
[results2]=get_results(dtable,results1.theta);
print_results(results2);
toc

    %For a fixed weighting matrix (W) this minimizes the GMM objective
    %this works much better if knitromatlab is installed
    function [res]=get_results(tableA,x0)
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
        
        % Put the results into structure
        [~,SEest]=getCovariance(that,dtable,draws);
        res.fval = fval; res.beta=beta; res.theta = that; 
        res.delta=delta; res.resid = resid; res.SE=full(SEest);
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
    end

    function [W2,SE]=getCovariance(theta,dtable,draws)
        % this extracts the parameters for your specification
        params1=get_params(theta,draws);
        [dstar,Jac]=solveAllShares(dtable,draws,params1,method);
        [bhat,uhat]=ivregression(dstar,X,Z,W);

        g=bsxfun(@times,uhat,Z);
        gstar=bsxfun(@minus,g,mean(g))./n;
        S=gstar'*gstar;
        G=Z'*[Jac X]./n;
        
        W2=inv(n.^2*S);
        SE=sqrt(diag(inv(G'*W*G)*G'*W*S*W*G*inv(G'*W*G)));
    end

end
