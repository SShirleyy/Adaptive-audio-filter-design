function [thetahat,xhat]=rls(x,y,N,lambda)

% [thetahat,xhat]=rls(x,y,N,lambda)
%
%	x			- Data sequence
%	y			- Data sequence
%	N			- Dimension of the parameter vector
%	lambda			- Forgetting factor
%	thetahat		- Matrix with estimates of theta. 
%				  Row n corresponds to time n-1
%	xhat			- Estimate of x for n=1
%
%
%
%  rls: Recursive Least-Squares Estimation
%
% 	Estimator: xhat(n)=Y^{T}(n)thetahat(n-1)
%
%	thetahat is estimated using RLS. 
%
%	Initalization:	P(0)=10000*I, thetahat(0)=0
%
%     
%     Author: 
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
M = length(y);
% Initialize P, xhat and thetahat
p = 10000*eye(N+1);
thetahat = zeros(M+1,N+1);
xhat = zeros(1,M+1);

% Loop

for n=1:M,

	% Generate Y(n). Set elements of Y that does not exist to zero
    Y = zeros(N+1,1);
    for k = 1:N+1
        if (n-k+1>0)
            Y(k) = y(n-k+1);
        else
            Y(k) = 0;
        end
    end

	% Estimate of x
    xhat(n+1) = Y'*thetahat(n,:)';

	% Update K
    k = p*Y/(lambda+Y'*p*Y);

	% Update P
    p = (p-k*Y'*p)/lambda;


	% Update the n+1 row in the matrix thetahat which in the 
	% notation in the Lecture Notes corresponds to thetahat(n)

	thetahat(n+1,:)= thetahat(n,:) + k'*(x(n)-xhat(n+1));
end

% Shift thetahat one step so that row n corresponds to time n

thetahat=thetahat(2:M+1,:);
xhat = xhat(2:end);
