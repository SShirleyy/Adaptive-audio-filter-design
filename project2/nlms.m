function [thetahat,xhat]=nlms(x,y,N,muu)

% [thetahat,xhat]=nlms(x,y,N,muu)
%
%	x			- Data sequence
%	y			- Data sequence
%	N			- Model order 
%	muu			- Step size
%	thetahat		- Matrix with estimates of theta. 
%				  Row n corresponds to the estimate thetahat(n)'
%	xhat			- Estimate of x
%
%
%
%  nlms: The Normalized Least-Mean Square Algorithm
%
% 	Estimator: xhat(n)=Y^{T}(n)thetahat(n-1)
%
%	thetahat is estimated using NLMS. 
%
%     
%     Author: 
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Initialize xhat and thetahat
M = length(y);
thetahat = zeros(M+1,N+1);
xhat = zeros(1,M+1);
c = 1;
% Loop

for n=1:M,

	% Generate Y. Set elements of Y that does not exist to zero
    Y = zeros(1,N+1);
    for k = 1:N+1
        if (n-k+1>0)
            Y(k) = y(n-k+1);
        else
            Y(k) = 0;
        end
    end
    
    xhat(n+1) = Y*thetahat(n,:)';

	% Estimate of x


	% Update the n+1 row in the matrix thetahat which in the notation in the Lecture Notes
	% corresponds to thetahat(n)

	thetahat(n+1,:)=thetahat(n,:) + muu/(c+norm(Y)^2)*Y*(x(n)-xhat(n+1));
end

% Shift thetahat one step so that row n corresponds to time n
xhat = xhat(2:end);
thetahat=thetahat(2:M+1,:);
