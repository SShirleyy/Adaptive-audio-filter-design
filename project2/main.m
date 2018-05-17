close all;
clear all;


play_sound = 0;

% Read audio file
[z,fs] = audioread('EQ2400project2data2017.wav');
if (play_sound)
    soundsc(z,fs);
end

d = 200;
z1 = z(d:end);

%% LMS
% optimal order 4 & optimal stepsize 0.37
[theta_lms,x_lms] = lms(z,z1,10,0.37);
z_lms = z(1:end-d+1) - x_lms;

%% NLMS
% optimal order 11 & optimal stepsize 2.2
[theta_nlms,x_nlms] = nlms(z,z1,11,2.2);
z_nlms = z(1:end-d+1) - x_nlms';

%% RLS
% optimal delay 200 & optimal order 3 & optimal lambda 0.82
[theta_rls,x_rls] = rls(z,z1,3,0.82);
z_rls = z(1:end-d + 1) - x_rls';



%% Kalman
R2 = 0.001;
N = 30;
M = length(z1);
P = eye(N+1);
theta_kal = zeros(M+1,N+1);
R1 = eye(N+1); 


for n=1:M,
    % Generate Y. Set elements of Y that does not exist to zero
    Y = zeros(N+1,1);
    for k = 1:N+1
        if ((n-k+1)>0)
            Y(k)=z1(n-k+1);
        else
            Y(k)=0;
        end
    end

    x_kal(n+1) = Y'*theta_kal(n,:)';
    K=P*Y/(Y'*P*Y+R2);
    P=P-K*Y'*P+R1;
    theta_kal(n+1,:)= theta_kal(n,:)+K'*(z(n)-x_kal(n+1));
end
x_kal = x_kal(2:M+1);
z_kal = z(1:end-d+1)-x_kal';

pause;
soundsc(z,fs);
pause;
soundsc(z_lms,fs);

pause;
soundsc(z,fs);
pause;
soundsc(z_nlms,fs);

pause;
soundsc(z,fs);
pause;
soundsc(z_rls,fs);

pause;
soundsc(z,fs);
pause;
soundsc(z_kal,fs);

% Below are codes we used to generate plot and determine parameters, which
% are not needed in presentation



% % determine delay
% for d = 1:300
%     z1 = z(d:end);
%     [~,x_lms] = lms(z,z1,4,0.33);
%     mse_lms_delay(d) = mean((z(1:end-d+1)-x_lms).^2);
% end
% figure(9);
% plot(1:300,mse_lms_delay);
% title('Mean Square Error of LMS filter with delay');
% xlabel('delay')
% ylabel('MSE')
% d = 200;
% z1 = z(d:end);
% 
% % detemine step size of lms
% m = linspace(0,0.6,61);
% for in = 1:length(m)
% [~,x_lms] = lms(z,z1,10,m(in));
% mse_lms_mu(in) = mean((z(1:end-d+1)-x_lms).^2);
% end
% figure(1);
% plot(m,mse_lms_mu);
% title('Mean Square Error of LMS filter with stepsize');
% xlabel('stepsize')
% ylabel('MSE')
% % 
% % detemine order of lms
% m = linspace(0,15,16);
% for in = 1:length(m)
% [~,x_lms] = lms(z,z1,m(in),0.37);
% mse_lms_or(in) = mean((z(1:end-d+1)-x_lms).^2);
% end
% figure(2);
% plot(m,mse_lms_or);
% title('Mean Square Error of LMS filter with order');
% xlabel('order')
% ylabel('MSE')
% 
% 





% % determine delay
% for d = 1:300
%     z1 = z(d:end);
%     [~,x_nlms] = nlms(z,z1,4,1.6);
%     mse_nlms_delay(d) = mean((z(1:end-d+1)-x_nlms').^2);
% end
% figure(8);
% plot(1:300,mse_nlms_delay);
% title('Mean Square Error of NLMS filter with delay');
% xlabel('delay')
% ylabel('MSE')
% d = 200;
% z1 = z(d:end);
% % detemine step size of nlms
% m = linspace(0,10,50);
% for in = 1:length(m)
% [~,x_nlms] = nlms(z,z1,4,m(in));
% mse_nlms_mu(in) = mean((z(1:end-d+1)-x_nlms').^2);
% end
% figure(3);
% plot(m,mse_nlms_mu);
% title('Mean Square Error of NLMS filter with stepsize');
% xlabel('stepsize')
% ylabel('MSE')
% 
% % detemine order of nlms
% m = linspace(0,100,101);
% for in = 1:length(m)
% [~,x_nlms] = nlms(z,z1,m(in),2.2);
% mse_nlms_or(in) = mean((z(1:end-d+1)-x_nlms').^2);
% end
% figure(4);
% plot(m,mse_nlms_or);
% title('Mean Square Error of NLMS filter with order');
% xlabel('order')
% ylabel('MSE')
% 
% 



% lambda = linspace(0,1,101);
% for in = 1:length(lambda)
% [~,x_rls] = rls(z,z1,4,lambda(in));
% mse_rls_lambda(in) = mean((z(1:end-d+1)-x_rls').^2);
% end
% figure(5);
% plot(lambda,mse_rls_lambda);
% title('Mean Square Error of RLS filter with lambda');
% xlabel('lambda')
% ylabel('MSE')
% 
% %detemine order of rls
% m = linspace(0,100,101);
% for in = 1:length(m)
% [~,x_rls] = rls(z,z1,m(in),0.82);
% mse_rls_or(in) = mean((z(1:end-d+1)-x_rls').^2);
% end
% figure(6);
% plot(m,mse_rls_or);
% title('Mean Square Error of RLS filter with order');
% xlabel('order')
% ylabel('MSE')
% 
% 
% % detemine step
% for d = 1:300
%     z1 = z(d:end);
%     [~,x_rls] = rls(z,z1,4,0.99);
%     mse_rls_delay(d) = mean((z(1:end-d+1)-x_rls').^2);
% end
% figure(7);
% plot(1:300,mse_rls_delay);
% title('Mean Square Error of RLS filter with delay');
% xlabel('order')
% ylabel('MSE')
% d = 200;
% z1 = z(d:end);




% N = 5;
% M = length(z1);
% theta_kal = zeros(M+1,N+1);
% 
% noise_ratio = linspace(0,0.1,50); 
% mse_kal_ratio = zeros(length(noise_ratio),1);
% R1 = eye(N+1); 
% 
% for in = 1:length(noise_ratio)
%     P = eye(N+1);
%     R2 = noise_ratio(in);
% 
%     for n=1:M,
%         % Generate Y. Set elements of Y that does not exist to zero
%         Y = zeros(N+1,1);
%         for k = 1:N+1
%             if ((n-k+1)>0)
%                 Y(k)=z1(n-k+1);
%             else
%                 Y(k)=0;
%             end
%         end
% 
%         x_kal(n+1) = Y'*theta_kal(n,:)';
%         K=P*Y/(Y'*P*Y+R2);
%         P=P-K*Y'*P+R1;
%         theta_kal(n+1,:)= theta_kal(n,:)+K'*(z(n)-x_kal(n+1));
%     end
%     x_kal = x_kal(2:M+1);
%     mse_kal_ratio(in) = mean((z(1:end-d+1)-x_kal').^2);
% end
% 
% figure(10)
% plot(noise_ratio,mse_kal_ratio);
% title('Mean Square Error of Kalman filter with ratio');
% xlabel('noise ratio')
% ylabel('MSE')
% 
% 
% 
% 
% R2 = 0.01;
% M = length(z1);
% 
% order = linspace(0,50,51); 
% mse_kal_order = zeros(length(order),1);
% 
% for in = 1:length(order)
%     N = order(in);
%     P = eye(N+1);
%     theta_kal = zeros(M+1,N+1);
%     R1 = eye(N+1); 
% 
% 
%     for n=1:M,
%         % Generate Y. Set elements of Y that does not exist to zero
%         Y = zeros(N+1,1);
%         for k = 1:N+1
%             if ((n-k+1)>0)
%                 Y(k)=z1(n-k+1);
%             else
%                 Y(k)=0;
%             end
%         end
% 
%         x_kal(n+1) = Y'*theta_kal(n,:)';
%         K=P*Y/(Y'*P*Y+R2);
%         P=P-K*Y'*P+R1;
%         theta_kal(n+1,:)= theta_kal(n,:)+K'*(z(n)-x_kal(n+1));
%     end
%     x_kal = x_kal(2:M+1);
%     mse_kal_order(in) = mean((z(1:end-d+1)-x_kal').^2);
% end
% 
% figure(11)
% plot(order,mse_kal_order);
% title('Mean Square Error of Kalman filter with order');
% xlabel('order')
% ylabel('MSE')