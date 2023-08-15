clc
clear all %#ok
close all
load("Data2_for_calibration.mat");

% Nominal Value
l__1 = 0.006;
l__2 = 0.176;
l__3 = 0.050;
l__4 = 0.115;
%% nominal 
L = l__1 + l__2 + l__3 + l__4 ;  %% 0.347
theta=42*(pi/180); %% 0.7333 Rad
alpha=50*(pi/180); %% 0.8727 Rad   
%%% rms = 8.6434
%% factorGrapgh 
L = 0.347019 ;  
theta = 42.0017*(pi/180); 
alpha = (180-130.026)*(pi/180); %% 49.9740 deg   

%% factorGrapgh Fk fix 
% L = 0.346967 ;  %% 0.347
% theta = 40.002*(pi/180); %% 0.7333 Rad-
% alpha = (49.9908)*(pi/180); %% 0.8727 Rad   

% %factorGrapgh RCM fix 
L = 0.347027 ;  %% 0.347
theta = 42.0003*(pi/180); %% 0.7333 Rad-
alpha = (48.7512)*(pi/180); %% 0.8727 Rad   

%% fmincon 
L = 0.347 ;
theta = 0.733;
alpha =0.8722;
%%% rms = 1.5002
%% GA
L = 0.347 ;
theta = 0.7331; % 42.0035 deg
% alpha =0.8719;
alpha =0.8722; % 49.9734 deg
%%%% rms = 1.4951

%% PSO
L = 0.3468 ;
theta = 0.7325;
alpha = 0.8802;
%%%% rms = 1.5002
%% Pattern search
L = 0.3468 ;
theta = 0.7325;
alpha = 0.8802;
%%%% rms = 1.5002

%%
phi = out.q(:,1).*(pi/180); %% Rad
psi = out.q(:,2).*(pi/180); %% Rad
d = (-out.q(:,3))/1000;     %% m

x_fk = zeros(size(d,1),1);  
y_fk = zeros(size(d,1),1);
z_fk = zeros(size(d,1),1);

for i = 1 : size(d,1)
    x_fk(i,1) = (sin(theta) * cos(psi(i) + alpha) * d(i) + cos(theta) * cos(phi(i)) * sin(psi(i) + alpha) * d(i) + sin(theta) * L)*1000;
    y_fk(i,1) = (sin(phi(i))* sin(psi(i) + alpha) * d(i))*1000;
    z_fk(i,1) = (cos(theta) * cos(psi(i) + alpha) * d(i) - cos(phi(i)) * sin(theta) * sin(psi(i) + alpha) * d(i) + cos(theta) * L)*1000;
end

x_rfk = squeeze(out.Position(1,:,:)).*1000;
y_rfk = squeeze(out.Position(2,:,:)).*1000;
z_rfk = squeeze(out.Position(3,:,:)).*1000;

x_GT = out.X ; %%mm
y_GT = out.Y ; %%mm
z_GT = out.Z ; %%mm

error_x = x_GT - x_fk ; %%mm
error_y = y_GT - y_fk ; %%mm
error_z = z_GT - z_fk ; %%mm

%error = 100*rms(error_x.^2 + error_y.^2 + error_z.^2)

[rms(error_x);rms(error_y);rms(error_z);]
subplot(3,1,1)
plot(error_x)
ylabel('mm');
legend('Error x')
grid minor
hold on

subplot(3,1,2)
plot(error_y)
ylabel('mm');
legend('Error Y')
grid minor 
hold on

subplot(3,1,3)
plot(error_z)
ylabel('mm');
legend('Error Z')
grid minor
hold on

% figure
% subplot(3,1,1)
% plot(x_fk)
% hold on 
% plot(x_rfk)
% ylabel('mm');
% grid minor
% 
% subplot(3,1,2)
% plot(y_fk)
% hold on 
% plot(y_rfk)
% ylabel('mm');
% grid minor
% 
% subplot(3,1,3)
% plot(z_fk)
% hold on 
% plot(z_rfk)
% ylabel('mm');
% grid minor

%%
for e = 1:1
    phi_Matrix = zeros(3*size(d,1),3);
    delta_r = zeros(3*size(d,1),1);
    for j=1:size(d,1)
        x_fk(j,1) = (sin(theta) * cos(psi(j) + alpha) * d(j) + cos(theta) * cos(phi(j)) * sin(psi(j) + alpha) * d(j) + sin(theta) * L)*1000;
        y_fk(j,1) = (sin(phi(j))* sin(psi(j) + alpha) * d(j))*1000;
        z_fk(j,1) = (cos(theta) * cos(psi(j) + alpha) * d(j) - cos(phi(j)) * sin(theta) * sin(psi(j) + alpha) * d(j) + cos(theta) * L)*1000;
       
        phi_Matrix ((j-1)*3+1:j*3,:) = [ -sin(theta)*sin(psi(j)+alpha)*d(j) + cos(theta)*cos(phi(j))*cos(psi(j) + alpha)*d(j) , -sin(theta)*cos(phi(j))*sin(psi(j) + alpha)*d(j) + cos(theta)*cos(psi(j) + alpha)*d(j) + L*cos(theta) ,   sin(theta)   ;...
            d(j)*sin(phi(j))*cos(psi(j) + alpha),     0    , 0   ;...
            -sin(theta)*cos(phi(j))*cos(psi(j) + alpha)*d(j) - cos(theta)*sin(psi(j) + alpha)*d(j)  ,   -cos(theta)*cos(phi(j))*sin(psi(j) + alpha)*d(j) - sin(theta)*cos(psi(j) + alpha)*d(j) - L*sin(theta)   ,  -cos(theta)*cos(phi(j))*sin(psi(j) + alpha)*d(j) - sin(theta)*cos(psi(j) + alpha)*d(j) - L*sin(theta) ];

        delta_r((j-1)*3+1:j*3,:) = [(x_GT(j)-x_fk(j));
                                    (y_GT(j)-y_fk(j));
                                    (z_GT(j)-z_fk(j))];
    end

    delat_phi = phi_Matrix'*pinv(phi_Matrix*phi_Matrix')*delta_r;

    alpha = alpha + delat_phi(1);
    theta = theta + delat_phi(2);
    L = L + delat_phi(3);
end