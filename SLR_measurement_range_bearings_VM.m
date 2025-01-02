function [A_l,b_l,Omega_l]=SLR_measurement_range_bearings_VM(meank,Pk,weights,W0,Nx,Nz,kappa,Ap_kappa,x_s)
%We make SLR of conditional moments for Von Mises distribution
% Á. F. García-Fernández, F. Tronarp and S. Särkkä, "Gaussian Target Tracking With Direction-of-Arrival von Mises–Fisher Measurements," 
% in IEEE Transactions on Signal Processing, vol. 67, no. 11, pp. 2960-2972, 1 June1, 2019, doi: 10.1109/TSP.2019.2911258.


%First we compute the moments of h(x), which projects to the unit circle
%using sigma-points

Nz_angle=Nz-1;

%Sigma-point generation
chol_var_mult=chol((Nx/(1-W0)*Pk));
sigma_points=[zeros(Nx,1),chol_var_mult',-chol_var_mult'];
sigma_points=repmat(meank,1,length(weights))+sigma_points;


%First the angle
%Transformation through h
dist_sigma_points=sqrt((sigma_points(1,:)-x_s(1)).^2+(sigma_points(3,:)-x_s(2)).^2);
sigma_points_h=[sigma_points(1,:)-x_s(1);sigma_points(3,:)-x_s(2)]./repmat(dist_sigma_points,2,1);
mean_h=sigma_points_h*weights';
mean_hh=zeros(Nz_angle);%E[h*h']
cov_xh=zeros(Nx,Nz_angle);

for j=1:length(weights)
    mean_hh=mean_hh+weights(j)*(sigma_points_h(:,j)*sigma_points_h(:,j)');
    cov_xh=cov_xh+weights(j)*((sigma_points(:,j)-meank)*(sigma_points_h(:,j)-mean_h)');
end

%Now we can compute the moments of the transformed variable (using the VM
%adjustments)
p=2; %p=2 as we are on 2D
z_pred_ukf=Ap_kappa*mean_h;
% var_pred_ukf=Ap_kappa/kappa*(eye(Nz)-(p+kappa*Ap_kappa-kappa/Ap_kappa)*mean_hh)...
%     +(Ap_kappa)^2*(mean_hh-mean_h*mean_h');
var_pred_ukf=(Ap_kappa)^2*(mean_hh-mean_h*mean_h')+...
    +Ap_kappa/kappa*eye(Nz_angle)+(1-(Ap_kappa)^2-p*Ap_kappa/kappa)*mean_hh';
var_xz_ukf=Ap_kappa*cov_xh;



%Statistical linearisaltion
A_l_angle=var_xz_ukf'/Pk;
b_l_angle=z_pred_ukf-A_l_angle*meank;
Omega_l_angle=var_pred_ukf-A_l_angle*Pk*A_l_angle';


%Now the range
z_pred_ukf=dist_sigma_points*weights';

var_pred_ukf=zeros(1);
var_xz_ukf=zeros(Nx,1);

for j=1:length(weights)
    sub_z_j=dist_sigma_points(:,j)-z_pred_ukf;
    var_pred_ukf=var_pred_ukf+weights(j)*(sub_z_j*sub_z_j');
    var_xz_ukf=var_xz_ukf+weights(j)*(sigma_points(:,j)-meank)*sub_z_j';
end

%Statistical linearisaltion

A_l_r=var_xz_ukf'/Pk;
b_l_r=z_pred_ukf-A_l_r*meank;
Omega_l_r=var_pred_ukf-A_l_r*Pk*A_l_r';


A_l=[A_l_angle;A_l_r];
b_l=[b_l_angle;b_l_r];
Omega_l=zeros(3:3);
Omega_l(1:2,1:2)=Omega_l_angle;
Omega_l(3,3)=Omega_l_r;




