function z_t=Measurement_sampling_all_times(X_truth,Nsteps,Nsensors,Sensors_pos,kappa,R_range)
Nz=3;

z_t=zeros(Nz,Nsteps,Nsensors);
for q=1:Nsensors
    x_s=Sensors_pos(:,q); %Position of the q-th sensor
    for k=1:Nsteps
        state_k=X_truth(:,k);

        %Angular measurement
        mean_z=[state_k(1)-x_s(1);state_k(3)-x_s(2)];
        %mean_z_aux=mean_z/sqrt(sum(mean_z.^2));
        mean_theta=atan2(mean_z(2),mean_z(1));

        %We generate von Mises meausurement
        z_real_angle = circ_vmrnd(mean_theta, kappa, 1);
        z_real(1:2)=[cos(z_real_angle);sin(z_real_angle)];


        mean_z_range=sqrt((state_k(1)-x_s(1))^2+(state_k(3)-x_s(2))^2);
        z_real(3)=mean_z_range+sqrt(R_range)*randn(1);

        z_t(:,k,q)=z_real;

    end
end
end