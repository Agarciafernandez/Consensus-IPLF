%Implementation of the centralised iterated posterior linearisation filter
%(IPLF) used for the paper 

% A. F. García-Fernández, G. Battistelli, "Consensus iterated posterior
% linearisation filter for distributed state estimation" in IEEE Signal
% Processing Letters, 2025.


clear
rand('seed',8)
randn('seed',8)

Consensus_scenario;


%Sigma-point parameters (according to the unscented transform)
Nx=4;
W0=1/3; %Weight at the mean
Wn=(1-W0)/(2*Nx);
weights=[W0,Wn*ones(1,2*Nx)];

%Ap_kappa (correction factor required for the von Mises distribution)
p=2; %p=2 as we are on the circle
Ap_kappa=besseli(p/2,kappa)./besseli(p/2-1,kappa);
index=isnan(Ap_kappa);
Ap_kappa(index)=1-(p-1)./(2*kappa(index)); %For large kappa we apply Filip´s expansion

%Additive noise of the KF in the range measurement
R_kf=zeros(3);
R_kf(3,3)=R_range;

%Dimension of the measurement

Nz=3;

%Note that Lambda_m is usually called Omega in IPLF literature (Resulting
%MSE matrix), but here it is called Lambda, as in consensus papers, Omega
%is the information matrix

rand('seed',8)
randn('seed',8)

%Square error of position elements considering all sensors at each time
%step
square_pos_error_t=zeros(1,Nsteps);

%Monte Carlo runs

for i=1:Nmc

    %We generate the measurements at all time steps for this Monte Carlo
    %run
    z_t=Measurement_sampling_all_times(X_truth,Nsteps,Nsensors,Sensors_pos,kappa,R_range);

    tic

    %Propagated mean at each time step for each node filter (all of them start
    %with the same prior at time step 1)
    meank_pred=repmat(mean_ini,1,Nnodes);
    Pk_pred=repmat(P_ini,1,1,Nnodes);


    for k=1:Nsteps
        %Update step

        %Measurement at time step k
        z_k=squeeze(z_t(:,k,:));

        %We obtain the predicted information matrix and vector
        Yk_pred=zeros(Nx,Nx,Nnodes);
        yk_pred=zeros(Nx,Nnodes);
        for q=1:Nnodes
            Yk_pred(:,:,q)=inv(Pk_pred(:,:,q));
            yk_pred(:,q)=Yk_pred(:,:,q)*meank_pred(:,q);
        end


        %Iterated posterior linearisation iterations
        meank_j=meank_pred;
        Pk_j=Pk_pred;

        delta_q=zeros(Nx,Nnodes);
        delta_Omega=zeros(Nx,Nx,Nnodes);

        for p=1:N_it_iplf

            %These variables will contain the sum of the information
            %vectors and matrices, which are required for the centralised
            %solution
            delta_q_tot=zeros(Nx,1);
            delta_Omega_tot=zeros(Nx,Nx);


            %Statistical linear regression done by each sensor
            for q=1:Nsensors
                x_s=Sensors_pos(:,q); %Position of the q-th sensor
                %Index of node corresponding to this sensor
                index_node=index_sensors(q);

                %We calculate the SLR (Note that Lambda is usually called Omega in IPLF literature (Resulting
                %MSE matrix), but here it is called Lambda, as in consensus papers, Omega
                %is the information matrix
                
                [A,b,Lambda]=SLR_measurement_range_bearings_VM(meank_j(:,index_node),Pk_j(:,:,index_node),weights,W0,Nx,Nz,kappa,Ap_kappa,x_s);

                Lambda=Lambda+R_kf;

                %We calculate delta_q and delta_Omega
                Inf_gain=A'/Lambda;
                delta_Omega(:,:,index_node)=Inf_gain*A;
                delta_q(:,index_node)=Inf_gain*(z_k(:,q)-b);


                delta_q_tot=delta_q_tot+delta_q(:,index_node);
                delta_Omega_tot=delta_Omega_tot+delta_Omega(:,:,index_node);
            end

            %Update for each node using the centralised solution
            meank_j=zeros(Nx,Nnodes);
            Pk_j=zeros(Nx,Nx,Nnodes);
            for q=1:Nnodes
                Yk_upd_q=Yk_pred(:,:,q)+delta_Omega_tot;
                yk_upd_q=yk_pred(:,q)+delta_q_tot;

                Pk_j(:,:,q)=inv(Yk_upd_q);
                meank_j(:,q)=Pk_j(:,:,q)*yk_upd_q;


            end

        end


        meank=meank_j;
        Pk=Pk_j;

        %Square error calculation
        for q=1:Nsensors
            square_pos_error_t(k)=square_pos_error_t(k)+(meank(1,q)-X_truth(1,k))^2+(meank(3,q)-X_truth(3,k))^2;
        end


        %Prediction step
        [meank_pred,Pk_pred]=Prediction(meank,Pk,F,Q,Nnodes);

    end
    t=toc;
    display(['Completed iteration number ', num2str(i),' time ', num2str(t), ' sec'])


end

%RMS position error at each time step
rms_position_error_t=sqrt(square_pos_error_t/(Nmc*Nnodes));

%RMS position error across all time steps
rms_position_error=sqrt(sum(square_pos_error_t)/(Nsteps*Nmc*Nnodes))

save(['Centralised_iplf_N_it_iplf',int2str(N_it_iplf),'_kappa',int2str(kappa),'_R_range',int2str(R_range)])


figure(1)
plot(1:Nsteps,rms_position_error_t,'Linewidth',1.3)
grid on
xlabel('Time step')
ylabel('RMS position error')

