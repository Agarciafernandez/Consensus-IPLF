%Implementation of the consensus iterated posterior linearisation filter
%(IPLF) described in 

% A. F. García-Fernández, G. Battistelli, "Consensus iterated posterior
% linearisation filter for distributed state estimation" in IEEE Signal
% Processing Letters, 2925.

%The consensus for linear models and graph description parts are based on a code provided by Prof. Giorgio
%Battistelli and written by Dr. Xiangxiang Dong

%Author: Ángel F. García-Fernández




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

%Consensus_type:
%1 -> Consensus on information (CI)
%2 -> Consensus on measurements (CM) with fusion coefficient given by Nnodes
%3 -> CM with fusion coefficient obtained by running consensus to
%approximate parameter b, see [2].
%4 -> Hybrid consensus on measurements and consensus on information (HCMCI)
%with fusion coefficient given by Nnodes.
%5 -> HCMCI with fusion coefficients obtained by running consensus to
%approximate parameter b, see [2].

%[2] G. Battistelli, L. Chisci, G. Mugnai, A. Farina and A. Graziano, "Consensus-Based Linear and Nonlinear Filtering," in IEEE Transactions on Automatic Control, vol. 60, no. 5, pp. 1410-1415, May 2015

Consensus_type=5;

if(Consensus_type==1)
    w_fusion=ones(1,Nnodes);
elseif(or(Consensus_type==2,Consensus_type==4))
    w_fusion=Nnodes*ones(1,Nnodes);
else
    w_fusion=Consensus_fusion_coefficient(Nodes_type,W_c,N_it_c);
end



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

        %We perform consensus on priors (CI and HCMCI methods) outside the loop over SLRs, as this
        %part does not change with the iterations.
        if(Consensus_type==1 || Consensus_type==4 || Consensus_type==5)
            [yk_pred,Yk_pred]=Consensus_mean_cov(yk_pred,Yk_pred,W_c,N_it_c,Nnodes,Nx);
        end


        %Iterated posterior linearisation iterations
        meank_j=meank_pred;
        Pk_j=Pk_pred;

        delta_q=zeros(Nx,Nnodes);
        delta_Omega=zeros(Nx,Nx,Nnodes);

        for p=1:N_it_iplf

            %Statistical linear regression done by each sensor
            for q=1:Nsensors
                x_s=Sensors_pos(:,q); %Position of the q-th sensor
                %Index of node corresponding to this sensor
                index_node=index_sensors(q);

                %We calculate the statistical linear regression (SLR). Note that Lambda is usually called Omega in IPLF literature (the resulting
                %MSE matrix), but here it is called Lambda, as in consensus papers, Omega
                %is the information matrix
                [A,b,Lambda]=SLR_measurement_range_bearings_VM(meank_j(:,index_node),Pk_j(:,:,index_node),weights,W0,Nx,Nz,kappa,Ap_kappa,x_s);

                Lambda=Lambda+R_kf;

                %We calculate delta_q and delta_Omega
                Inf_gain=A'/Lambda;
                delta_Omega(:,:,index_node)=Inf_gain*A;
                delta_q(:,index_node)=Inf_gain*(z_k(:,q)-b);

                %If no consensus is applied, we can use uncomment this line
                %and comment the one for consensus
                %                 [meank_j(:,index_node),Pk_j(:,:,index_node)]=linear_kf_update(meank_pred(:,index_node),Pk_pred(:,:,index_node),...
                %                     A,b,Lambda,0,z_k(:,q));

            end

            %Consensus for the linearised model
            [meank_j,Pk_j]=Consensus_linear(yk_pred,Yk_pred,delta_q, delta_Omega,W_c,w_fusion,N_it_c,Nnodes,Nx);



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

save(['Consensus_iplf_N_it_c',int2str(N_it_c),'_N_it_iplf',int2str(N_it_iplf),'_Consensus_type',...
    int2str(Consensus_type),'_kappa',int2str(kappa),'_R_range',int2str(R_range)])


figure(1)
plot(1:Nsteps,rms_position_error_t,'Linewidth',1.3)
grid on
xlabel('Time step')
ylabel('RMS position error')

