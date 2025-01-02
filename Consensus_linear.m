function [meank_upd,Pk_upd]=Consensus_linear(yk_pred,Yk_pred,delta_q, delta_Omega,W_c,w_fusion,N_it_c,Nnodes,Nx)

%Consensus algorithms when the measurement model is linearised

%yk_pred and Yk_pred are the output of the consensus on priors if consensur on priors is
%used. Otherwise, they correspond to the predicted information matrix and
%vector.

%Consensus on measurements (consensus on delta_q and delta_Omega)
[delta_q,delta_Omega]=Consensus_mean_cov(delta_q,delta_Omega,W_c,N_it_c,Nnodes,Nx);

%We calculate the updated mean and covariance matrix for each node
meank_upd=zeros(Nx,Nnodes);
Pk_upd=zeros(Nx,Nx,Nnodes);


%Fusion HCMI
for i=1:Nnodes
    y=yk_pred(:,i)+w_fusion(i)*delta_q(:,i);
    Y=Yk_pred(:,:,i)+w_fusion(i)*delta_Omega(:,:,i);
    Pk_upd(:,:,i)=inv(Y);
    meank_upd(:,i)=Pk_upd(:,:,i)*y;
end



end


