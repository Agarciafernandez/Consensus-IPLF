function [yk_pred,Yk_pred]=Consensus_mean_cov(yk_pred,Yk_pred,W_c,N_it_c,Nnodes,Nx)

%We run consensus for a mean and a covariance/information matrix


for l=1:N_it_c
    y_backup=yk_pred;
    Y_backup=Yk_pred;
    for i=1:Nnodes
        yk_pred(:,i)=zeros(Nx,1);
        Yk_pred(:,:,i)=zeros(Nx,Nx);
        for j=1:Nnodes
            if(W_c(i,j)>0)
                yk_pred(:,i)=yk_pred(:,i)+W_c(i,j)*y_backup(:,j);
                Yk_pred(:,:,i)=Yk_pred(:,:,i)+W_c(i,j)*Y_backup(:,:,j);
            end
        end
    end
end
end