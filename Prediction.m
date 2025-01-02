function [meank_pred,Pk_pred]=Prediction(meank,Pk,F,Q,Nnodes)

%Kalman filter prediction step for the filter of each node
meank_pred=zeros(size(meank));
Pk_pred=zeros(size(Pk));

for n=1:Nnodes
    meank_n=meank(:,n);
    Pk_n=Pk(:,:,n);

    meank_pred(:,n)=F*meank_n;
    Pk_pred(:,:,n)=F*Pk_n*F'+Q;    

end
