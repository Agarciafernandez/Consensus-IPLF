%Number of Monte Carlo runs
Nmc=100;
%Number of time steps
Nsteps=50;

%Dynamic model
T=2;
sigmaQ=4;

F=[1 T 0 0;0 1 0 0; 0 0 1 T;0 0 0 1];
Q=sigmaQ^2*kron(eye(2),[T^3/3 T^2/2; T^2/2 T]);
chol_Q=chol(Q)';


%Initial mean and covariance at time step 0
mean_ini=[2000,15,4000,5]';

%P_ini=diag([1500^2 16 1500^2 4]);

var_pos_x=1500^2;
var_pos_y=1500^2;
rho_pos=0.6;
cross_pos_xy=rho_pos*sqrt(var_pos_x)*sqrt(var_pos_y);

P_ini=[1500^2,0, cross_pos_xy, 0;...
    0, 4, 0,0;...
     cross_pos_xy,0, 1500^2,0;...
     0,0,0,4];

%P_ini=eye(4)

chol_P_ini=chol(P_ini)';

%We generate a ground truth trajectory
X_truth=zeros(4,Nsteps);
xk=mean_ini+chol_P_ini*randn(4,1);
X_truth(:,1)=xk;
for k=1:Nsteps-1
    xk_pred=F*xk+chol_Q*randn(4,1);
    X_truth(:,k+1)=xk_pred;
    xk=xk_pred;
end

%Measurement model parameters
kappa=600;
R_range=100;

%Network parameters. We load the same network as in 
%X. Dong, G. Battistelli, L. Chisci, and Y. Cai, "An event-triggered hybrid
%consensus filter for distributed sensor network," IEEE Signal Processing
%Letters, vol. 29, pp. 1472–1476, 2022.

load Network_1.mat
Nnodes= size(Network.Nodes.loc,1);  %Number of nodes
Nodes_pos=(Network.Nodes.loc)'; %Node positions

%Node type: Sensor node (type=1) and communication node(type=0)
Nodes_type=zeros(1,Nnodes); 
Nodes_type(2)=1;Nodes_type(38)=1;Nodes_type(16)=1;Nodes_type(27)=1;Nodes_type(11)=1; 
Nodes_type(20)=1;Nodes_type(36)=1;Nodes_type(9)=1;Nodes_type(42)=1;Nodes_type(28)=1; 

Nsensors= sum(Nodes_type>0);
Sensors_pos=Nodes_pos(:,Nodes_type>0);
index_sensors=find(Nodes_type>0);


%  Node_types(:)=1;  % All nodes are sensor nodes

%We calculate adjacency and degree matrix.
[A_matrix, D_matrix]=Calculate_adjacency(Network,Nnodes);

%Consensus weight matrix (using Metropolis weights)
W_c=zeros(Nnodes,Nnodes);
for i=1:Nnodes
    sum_w=0;
    for j=1:Nnodes
        if A_matrix(i,j)==1
            W_c(i,j)=(1+max(D_matrix(i,i),D_matrix(j,j)))^(-1);
            sum_w=sum_w+W_c(i,j);
        end
    end
    W_c(i,i)=1-sum_w;
end


% Consensus parameters
N_it_c=10; %Number of consensus iterations (L in the paper)

%IPLF parameters
N_it_iplf=3; %Number of IPLF iterations (J in the paper)



% Plot the scenario
% figure(2)
% clf
% hold on;
% plot(X_truth(1,:),X_truth(3,:),'LineWidth',0.8,'Color','k');% True trajectory
% plot(X_truth(1,1),X_truth(3,1),'g>','MarkerSize',6,'MarkerFaceColor','g'); %Start Position
% plot(X_truth(1,end),X_truth(3,end),'k>','MarkerSize',6,'MarkerFaceColor','k'); %End Position
% xlabel('x (m)');ylabel('y (m)');
% 
% % Plot Network Nodes (Sensor Node & Communication Node)
% 
% for k=1:Nnodes %plot Sensor Nodes & Communication Nodes 
%     if Nodes_type(k)==1
%         plot(Nodes_pos(1,k),Nodes_pos(2,k),'rs','MarkerSize',10,'MarkerFaceColor','r');
%     else
%         plot(Nodes_pos(1,k),Nodes_pos(2,k),'mo','MarkerSize',6,'MarkerFaceColor','m');
%     end
%    % text(posOfSensors(1,k),posOfSensors(2,k)-300,num2str(k)); %Node type
% end
% 
% % Plot Network Edges 
% for i=1:Nnodes 
%     for j=1:i
%         if (A_matrix(i,j)>0 && i~=j)
%             plot([Nodes_pos(1,i) Nodes_pos(1,j)],[Nodes_pos(2,i) Nodes_pos(2,j)],'LineWidth',0.8,'Color','b');
%         end
%     end
% end
% legend('Target','Initial position','End position','Com. node','Sensor node');
% grid on



function [A_matrix,D_matrix] = Calculate_adjacency(Network,Nnodes)
%This function calculates adjacency matrix and degree matrix
%The degree matrix is required to calculate the Metropolis weights

lap_matrix = zeros(Nnodes);  %Laplace Matrix: diag(d1,d2,...dn) with di the degree of Node i
D_matrix = zeros(Nnodes);% Degree matrix 
for i=1:Nnodes
    idx = Network.Nodes.neighbors{i};  % the ID of joint node 
    D_matrix(i,i) = size(idx,2); %Degree matrix, the diagonal elements 
    lap_matrix(i,idx) = -1;   
    lap_matrix(i,i) = length(idx);
end
A_matrix = D_matrix-lap_matrix; %Adjacency matrix
end
