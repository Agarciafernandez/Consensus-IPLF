function w_fusion=Consensus_fusion_coefficient(Nodes_type,W_c,N_it_c)
%Consensus on the fusion coefficient, parameter b ine Eq. (3) and (4) in
%G. Battistelli, L. Chisci, G. Mugnai, A. Farina and A. Graziano, "Consensus-Based Linear and Nonlinear Filtering,"
% in IEEE Transactions on Automatic Control, vol. 60, no. 5, pp. 1410-1415, May 2015

b=double(Nodes_type>0)';
for l=1:N_it_c
    b_backup=b;
    b=W_c*b_backup;
end

w_fusion=ones(1,length(b));
w_fusion(b>1e-6)=1./b(b>1e-6);