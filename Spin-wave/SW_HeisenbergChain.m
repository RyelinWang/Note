%%%%% Two-sublattice spin-wave calculation for spin-1/2 Heisenberg chain %%%%%
%%% Output: A figure, showing the E(k) in first magnetic BZ.

%% Model parameter setup
J=1.0; %Strength of exchange interaction, and positive characterize AFM type.

%% Gell-mann matrix setup
S=zeros(2,2,3);
S_A=zeros(2,2,3);
S_B=zeros(2,2,3);
S(:,:,1)=[0,1;1,0]/2;
S(:,:,2)=[0,-1i;1i,0]/2;
S(:,:,3)=[1,0;0,-1]/2;
R_A=eye(2,2);
R_B=[0,1;1,0];
for m=1:3
    S_A(:,:,m)=R_A'*S(:,:,m)*R_A;
    S_B(:,:,m)=R_B'*S(:,:,m)*R_B;
end
%% Hamiltonian cell setup
V_AB=zeros(1,1);
V_BA=zeros(1,1);
T_AB=zeros(1,1);
Delta_AB=zeros(1,1);
Z=2;
for m=1:3
    V_AB=V_AB+J*(S_A(2,2,m)-S_A(1,1,m)*eye(1,1))*S_B(1,1,m);
    V_BA=V_BA+J*(S_B(2,2,m)-S_B(1,1,m)*eye(1,1))*S_A(1,1,m);
    for alpha=1:1
        for beta=1:1
            T_AB(alpha,beta)=T_AB(alpha,beta)+J*S_A(alpha+1,1,m)*S_B(1,beta+1,m);
            Delta_AB(alpha,beta)=Delta_AB(alpha,beta)+J*S_A(alpha+1,1,m)*S_B(beta+1,1,m);
        end
    end
end
delta=[1;-1]'/2;
%delta=[1,0;-1,1;0,-1]';
%% Compute
Size=500;
K_grid_x=zeros(1,Size+1); %Storage k_x for plot
E=zeros(Size+1,2); %Storage energy for magnons
for i=1:Size+1
    %Compute location for k
    k=1/Size*(i-1)-1/2;
    K_grid_x(i)=k;
    %Structure factor
    gamma=0.0; 
    for nn=1:2
        gamma=gamma+exp(-1i*2*pi*k*delta(:,nn));
    end
    %Hamiltonian construct
    H_1=[Z*V_AB,T_AB*gamma;T_AB'*conj(gamma),Z*V_BA];
    H_2=[Z*V_AB.',T_AB'*gamma;T_AB*conj(gamma),Z*V_BA.'];
    Delta=[zeros(1,1),Delta_AB*gamma;Delta_AB.'*conj(gamma),zeros(1,1)];
    H=[H_1,Delta;Delta',H_2];
    Sigma=[eye(2,2),zeros(2,2);zeros(2,2),-eye(2,2)];
    [Vec,EigenAll]=eig(Sigma*H);
    E_eigen=[];
    for m=1:4
        if Vec(:,m)'*Sigma*Vec(:,m)>0
            E_eigen=[E_eigen,EigenAll(m,m)];
        end
    end
    if(numel(E_eigen)~=2)
        fprintf("k=(%d),no solution\n",i);
        E_eigen=nan*ones(2,1);
    end
    if(sum(abs(imag(E_eigen)))>sum(abs(real(E_eigen)))*1e-3&&sum(abs(E_eigen))>1e-3)
        fprintf("k=(%d), energy instable\n",i);
        E(i,:)=nan*ones(2,1);
    else
        E_eigen=sort(real(E_eigen));
        %E(i,j,:)=E_eigen(numel(E_eigen)/2+1:numel(E_eigen));
        %E(i,j,:)=E_eigen(1:numel(E_eigen)/2);
        E(i,:)=E_eigen;
    end
end
%% Plot
figure
hold on
plot(K_grid_x,E(:,1),"LineWidth",2)
plot(K_grid_x,E(:,2),"LineWidth",2);
title("SU(2) spin wave","FontSize",20)
