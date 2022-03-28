function [T,S,E,D,nD,U,Ratios,Age_structure, FinalState] = VoC_Vacc_ODEs2(M_from_to_H, M_from_to_O, alpha, gamma, sigma, d, Od, tau, nV_Beta, nV_Speed, HHQ, N0, V1, V2, VET1, VET2, VoC_Sus, Transmission_Reduction, Symptom_Reduction, RatioPf, MaxTime, InitialState) %#codegen
%%
%
%  Extended_ODEs( Mixing_Matrix, 1/Latent_Per, 1/Inf_Per, Susceptibility, Detectibility, Transmission,  Pop_Size, MaxTime, InitialState)   
%
%%

%
%   nV_Beta is the multiplicative on R (so multi on beta is nV_Speed*nV_Beta)
%   nV_Speed is the multiplicative on alpha & gamma.
%

%fprintf(1,'ODEs gamma=%g\n',gamma);
% Sets up parameters if necessary.

number_E_states=3;

L=length(N0);

if length(sigma)==1
    sigma=sigma+0*N0;
else
    sigma=reshape(sigma,L,1);
end
if length(d)==1
    d=d+0*N0;
else
    d=reshape(d,L,1);
end
if length(Od)==1
    Od=Od+0*N0;
else
    Od=reshape(Od,L,1);
end
if length(tau)==1
    tau=tau+0*N0;
else
    tau=reshape(tau,L,1);
end

if length(HHQ)==1
    hhq=HHQ*ones(size(N0,1),size(N0,2));
end


M_from_to=M_from_to_H+M_from_to_O;

M_from_to_HAT=zeros(L,L);
for f=1:L
    for t=1:L
        M_from_to_HAT(f,t)=M_from_to(f,t)*d(t)*sigma(t)*(1+tau(f)*(1-d(f))/d(f));
    end
end

[V,D]=eig(M_from_to_HAT);
[R0 i]=max(abs(diag(D))); R0=R0/gamma;
Age_structure=abs(V(:,i));

if InitialState(1)<0
   E0=Age_structure; D0=Age_structure; U0=Age_structure; S0=N0;
   InitialState=zeros(1,441); L=21;
   InitialState(1:L)=S0'; m=L;
   for i=1:number_E_states
       InitialState(m+[1:3*L])=[E0'/number_E_states 0*E0' 0*E0']; m=m+3*L;
   end
   InitialState(m+[1:4*L])=[D0' zeros(1,length(E0)*2) U0'];
   InitialState=[InitialState 0*InitialState((L+1):end) 0*InitialState((L+1):end) 0*InitialState(1:(7*L))]; % extra for newVariant, extra for VoC and 5 more for vaccination status, extra S classes and 3 immune class (V1,V2 & R).
end

options = odeset('RelTol', 1e-8);

m=number_E_states;

SP=3*(8*L+4*m*L) + L;
NVacc1=InitialState(SP+[1:L]);
NVacc2=InitialState(SP+L+[1:L]);
NVacc0=N0'-NVacc1-NVacc2;
V1(V1>NVacc0)=NVacc0(V1>NVacc0);
V2(V2>NVacc1)=NVacc1(V2>NVacc1);


[t, pop]=ode45(@Diff_3_4,[0:1:MaxTime],[InitialState],options,[L number_E_states reshape(M_from_to_H,1,[]) reshape(M_from_to_O,1,[]) sigma' d' Od' tau' alpha gamma N0' hhq' nV_Beta nV_Speed V1 V2 VET1 VET2 VoC_Sus ...
    Transmission_Reduction, Symptom_Reduction, RatioPf]);

T=t; S=pop(:,1:L); 

S=pop(:,1:L);  
EF=zeros(size(pop,1),3*L); ES1=EF; ES2=EF; EQ=EF;
DF=EF; DS1=EF; DS2=EF;
UF=EF; US=EF;
DQF=EF; DQS=EF; UQ=EF;
a=[1:L];


% DO THE VoC to define the size
SP=2*(8*L+4*m*L);   % offset for VoC
LL=2*L;
for i=1:m
    EF(:,a+LL)=EF(:,a+LL)+pop(:,SP+L+a+3*(i-1)*L); 
    ES1(:,a+LL)=ES1(:,a+LL)+pop(:,SP+2*L+a+3*(i-1)*L); 
    ES2(:,a+LL)=ES2(:,a+LL)+pop(:,SP+3*L+a+3*(i-1)*L);
    EQ(:,a+LL)=EQ(:,a+LL)+pop(:,SP+6*L+a+3*m*L+(i-1)*L);
end
DF(:,a+LL)=pop(:,SP+1*L+a+3*m*L); DS1(:,a+LL)=pop(:,SP+2*L+a+3*m*L); DS2(:,a+L)=pop(:,SP+3*L+a+3*m*L); 
UF(:,a+LL)=pop(:,SP+4*L+a+3*m*L); US(:,a+LL)=pop(:,SP+5*L+a+3*m*L); 
DQF(:,a+LL)=pop(:,SP+6*L+a+4*m*L); DQS(:,a+LL)=pop(:,SP+7*L+a+4*m*L); UQ(:,a+L)=pop(:,SP+8*L+a+4*m*L);

% DO THE B117 VARIANT
SP=(8*L+4*m*L);   % offset for B117
for i=1:m
    EF(:,a+L)=EF(:,a+L)+pop(:,SP+L+a+3*(i-1)*L); 
    ES1(:,a+L)=ES1(:,a+L)+pop(:,SP+2*L+a+3*(i-1)*L); 
    ES2(:,a+L)=ES2(:,a+L)+pop(:,SP+3*L+a+3*(i-1)*L);
    EQ(:,a+L)=EQ(:,a+L)+pop(:,SP+6*L+a+3*m*L+(i-1)*L);
end
DF(:,a+L)=pop(:,SP+1*L+a+3*m*L); DS1(:,a+L)=pop(:,SP+2*L+a+3*m*L); DS2(:,a+L)=pop(:,SP+3*L+a+3*m*L); 
UF(:,a+L)=pop(:,SP+4*L+a+3*m*L); US(:,a+L)=pop(:,SP+5*L+a+3*m*L); 
DQF(:,a+L)=pop(:,SP+6*L+a+4*m*L); DQS(:,a+L)=pop(:,SP+7*L+a+4*m*L); UQ(:,a+L)=pop(:,SP+8*L+a+4*m*L);

%THE DO ORIGINAL VARIANT
for i=1:m
    EF(:,a)=EF(:,a)+pop(:,L+a+3*(i-1)*L); 
    ES1(:,a)=ES1(:,a)+pop(:,2*L+a+3*(i-1)*L); 
    ES2(:,a)=ES2(:,a)+pop(:,3*L+a+3*(i-1)*L);
    EQ(:,a)=EQ(:,a)+pop(:,6*L+a+3*m*L+(i-1)*L);
end
DF(:,a)=pop(:,1*L+a+3*m*L); DS1(:,a)=pop(:,2*L+a+3*m*L); DS2(:,a)=pop(:,3*L+a+3*m*L); 
UF(:,a)=pop(:,4*L+a+3*m*L); US(:,a)=pop(:,5*L+a+3*m*L); 
DQF(:,a)=pop(:,6*L+a+4*m*L); DQS(:,a)=pop(:,7*L+a+4*m*L); UQ(:,a)=pop(:,8*L+a+4*m*L);

%E=EF+ES1+ES2+EQ;  D=DF+DS1+DS2+DQF+DQS; U=UF+US+UQ;


E=EF+ES1+ES2+EQ;  D=DF+DS1+DS2+DQF+DQS; U=UF+US+UQ;

%Ratios is a list of all those that can potentially be infected:  S, SV1, SV2, RV1, RV2 and R.
SPV=3*(8*L+4*m*L)+L;
Ratios=zeros(size(pop,1),21,6);
Ratios(:,a,1)=pop(:,a);  
Ratios(:,a,2)=pop(:,SPV+2*L+[1:L]);
Ratios(:,a,3)=pop(:,SPV+3*L+[1:L]);
Ratios(:,a,4)=pop(:,a+SPV+4*L);
Ratios(:,a,5)=pop(:,a+SPV+5*L); 
Ratios(:,a,6)=pop(:,a+SPV+6*L); 
%Ratio(a+L)=(pop(:,a))./(pop(:,a)+pop(:,SPV+2*L+[1:L])+pop(:,SPV+3*L+[1:L]) + pop(:,a+SPV+4*L) + pop(:,a+SPV+5*L) + pop(:,a+SPV+6*L));  


% Find all new detectable / symptomatic infections
nD=0*EF;
a=[1:L];
nD(:,a)=(ones(size(pop,1),1)*d(a)').*(alpha*m*(pop(:,L+a+3*(m-1)*L)+pop(:,2*L+a+3*(m-1)*L)+pop(:,3*L+a+3*(m-1)*L)+pop(:,6*L+a+3*m*L+(m-1)*L)));
SP=8*L+4*m*L;
nD(:,a+L)=(ones(size(pop,1),1)*d(a)').*(alpha*m*(pop(:,SP+L+a+3*(m-1)*L)+pop(:,SP+2*L+a+3*(m-1)*L)+pop(:,SP+3*L+a+3*(m-1)*L)+pop(:,SP+6*L+a+3*m*L+(m-1)*L)));
SP=2*(8*L+4*m*L);
op=ones(size(pop,1),1); OD=(op*Od');
nd=(squeeze(Ratios(:,:,1)).*OD + squeeze(Ratios(:,:,2)).*OD*Symptom_Reduction(1) + squeeze(Ratios(:,:,3)).*OD*Symptom_Reduction(2) ...
    + squeeze(Ratios(:,:,4)).*OD*Symptom_Reduction(3) + squeeze(Ratios(:,:,5)).*OD*Symptom_Reduction(4) + squeeze(Ratios(:,:,6)).*OD*Symptom_Reduction(5))./(squeeze(sum(Ratios,3)));
nD(:,a+LL)=(nd(:,a)).*(alpha*m*(pop(:,SP+L+a+3*(m-1)*L)+pop(:,SP+2*L+a+3*(m-1)*L)+pop(:,SP+3*L+a+3*(m-1)*L)+pop(:,SP+6*L+a+3*m*L+(m-1)*L)));


FinalState=pop(end,:);

end

% Calculates the differential rates used in the integration.
function dPop=Diff_3_4(t, pop, parameter)

L=parameter(1);
m=parameter(2);

s=2; 
l=L*L;  M_from_toH=reshape(parameter(s+[1:l]),L,L);    s=s+l;
l=L*L;  M_from_toO=reshape(parameter(s+[1:l]),L,L);    s=s+l;
l=L;    sigma = parameter(s+[1:l]);           s=s+l;
l=L;    d = parameter(s+[1:l]);               s=s+l;
l=L;    Od = parameter(s+[1:l]);              s=s+l;
l=L;    tau = parameter(s+[1:l]);             s=s+l;

l=1;    alpha=parameter(s+[1:l]);                 s=s+l;
l=1;    gamma=parameter(s+[1:l]);                 s=s+l; 
l=L;    N=parameter(s+[1:l]);                     s=s+l;
l=L;    hhq=parameter(s+[1:l]);                   s=s+l;
l=2;    nV_Beta=parameter(s+[1:l]);               s=s+l;
l=2;    nV_Speed=parameter(s+[1:l]);              s=s+l;

l=L;    V1=parameter(s+[1:l]);                    s=s+l;
l=L;    V2=parameter(s+[1:l]);                    s=s+l;
l=L;    VET1=parameter(s+[1:l]);                  s=s+l;
l=L;    VET2=parameter(s+[1:l]);                  s=s+l;

l=5;    VoC_Sus=parameter(s+[1:l]);               s=s+l;
l=5;    Transmission_Reduction=parameter(s+[1:l]);               s=s+l;
l=5;    Symptom_Reduction=parameter(s+[1:l]);               s=s+l;
l=L;    RatioPf=parameter(s+[1:l]);               s=s+l;


% THIS IS ALL FOR THE WT VARIANT

%DS1 are infected by a detected and cannot cause lock-down, DS2 are
%infected by undected and can. In the write-up these are D^SD and D^SU
S=pop(1:L)'; 
SPV=3*(8*L+4*m*L)+L;
SV1=pop(SPV+2*L+[1:L])';  SV2=pop(SPV+3*L+[1:L])';  % these are vaccinated but still susceptible to B117 or WT.
RV1=pop(SPV+4*L+[1:L])';  RV2=pop(SPV+5*L+[1:L])';  % these are vaccinated but susceptible to VoC.
REC=pop(SPV+6*L+[1:L])';  % recovered from B117 or WT, but still Sus to VoC.

SALL=S+SV1+SV2;

EF=zeros(m,L); ES1=zeros(m,L); ES2=zeros(m,L); EQ=zeros(m,L);
for i=1:m
    EF(i,:)=pop(L+[1:L]+3*(i-1)*L)';  ES1(i,:)=pop(2*L+[1:L]+3*(i-1)*L)'; ES2(i,:)=pop(3*L+[1:L]+3*(i-1)*L)';
end
DF=pop(1*L+[1:L]+3*m*L)'; DS1=pop(2*L+[1:L]+3*m*L)'; DS2=pop(3*L+[1:L]+3*m*L)'; 
UF=pop(4*L+[1:L]+3*m*L)'; US=pop(5*L+[1:L]+3*m*L)'; 
for i=1:m
    EQ(i,:)=pop(6*L+[1:L]+3*m*L+(i-1)*L)';
end
DQF=pop(6*L+[1:L]+4*m*L)'; DQS=pop(7*L+[1:L]+4*m*L)'; UQ=pop(8*L+[1:L]+4*m*L)'; 

% if hhq is in place - put individals straight into Q classes.

dPop=zeros(length(pop),1);  RecPos=SPV+6*L;

IF=DF + tau.*UF;   IS=(DS1+DS2) + tau.*US; 

a=[1:L];
Trans_Scaling=(1*S + Transmission_Reduction(1)*SV1 + Transmission_Reduction(2)*SV2)./(S+SV1+SV2);
InfF=Trans_Scaling(a).*(sigma(a).*((IF+ IS)*M_from_toO(:,a)));
InfS1=Trans_Scaling(a).*(sigma(a).*((DF)*M_from_toH(:,a)));
InfS2=Trans_Scaling(a).*(sigma(a).*((tau.*UF)*M_from_toH(:,a)));
InfSQ=Trans_Scaling(a).*(sigma(a).*((DQF)*M_from_toH(:,a)));

s=0;
%dS
dPop(a+s)= - (InfF + InfS1 + InfS2 + InfSQ).*S(a)./N(a);             s=s+L;
dPop(a+SPV+2*L)= - (InfF + InfS1 + InfS2 + InfSQ).*SV1(a)./N(a);       
dPop(a+SPV+3*L)= - (InfF + InfS1 + InfS2 + InfSQ).*SV2(a)./N(a);       
    
    
    %dEF   dES1  and dES2 (first ones)
    dPop(a+s)=InfF.*SALL(a)./N(a)  - m*alpha*EF(1,a);                s=s+L;
    dPop(a+s)=InfS1.*SALL(a)./N(a) - m*alpha*ES1(1,a);               s=s+L;
    dPop(a+s)=InfS2.*SALL(a)./N(a) - m*alpha*ES2(1,a);               s=s+L;
    
    %dEF  dES1  & dES2 (subsequent ones)
    for i=2:m
        dPop(a+s)=m*alpha*EF(i-1,a) - m*alpha*EF(i,a);      s=s+L;
        dPop(a+s)=m*alpha*ES1(i-1,a) - m*alpha*ES1(i,a);    s=s+L;
        dPop(a+s)=m*alpha*ES2(i-1,a) - m*alpha*ES2(i,a);    s=s+L;
    end
    
    %dDF  dDS1  and dDS2
    dPop(a+s)=d(a) .* (1-hhq(a)) .* alpha .* EF(m,a) * m - gamma*DF(a);       s=s+L;
    dPop(a+s)=d(a) .* alpha .* ES1(m,a) * m - gamma*DS1(a);                   s=s+L;
    dPop(a+s)=d(a) .* (1-hhq(a)) .* alpha .* ES2(m,a) * m - gamma*DS2(a);     s=s+L;
    
    %dUF and dUS
    dPop(a+s)=(1-d(a)) .* alpha .* EF(m,a) * m - gamma*UF(a);                   s=s+L;
    dPop(a+s)=(1-d(a)) .* alpha .* (ES1(m,a) + ES2(m,a)) * m - gamma*US(a);     s=s+L;
    
    %dEQ
    dPop(a+s) = InfSQ - alpha .* EQ(1,a) * m;       s=s+L;
    for i=2:m
        dPop(a+s) = alpha .* EQ(i-1,a) * m - alpha .* EQ(i,a) * m;       s=s+L;
    end
    
    %dDQF
    dPop(a+s) = d(a) .* hhq(a) .* alpha .* EF(m,a) * m - gamma .* DQF(a);       s=s+L;

    %dDQS
    dPop(a+s) = d(a) .* hhq(a) .* alpha .* ES2(m,a) * m  +  d(a) .* alpha .* EQ(m,a) * m - gamma*DQS(a);       s=s+L;
    
    %dUQ
    dPop(a+s) = (1-d(a)) .* alpha .* EQ(m,a) * m - gamma*UQ(a);       s=s+L;
    
    dPop(a+RecPos) = VoC_Sus(5)*gamma*(DF(a)+DS1(a)+DS2(a)+UF(a)+US(a)+DQF(a)+DQS(a)+UQ(a));
   
% REPEAT FOR THE B117 VARIANT

SP=(8*L+4*m*L);  % 8 states + 4 with m classes. [S class is not included].

EF=zeros(m,L); ES1=zeros(m,L); ES2=zeros(m,L); EQ=zeros(m,L);
for i=1:m
    EF(i,:)=pop(SP+L+[1:L]+3*(i-1)*L)';  
    ES1(i,:)=pop(SP+2*L+[1:L]+3*(i-1)*L)'; 
    ES2(i,:)=pop(SP+3*L+[1:L]+3*(i-1)*L)';
    EQ(i,:)=pop(SP+6*L+[1:L]+3*m*L+(i-1)*L)';
end
DF=pop(SP+1*L+[1:L]+3*m*L)'; DS1=pop(SP+2*L+[1:L]+3*m*L)'; DS2=pop(SP+3*L+[1:L]+3*m*L)'; 
UF=pop(SP+4*L+[1:L]+3*m*L)'; US=pop(SP+5*L+[1:L]+3*m*L)'; 
DQF=pop(SP+6*L+[1:L]+4*m*L)'; DQS=pop(SP+7*L+[1:L]+4*m*L)'; UQ=pop(SP+8*L+[1:L]+4*m*L)'; 

IF=DF + tau.*UF;   IS=(DS1+DS2) + tau.*US;

a=[1:L]; nBeta=nV_Beta(1)*nV_Speed(1);
Trans_Scaling=(1*S + Transmission_Reduction(1)*SV1 + Transmission_Reduction(2)*SV2)./(S+SV1+SV2);
InfF=nBeta*Trans_Scaling.*(sigma(a).*((IF+ IS)*M_from_toO(:,a)));
InfS1=nBeta*Trans_Scaling.*(sigma(a).*((DF)*M_from_toH(:,a)));
InfS2=nBeta*Trans_Scaling.*(sigma(a).*((tau.*UF)*M_from_toH(:,a)));
InfSQ=nBeta*Trans_Scaling.*(sigma(a).*((DQF)*M_from_toH(:,a)));

%dS
dPop(a)= dPop(a) - ((InfF + InfS1 +InfS2 + InfSQ).*S(a)./N(a))';       %
dPop(a+SPV+2*L)= dPop(a+SPV+2*L) - ((InfF + InfS1 + InfS2 + InfSQ).*SV1(a)./N(a))';       
dPop(a+SPV+3*L)= dPop(a+SPV+3*L) - ((InfF + InfS1 + InfS2 + InfSQ).*SV2(a)./N(a))';       
    
nalpha=alpha*nV_Speed(1);
%dEF   dES1  and dES2 (first ones)
dPop(a+s)=InfF.*SALL(a)./N(a)  - m*nalpha*EF(1,a);                s=s+L;
dPop(a+s)=InfS1.*SALL(a)./N(a) - m*nalpha*ES1(1,a);               s=s+L;
dPop(a+s)=InfS2.*SALL(a)./N(a) - m*nalpha*ES2(1,a);               s=s+L;

%dEF  dES1  & dES2 (subsequent ones)
for i=2:m
    dPop(a+s)=m*nalpha*EF(i-1,a) - m*nalpha*EF(i,a);      s=s+L;
    dPop(a+s)=m*nalpha*ES1(i-1,a) - m*nalpha*ES1(i,a);    s=s+L;
    dPop(a+s)=m*nalpha*ES2(i-1,a) - m*nalpha*ES2(i,a);    s=s+L;
end

ngamma=gamma*nV_Speed(1);
%dDF  dDS1  and dDS2
dPop(a+s)=d(a) .* (1-hhq(a)) .* nalpha .* EF(m,a) * m - ngamma*DF(a);       s=s+L;
dPop(a+s)=d(a) .* nalpha .* ES1(m,a) * m - ngamma*DS1(a);                   s=s+L;
dPop(a+s)=d(a) .* (1-hhq(a)) .* nalpha .* ES2(m,a) * m - ngamma*DS2(a);     s=s+L;

%dUF and dUS
dPop(a+s)=(1-d(a)) .* nalpha .* EF(m,a) * m - ngamma*UF(a);                   s=s+L;
dPop(a+s)=(1-d(a)) .* nalpha .* (ES1(m,a) + ES2(m,a)) * m - ngamma*US(a);     s=s+L;

%dEQ
dPop(a+s) = InfSQ - nalpha .* EQ(1,a) * m;       s=s+L;
for i=2:m
    dPop(a+s) = nalpha .* EQ(i-1,a) * m - nalpha .* EQ(i,a) * m;       s=s+L;
end

%dDQF
dPop(a+s) = d(a) .* hhq(a) .* nalpha .* EF(m,a) * m - ngamma .* DQF(a);       s=s+L;

%dDQS
dPop(a+s) = d(a) .* hhq(a) .* nalpha .* ES2(m,a) * m  +  d(a) .* alpha .* EQ(m,a) * m - ngamma*DQS(a);       s=s+L;

%dUQ
dPop(a+s) = (1-d(a)) .* nalpha .* EQ(m,a) * m - ngamma*UQ(a);       s=s+L;

% Recovereds
dPop(a+RecPos) = dPop(a+RecPos) + VoC_Sus(5)*ngamma*(DF(a)+DS1(a)+DS2(a)+UF(a)+US(a)+DQF(a)+DQS(a)+UQ(a))';


% AND AGAIN FOR THE VoC

SP=2*(8*L+4*m*L);  % Take off an L as no longer have S's

EF=zeros(m,L); ES1=zeros(m,L); ES2=zeros(m,L); EQ=zeros(m,L);
for i=1:m
    EF(i,:)=pop(SP+L+[1:L]+3*(i-1)*L)';  
    ES1(i,:)=pop(SP+2*L+[1:L]+3*(i-1)*L)'; 
    ES2(i,:)=pop(SP+3*L+[1:L]+3*(i-1)*L)';
    EQ(i,:)=pop(SP+6*L+[1:L]+3*m*L+(i-1)*L)';
end
DF=pop(SP+1*L+[1:L]+3*m*L)'; DS1=pop(SP+2*L+[1:L]+3*m*L)'; DS2=pop(SP+3*L+[1:L]+3*m*L)'; 
UF=pop(SP+4*L+[1:L]+3*m*L)'; US=pop(SP+5*L+[1:L]+3*m*L)'; 
DQF=pop(SP+6*L+[1:L]+4*m*L)'; DQS=pop(SP+7*L+[1:L]+4*m*L)'; UQ=pop(SP+8*L+[1:L]+4*m*L)'; 

IF=DF + tau.*UF;   IS=(DS1+DS2) + tau.*US;

a=[1:L]; nBeta=nV_Beta(2)*nV_Speed(2);
Trans_Scaling=(1*S + Transmission_Reduction(1)*SV1 + Transmission_Reduction(2)*SV2 + Transmission_Reduction(3)*RV1 + Transmission_Reduction(4)*RV2 + Transmission_Reduction(5)*REC)./(S+SV1+SV2+RV1+RV2+REC);
InfF=nBeta*Trans_Scaling.*(sigma(a).*((IF+ IS)*M_from_toO(:,a)));
InfS1=nBeta*Trans_Scaling.*(sigma(a).*((DF)*M_from_toH(:,a)));
InfS2=nBeta*Trans_Scaling.*(sigma(a).*((tau.*UF)*M_from_toH(:,a)));
InfSQ=nBeta*Trans_Scaling.*(sigma(a).*((DQF)*M_from_toH(:,a)));

%dS
dPop(a)= dPop(a) - ((InfF + InfS1 +InfS2 + InfSQ).*S(a)./N(a))';       %
dPop(a+SPV+2*L)= dPop(a+SPV+2*L) - ((InfF + InfS1 + InfS2 + InfSQ).*SV1(a)./N(a))';       
dPop(a+SPV+3*L)= dPop(a+SPV+3*L) - ((InfF + InfS1 + InfS2 + InfSQ).*SV2(a)./N(a))';

% NEW !! and impacting vaccine derived immunity and natural immunity.
SPV=3*(8*L+4*m*L)+L;
dPop(a+SPV+4*L)= dPop(a+SPV+4*L) - ((InfF + InfS1 + InfS2 + InfSQ)'.*pop(a+SPV+4*L)./N(a)');       
dPop(a+SPV+5*L)= dPop(a+SPV+5*L) - ((InfF + InfS1 + InfS2 + InfSQ)'.*pop(a+SPV+5*L)./N(a)');
dPop(a+SPV+6*L)= dPop(a+SPV+6*L) - ((InfF + InfS1 + InfS2 + InfSQ)'.*pop(a+SPV+6*L)./N(a)');

oSALL=SALL;
SALL(a)=SALL(a) + pop(a+SPV+4*L)' + pop(a+SPV+5*L)' + pop(a+SPV+6*L)';

nd=(S.*Od + SV1.*Od*Symptom_Reduction(1) + SV2.*Od*Symptom_Reduction(2) + RV1.*Od*Symptom_Reduction(3) + RV2.*Od*Symptom_Reduction(4) + REC.*Od*Symptom_Reduction(5))./(S+SV1+SV2+RV1+RV2+REC);

nalpha=alpha*nV_Speed(2);
%dEF   dES1  and dES2 (first ones)
dPop(a+s)=InfF.*SALL(a)./N(a)  - m*nalpha*EF(1,a);                s=s+L;
dPop(a+s)=InfS1.*SALL(a)./N(a) - m*nalpha*ES1(1,a);               s=s+L;
dPop(a+s)=InfS2.*SALL(a)./N(a) - m*nalpha*ES2(1,a);               s=s+L;

%dEF  dES1  & dES2 (subsequent ones)
for i=2:m
    dPop(a+s)=m*nalpha*EF(i-1,a) - m*nalpha*EF(i,a);      s=s+L;
    dPop(a+s)=m*nalpha*ES1(i-1,a) - m*nalpha*ES1(i,a);    s=s+L;
    dPop(a+s)=m*nalpha*ES2(i-1,a) - m*nalpha*ES2(i,a);    s=s+L;
end

ngamma=gamma*nV_Speed(2);
%dDF  dDS1  and dDS2
dPop(a+s)=nd(a) .* (1-hhq(a)) .* nalpha .* EF(m,a) * m - ngamma*DF(a);       s=s+L;
dPop(a+s)=nd(a) .* nalpha .* ES1(m,a) * m - ngamma*DS1(a);                   s=s+L;
dPop(a+s)=nd(a) .* (1-hhq(a)) .* nalpha .* ES2(m,a) * m - ngamma*DS2(a);     s=s+L;

%dUF and dUS
dPop(a+s)=(1-nd(a)) .* nalpha .* EF(m,a) * m - ngamma*UF(a);                   s=s+L;
dPop(a+s)=(1-nd(a)) .* nalpha .* (ES1(m,a) + ES2(m,a)) * m - ngamma*US(a);     s=s+L;

%dEQ
dPop(a+s) = InfSQ - nalpha .* EQ(1,a) * m;       s=s+L;
for i=2:m
    dPop(a+s) = nalpha .* EQ(i-1,a) * m - nalpha .* EQ(i,a) * m;       s=s+L;
end

%dDQF
dPop(a+s) = nd(a) .* hhq(a) .* nalpha .* EF(m,a) * m - ngamma .* DQF(a);       s=s+L;

%dDQS
dPop(a+s) = nd(a) .* hhq(a) .* nalpha .* ES2(m,a) * m  +  d(a) .* alpha .* EQ(m,a) * m - ngamma*DQS(a);       s=s+L;

%dUQ
dPop(a+s) = (1-nd(a)) .* nalpha .* EQ(m,a) * m - ngamma*UQ(a);       s=s+L;

% Recovereds are taken out of the picture.
% dPop(a+RecPos) = dPop(a+RecPos) + gamma*(DF(a)+DS1(a)+DS2(a)+UF(a)+US(a)+DQF(a)+DQS(a)+UQ(a));


% NOW DO VACCINATION
SPV=3*(8*L+4*m*L)+L; s=SPV;
Vacc1=pop(SPV+a)'; 
Vacc2=pop(SPV+L+a)'; 
SET1=1-VET1;  SET2=1-VET2;

%%% Subtract vaccination from S.
dPop(a) = dPop(a) - (V1.*S./(1+N-Vacc1-Vacc2))' ;

%dVacc1 -- everyone that has been vaccinated
dPop(a+s) = V1 - V2;          s=s+L;
%dVacc2 -- everyone that has been vaccinated
dPop(a+s) = V2;          s=s+L;

%dSVacc1 -- vaccinated once but still susceptible
dPop(a+s) = dPop(a+s) +  (V1.*S.*SET1./(1+N-Vacc1-Vacc2))'  - (V2.*SV1./(1+Vacc1))' ;        s=s+L;
%dSVacc2 -- vaccinated twice but still susceptible
dPop(a+s) = dPop(a+s) +  (V2.*(SV1./(1+Vacc1)).*(SET2./SET1))';          s=s+L;

%dSVacc1 -- vaccinated once and now "resistant"
% Calculate this as a faction VoC_Sus(1) of those that would be resistant
% then take off all 2nd doses
VocS1=RatioPf(a)*VoC_Sus(1)+(1-RatioPf(a))*VoC_Sus(3);
dPop(a+s) = dPop(a+s) +  (VocS1.*V1.*VET1.*S./(1+N-Vacc1-Vacc2))'  - (V2.*RV1./(1+Vacc1))' ;         s=s+L;
%dSVacc2 -- vaccinated twice and now "resistant"
OUTfromSV1=(V2.*(SV1./(1+Vacc1)).*(1-(SET2./SET1)))';  % those leaving SV1 that do not go into SV2
OUTfromRV1=(V2.*RV1./(1+Vacc1))';    % those leaving RV1
VocS2=RatioPf(a)*VoC_Sus(2)+(1-RatioPf(a))*VoC_Sus(4);
dPop(a+s) = dPop(a+s) +  (OUTfromSV1+OUTfromRV1).*((VocS2.*VET2)./((VET2-VET1)+VocS1.*VET1))';    s=s+L;

% and finally recovereds which are done above !!
end


