function [nDC, nDCH, nDCD, nHospital, inHospital, nICU, inICU, nDeaths, PD_Lockdown, INF]=Simulate_One_Region(Region, TAU, ALPHA, INC_P, S_Scale, Factor, h_factor, i_factor, d_factor, h_stretch, ...
    i_stretch, Lag, Start_Date,  WALES_FLAG, ComplianceT, Run_stop, nV_Beta, nV_Speed, Import_Day, Import_Level, ...
    V1, V2, VET1, VET2, Transmission_Reduction, RC, RCH, RCD, VoC_Sus, General_Reduction, ...
    RatioPf, Detection, Susceptibility, gamma, Seasonality, XTrans, XTransDate)


ODEs=@VoC_Vacc_ODEs2; % this function could be compiled to mex code to make everything run faster.
%ODEs=@VoC_Vacc_ODEs2_mex; % something like this

ODE_to_Obs=@ODE_to_Observables; % this function could be compiled to mex code to make everything run faster.
%ODE_to_Obs=@ODE_to_Observables_mex;  % something like this


Christmas=0;

if Run_stop(1)>1000
    Run_stop=Run_stop+1-datenum(2020,1,1);
end

S=size(Run_stop);
if S(2)==1
    Run_stop=Run_stop';
    ComplianceT=ComplianceT';
end

if length(h_factor)==1
    h_factor(2)=h_factor(1); i_factor(2)=i_factor(1); d_factor(2)=d_factor(1);
    h_stretch(2)=h_stretch(1); i_stretch(2)=i_stretch(1);
end


if length(VoC_Sus)==3
    VoC_Sus=VoC_Sus([1 2 1 2 3]);
end


load Other_Data

UK_PP=UK_PP'*sum(Region_PP(Region,:))/sum(UK_PP);

UK_from_toH = UK_from_toH .* ( ones(21,1) * (Region_PP(Region,:)./UK_PP));
UK_from_toW = UK_from_toW .* ( ones(21,1) * (Region_PP(Region,:)./UK_PP));
UK_from_toS = UK_from_toS .* ( ones(21,1) * (Region_PP(Region,:)./UK_PP));
UK_from_toO = UK_from_toO .* ( ones(21,1) * (Region_PP(Region,:)./UK_PP));

%Sort_Data_to_Match

Names={'ENGLAND','East of England','London','Midlands','North East and Yorkshire','North West','South East','South West','Wales','Scotland','Northern Ireland'};

tau=TAU;
%fprintf(1,'TEST_RUNS gamma=%g\n',gamma);

Early_R0=3.8;
Early_gamma=gamma/1.5;
Z=1.3;
a=INC_P;

%[T2,S,E,D,nD2,U , R0, Da, FinalState]=nV_Vacc_ODEs_mex(UK_from_toH*0 , gamma*(UK_from_toH + UK_from_toW + UK_from_toS + UK_from_toO), a, Early_gamma, Susceptibility, Detection, TAU, nV_Beta(1), nV_Speed(1), 0, Region_PP(Region,:)' , V1(1,:), V2(1,:), VET1, VET2, 2, -1);
tmp=[1 0.924 0.925 0.944 0.942 0.947 0.925 0.874 0.912 0.928 0.971];
Susceptibility2=S_Scale*Susceptibility*tmp(Region);
Susceptibility3=S_Scale*Susceptibility;

if length(Seasonality)==2
    TT=Seasonality(2); 
else
    TT=round(datenum(now)+1-datenum(2020,1,1));
end
T=[1:max(Run_stop)];
Extra_Seasonality=(1-0.5*Seasonality(1)-0.5*Seasonality(1)*cos((T-593)*2*pi/365)); % RAW
Extra_Seasonality(1:TT)=1; 

if length(Import_Day)==1
    Import_Day(2)=round(datenum(now)+1-datenum(2020,1,1));
end

Symptom_Reduction=General_Reduction(1,:);

[T2,S,E,D,nD2,U , Ratio, Da, FinalState]=ODEs(UK_from_toH*0 , gamma*(UK_from_toH + UK_from_toW + UK_from_toS + UK_from_toO), a, Early_gamma, Susceptibility, Detection, Detection, TAU, nV_Beta, nV_Speed, 0, Region_PP(Region,:)' , V1(1,:), V2(1,:), VET1, VET2, VoC_Sus, Transmission_Reduction, Symptom_Reduction, RatioPf,  2, -1);


%% EARLY SET-UP
%Run Hot to 12th March (day 71)
FinalState(22:end)=FinalState(22:end)*Factor;%/sum(FinalState(22:end));
[T,S,E,D,nD,U , Ratio, Da, FinalState]=ODEs(gamma*UK_from_toH*Z , gamma*(UK_from_toW + UK_from_toS + UK_from_toO), a, Early_gamma, Susceptibility2, Detection, Detection, tau, nV_Beta, nV_Speed, 0, Region_PP(Region,:)', V1(1,:), V2(1,:),VET1, VET2, VoC_Sus, Transmission_Reduction, Symptom_Reduction, RatioPf,71-Start_Date, [FinalState]);
T=T+Start_Date;

%Self Isolation for 4 days
[t,s,e,d,nd,u , ratio, Da, FinalState]=ODEs(gamma*UK_from_toH*Z, gamma*(UK_from_toW + UK_from_toS + UK_from_toO), a, gamma, Susceptibility3, Detection, Detection, tau, nV_Beta, nV_Speed, 0, Region_PP(Region,:)' , V1(1,:), V2(1,:),VET1, VET2, VoC_Sus, Transmission_Reduction, Symptom_Reduction, RatioPf, 4, [FinalState]);
T=[T; T(end)+t(2:end)]; S=[S; s(2:end,:)]; E=[E; e(2:end,:)]; D=[D; d(2:end,:)]; nD=[nD; nd(2:end,:)]; U=[U; u(2:end,:)]; Ratio=[Ratio; ratio(2:end,:,:)];

%Work from Home + Self Isolation for 4 days
[new_UK_from_toH, new_UK_from_toW, new_UK_from_toS, new_UK_from_toO] = Return_New_Matrices(0, 0.3, 0, ComplianceT(:,1), UK_from_toH, UK_from_toW, UK_from_toS, UK_from_toO);
[t,s,e,d,nd,u , ratio, Da, FinalState]=ODEs(gamma*new_UK_from_toH*Z, gamma*(new_UK_from_toW + new_UK_from_toS + new_UK_from_toO), a, gamma, Susceptibility3, Detection, Detection, tau, nV_Beta, nV_Speed, 0, Region_PP(Region,:)', V1(1,:), V2(1,:),VET1, VET2, VoC_Sus, Transmission_Reduction, Symptom_Reduction, RatioPf, 4, [FinalState]);
T=[T; T(end)+t(2:end)]; S=[S; s(2:end,:)]; E=[E; e(2:end,:)]; D=[D; d(2:end,:)]; nD=[nD; nd(2:end,:)]; U=[U; u(2:end,:)]; Ratio=[Ratio; ratio(2:end,:,:)];

%Some x (Social distancing + HHQ + School Closures + Work from Home + Self Isolatio)n for 3 days
[new_UK_from_toH, new_UK_from_toW, new_UK_from_toS, new_UK_from_toO] = Return_New_Matrices(0.8, 0.5, 0.3, ComplianceT(:,1), UK_from_toH, UK_from_toW, UK_from_toS, UK_from_toO);
[t,s,e,d,nd,u , ratio, Da, FinalState]=ODEs(gamma*new_UK_from_toH*Z, gamma*(new_UK_from_toW + new_UK_from_toS + new_UK_from_toO), a, gamma, Susceptibility3, Detection, Detection, tau, nV_Beta, nV_Speed, 0.3*ComplianceT(1,1), Region_PP(Region,:)', V1(1,:), V2(1,:), VET1, VET2, VoC_Sus, Transmission_Reduction, Symptom_Reduction, RatioPf, 3, [FinalState]);
T=[T; T(end)+t(2:end)]; S=[S; s(2:end,:)]; E=[E; e(2:end,:)]; D=[D; d(2:end,:)]; nD=[nD; nd(2:end,:)]; U=[U; u(2:end,:)]; Ratio=[Ratio; ratio(2:end,:,:)];


%% THEN RUN SENARIOS

PD_Lockdown = 0;

Flag=0; fixedaC=0;
UF=1; DF=1;

Multi_Factor=1.0;

for W=1:length(Run_stop)
     
    mmC=max(max(abs(ComplianceT)));
    if mmC>0
        if length(ComplianceT(:,W))==1
            aC=ones(21,1)*abs(ComplianceT(:,W));
        end
        if length(ComplianceT(:,W))==4
            aC=abs(ComplianceT([1 2 2 2 3 3 3 3 3 3 3 3 3 3 4 4 4 4 4 4 4],W));
        end
        if length(ComplianceT(:,W))==21
            aC=reshape(abs(ComplianceT(:,W)),21,1);
        end
        PD_Lockdown=PD_Lockdown+(Run_stop(W)-T(end)).*(Region_PP(Region,:)*aC)./mmC;
    end
    
    % INCLUDE THIS FOR RE-INTRODUCTIONS !!
    if sum(FinalState(22:(end-84)))<1 & W>1
        FinalState(22:(end-84))=abs(FinalState(22:(end-84)))/abs(sum(FinalState(22:(end-84))));
    end
    
    % Add in B117
    if T(end)<Import_Day(1) & Run_stop(W)>=Import_Day(1)
        SP=9*21+4*3*21;
        FinalState(SP+[1:(SP-21)])=Import_Level(1)*FinalState(22:SP);
        FinalState(22:SP)=(1-Import_Level(1))*FinalState(22:SP);
    end
    % Add in VoC
    if T(end)<Import_Day(2) & Run_stop(W)>=Import_Day(2)
        SP=(9*21+4*3*21); SP2=SP-21;  tmp=FinalState(22:SP)+FinalState(SP+[1:SP2]);
        FinalState(SP+SP2+[1:SP2])=Import_Level(2)*tmp;
        FinalState(22:SP)=0*(1-Import_Level(2))*FinalState(22:SP);  % Get rid of pre-B117 wildetypes.
        FinalState(SP+[1:SP2])=(1-Import_Level(2))*FinalState(SP+[1:SP2]);
    end
    
    InSchoolFlag=0;
    if Region==10
        if T(end)>=230  % ie Mid-August put schools back in Scotland
            InSchoolFlag=1;
        end
        if T(end)>=292 && T(end)<=294
            InSchoolFlag=0;  % Unless its Half Term
        end
        if (T(end)>=355 && T(end)<=365)
            InSchoolFlag=0;   % ... or Christmas
        end
    else
        if T(end)>=244 % ie 1st Sept put schools back
            InSchoolFlag=1;
        end
        if T(end)>=299 && T(end)<=301
            InSchoolFlag=0;   % Unless its Half Term
        end
        if (T(end)>=355 && T(end)<=365)
            InSchoolFlag=0;   % ... or Christmas
        end
    end
    
    if Region==9 && T(end)>=347 && T(end)<=365
        InSchoolFlag=0;  % Weeks break in Wales
    end
    
    if Region<9 || Region==11   % In England
        if (T(end)>=355 && T(end)<=433-2)
            InSchoolFlag=0;   % ... EXTENDED CHRISTMAS BREAK & LOCKDOWN
        end
    else
        if (T(end)>=355 && T(end)<=440-2)
            InSchoolFlag=0;   % ... EXTENDED CHRISTMAS BREAK & LOCKDOWN
        end
    end
    
    
    if T(end)>=456 && T(end)<470 
        InSchoolFlag=0; % Easter '21!
    end
    
    if T(end)>=datenum(1,7,23)-2 && T(end)<datenum(1,9,6)-2 && Region<9
        InSchoolFlag=0; % Summer holidays '21 England!
    end
    if T(end)>=datenum(1,7,20)-2 && T(end)<datenum(1,9,2)-2 && Region==10
        InSchoolFlag=0; % Summer holidays '21 Wales!
    end
    if T(end)>=datenum(1,6,25)-2 && T(end)<datenum(1,8,16)-2 && Region==11
        InSchoolFlag=0; % Summer holidays '21 Scotland!
    end
    if T(end)>=datenum(1,7,1)-2 && T(end)<datenum(1,9,1)-2 && Region==12
        InSchoolFlag=0; % Summer holidays '21 N.Ireland!
    end
    
    
    if T(end)>=datenum(1,12,17)-2 && T(end)<datenum(2,1,4)-2
        InSchoolFlag=0; % Xmax holidays '21!
    end
    
    if T(end)>=datenum(2,4,8)-2 && T(end)<datenum(2,4,22)-2
        InSchoolFlag=0; % Easter holidays '22!
    end
    
    if T(end)>=datenum(2,7,25)-2 && T(end)<datenum(2,9,3)-2
        InSchoolFlag=0; % Summer holidays '22!
    end
     
    Extra_Beta = mean(Extra_Seasonality(T(end)+1:Run_stop(W))); 
    
    if T(end)>XTransDate
        nV_Beta(2)=nV_Beta(2)*XTrans;
        XTransDate=1e10;  % don't do this again.
    end
    
    oldaC=aC;
    
    if InSchoolFlag    % reduce the compliance between kids
        aC(2:4)=0.2*aC(2:4);
    end
    
    [new_UK_from_toH, new_UK_from_toW, new_UK_from_toS, new_UK_from_toO] = Return_New_Matrices(0.95, 0.8, 0.95, aC, UK_from_toH, UK_from_toW, UK_from_toS, UK_from_toO);
    mRC=mean(RC((T(end)+1):Run_stop(W),:),1);
    
    [t,s,e,d,nd,u , ratio, Da, FinalState]=ODEs(Extra_Beta*gamma*new_UK_from_toH*Z, Extra_Beta*gamma*(new_UK_from_toW + new_UK_from_toS + new_UK_from_toO), a, gamma, ...
        Susceptibility3, Detection.*mRC, Detection, tau, nV_Beta, nV_Speed, 0.8*abs(ComplianceT(:,W)), Region_PP(Region,:)' , mean(V1((T(end)+1):Run_stop(W),:),1), ...
        mean(V2((T(end)+1):Run_stop(W),:),1), VET1, VET2, VoC_Sus, Transmission_Reduction, Symptom_Reduction, RatioPf, Run_stop(W)-T(end), [FinalState]);
    T=[T; T(end)+t(2:end)]; S=[S; s(2:end,:)]; E=[E; e(2:end,:)]; D=[D; d(2:end,:)]; nD=[nD; nd(2:end,:)]; U=[U; u(2:end,:)]; Ratio=[Ratio; ratio(2:end,:,:)];
    
    if T(end)<Run_stop(W)
        fprintf(1,'Fails to complete, Region %d, Time %d=%s\n',Region, Run_stop(W)-7, datestr(Run_stop(W)-8+datenum(2020,1,1)));
    end
    
end

% Detectable Cases
clear nDC;
nDC(T+Lag,:)=nD;
Ratios(T+Lag,:,:)=Ratio; DRatios=Ratios; DRatios(DRatios==0)=1;


L=21; aa=[1:L]; LL=2*L;

% Less cases due to vaccination  Reduce for Hospitals & then again for
% Deaths.

nDCH(:,aa)=RCH(1:size(nDC,1),aa).*nDC(:,aa)./RC(1:size(nDC,1),aa);
nDCH(:,aa+L)=RCH(1:size(nDC,1),aa).*nDC(:,aa+L)./RC(1:size(nDC,1),aa);
Hospital_Reduction=General_Reduction(2,:);
tmpRatios=Ratios;
for qa=1:21
    tmpRatios(:,qa,2)=RatioPf(qa)*(Ratios(:,qa,2)+Ratios(:,qa,4));
    tmpRatios(:,qa,3)=RatioPf(qa)*(Ratios(:,qa,3)+Ratios(:,qa,5));
    tmpRatios(:,qa,4)=(1-RatioPf(qa))*(Ratios(:,qa,2)+Ratios(:,qa,4));
    tmpRatios(:,qa,5)=(1-RatioPf(qa))*(Ratios(:,qa,3)+Ratios(:,qa,5));
end
Scale=(squeeze(Ratios(:,:,1)) + squeeze(Ratios(:,:,2)).*Hospital_Reduction(1) + squeeze(Ratios(:,:,3)).*Hospital_Reduction(2) + squeeze(Ratios(:,:,4)).*Hospital_Reduction(3) + squeeze(Ratios(:,:,5)).*Hospital_Reduction(4) + squeeze(Ratios(:,:,6)).*Hospital_Reduction(5)) ...
    ./(squeeze(DRatios(:,:,1)) + squeeze(DRatios(:,:,2)).*Symptom_Reduction(1) + squeeze(DRatios(:,:,3)).*Symptom_Reduction(2) + squeeze(DRatios(:,:,4)).*Symptom_Reduction(3) + squeeze(DRatios(:,:,5)).*Symptom_Reduction(4) + squeeze(DRatios(:,:,6)).*Symptom_Reduction(5));
nDCH(:,aa+LL)=Scale(1:size(nDC,1),aa).*nDC(:,aa+LL);

nDCD(:,aa)=RCD(1:size(nDC,1),aa).*nDC(:,aa)./RC(1:size(nDC,1),aa);
nDCD(:,aa+L)=RCD(1:size(nDC,1),aa).*nDC(:,aa+L)./RC(1:size(nDC,1),aa);
Death_Reduction=General_Reduction(3,:);
Scale=(squeeze(Ratios(:,:,1)) + squeeze(Ratios(:,:,2)).*Death_Reduction(1) + squeeze(Ratios(:,:,3)).*Death_Reduction(2) + squeeze(Ratios(:,:,4)).*Death_Reduction(3) + squeeze(Ratios(:,:,5)).*Death_Reduction(4) + squeeze(Ratios(:,:,6)).*Death_Reduction(5)) ...
    ./(squeeze(DRatios(:,:,1)) + squeeze(DRatios(:,:,2)).*Symptom_Reduction(1) + squeeze(DRatios(:,:,3)).*Symptom_Reduction(2) + squeeze(DRatios(:,:,4)).*Symptom_Reduction(3) + squeeze(DRatios(:,:,5)).*Symptom_Reduction(4) + squeeze(DRatios(:,:,6)).*Symptom_Reduction(5));
nDCD(:,aa+LL)=Scale(1:size(nDC,1),aa).*nDC(:,aa+LL);

m2=min(length(h_factor),2);
m3=min(length(h_factor),3);

[NDC(:,aa), nHospital, inHospital, nICU, inICU, nDeaths, INF]=ODE_to_Obs(T, nDCH(:,aa), E(:,aa), Region, h_factor(1),i_factor(1),d_factor(1),h_stretch(1),i_stretch(1),a,...
    Assumed_Delay_Reporting_Deaths, Distribution_Hosp_Time', Distribution_HospICU_Time', Distribution_ICU_Time', Distribution_Symptoms_to_Hospital, Distribution_Symptoms_to_ICU, Distribution_Hopital_to_Death,...
    Hosp_2_Death, Sympt_2_critcal, Sympt_2_hosp,  WALES_FLAG);
[~, ~, ~, nICU, inICU, nDeaths, ~]=ODE_to_Obs(T, nDCD(:,aa), E(:,aa), Region, h_factor(1),i_factor(1),d_factor(1),h_stretch(1),i_stretch(1),a,...
    Assumed_Delay_Reporting_Deaths, Distribution_Hosp_Time', Distribution_HospICU_Time', Distribution_ICU_Time', Distribution_Symptoms_to_Hospital, Distribution_Symptoms_to_ICU, Distribution_Hopital_to_Death,...
    Hosp_2_Death, Sympt_2_critcal, Sympt_2_hosp,  WALES_FLAG);

[NDC(:,aa+L), nHospital(:,aa+L), inHospital(:,aa+L), nICU(:,aa+L), inICU(:,aa+L), nDeaths(:,aa+L), INF(:,aa+L)]=ODE_to_Obs(T, nDCH(:,aa+L), E(:,aa+L), Region, h_factor(m2),i_factor(m2),d_factor(m2),h_stretch(m2),i_stretch(m2),a,...
    Assumed_Delay_Reporting_Deaths, Distribution_Hosp_Time', Distribution_HospICU_Time', Distribution_ICU_Time', Distribution_Symptoms_to_Hospital, Distribution_Symptoms_to_ICU, Distribution_Hopital_to_Death,...
    Hosp_2_Death, Sympt_2_critcal, Sympt_2_hosp,  WALES_FLAG);
[~, ~, ~, nICU(:,aa+L), inICU(:,aa+L), nDeaths(:,aa+L), ~]=ODE_to_Obs(T, nDCD(:,aa+L), E(:,aa+L), Region, h_factor(m2),i_factor(m2),d_factor(m2),h_stretch(m2),i_stretch(m2),a,...
    Assumed_Delay_Reporting_Deaths, Distribution_Hosp_Time', Distribution_HospICU_Time', Distribution_ICU_Time', Distribution_Symptoms_to_Hospital, Distribution_Symptoms_to_ICU, Distribution_Hopital_to_Death,...
    Hosp_2_Death, Sympt_2_critcal, Sympt_2_hosp,  WALES_FLAG);

[~, nHospital(:,aa+LL), inHospital(:,aa+LL), nICU(:,aa+LL), inICU(:,aa+LL), nDeaths(:,aa+LL), INF(:,aa+LL)]=ODE_to_Obs(T, nDCH(:,aa+LL), E(:,aa+LL), Region, h_factor(m3),i_factor(m3),d_factor(m3),h_stretch(m3),i_stretch(m3),a,...
    Assumed_Delay_Reporting_Deaths, Distribution_Hosp_Time', Distribution_HospICU_Time', Distribution_ICU_Time', Distribution_Symptoms_to_Hospital, Distribution_Symptoms_to_ICU, Distribution_Hopital_to_Death,...
    Hosp_2_Death, Sympt_2_critcal, Sympt_2_hosp,  WALES_FLAG);
[~, ~, ~, nICU(:,aa+LL), inICU(:,aa+LL), nDeaths(:,aa+LL), ~]=ODE_to_Obs(T, nDCD(:,aa+LL), E(:,aa+LL), Region, h_factor(m3),i_factor(m3),d_factor(m3),h_stretch(m3),i_stretch(m3),a,...
    Assumed_Delay_Reporting_Deaths, Distribution_Hosp_Time', Distribution_HospICU_Time', Distribution_ICU_Time', Distribution_Symptoms_to_Hospital, Distribution_Symptoms_to_ICU, Distribution_Hopital_to_Death,...
    Hosp_2_Death, Sympt_2_critcal, Sympt_2_hosp,  WALES_FLAG);

end

%%


