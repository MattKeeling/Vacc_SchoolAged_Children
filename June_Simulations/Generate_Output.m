%%
clear all

Remove_New_Variant=0;

Matlab_Save=1;

DATE_STR='26_06_21';

load(['Data_File_' DATE_STR '.mat']);

KL=0;

ZZ=1;

Names={'England','East of England','London','Midlands','North East and Yorkshire','North West','South East','South West','Wales','Scotland','Northern Ireland','UK'};

for VACC_LOOP=0:2 % What type of childhood vaccine is deployed
    
    L1=[10 10 10 10 2 4 7]; L2=[1 4 7 10 2 4 7];
    
    for PreCBeh=length(L1):-1:1
        
        Loop=L1(PreCBeh)-1;
        
        Loop2=L2(PreCBeh)-1;
        
        
        KL=KL+1;
        
        PART=4;
        
        EndValue=0.1; EV2=0.1;   RollOutSlower=1.0;
        
        STEP2=0.7; STEP3=0.4;
        STEP4_TIME=datenum(2021,6,21)+1-datenum(2020,1,1);
        Seasonality=0.1;  Vspeed=1;  XTrans=1.0; XTransDate=300;
        RollOutSlower=1.0;
        
        STEP4_TIME=max(datenum(2021,7,19)+1-datenum(2020,1,1),datenum(DATE_STR,'dd_mm_yy')+1-datenum(2020,1,1));
        STEP4_LEVEL=Loop/10;
        STEP4_UNTIL=STEP4_TIME+3*7+3*7*Loop2;
        
        switch VACC_LOOP
            case 0
                UT=0.95*ones(1,21); UT(15:21)=0.95; UT(5:10)=0.8; UT(4)=0.4*UT(5); UT(1:3)=0.0; UT(ALREADY_VACC)=0;  % No Vacc in under 18's
            case 1
                UT=0.95*ones(1,21); UT(15:21)=0.95; UT(5:10)=0.8; UT(4)=UT(5); UT(3)=0.6*UT(5); UT(1:2)=0.0; UT(ALREADY_VACC)=0;  % Vacc 12-17
            case 2
                UT=0.95*ones(1,21); UT(15:21)=0.95; UT(5:10)=0.8; UT(2:4)=UT(5); UT(1)=0.0; UT(ALREADY_VACC)=0; % Vacc 5-17
        end
        
        
        [V1,V2,VET1,VET2,RC,RCH, RCD,Phase, RatioPf]=Vaccination_over_Time(UT,VIPfAZ,VSPfAZ,VHPfAZ,VDPfAZ,Vspeed,RollOutSlower,datenum(DATE_STR,'dd_mm_yy')+1-datenum(2020,1,1));
        
        
        %%
        clear nALL nDEATHS SYMPTOMS HOSPITALISED nHOSP_AD nCASES PILLAR2
        
        load Other_Data.mat
        
        Names={'England','East of England','London','Midlands','North East and Yorkshire','North West','South East','South West','Wales','Scotland','Northern Ireland','United Kingdom'};
        
        Comps = zeros(4,Num_Comp+1,11);
        
        tmpRS=RUN_STOP;  Z=(RUN_STOP(end-1)+7):7:(ceil(datenum(2023,1,14)+1-datenum(2020,1,1)));
        tmpRS(Num_Comp+[1:length(Z)]-1)=Z;
        
        RUN_STOPs = zeros(11,Num_Comp+1);
        for Region = 2:11
            RUN_STOPs(Region,1:length(tmpRS)) = tmpRS;
        end
        
        Keep_RUN_STOPs=RUN_STOPs(2,:);
        
        
        %%
        maxtime = max(RUN_STOPs(end,:))+30; % needs to be big enough to take account of any lags
        nALL = zeros(11,maxtime,63);
        nDEATHS = zeros(11,maxtime,63);
        nHOSP_AD = zeros(11,maxtime,63);
        nHOSP_OCC = zeros(11,maxtime,63);
        nICU_AD = zeros(11,maxtime,63);
        nICU_OCC = zeros(11,maxtime,63);
        
        clear Comps
        
        RUN_STARTs=[82*ones(11,1) RUN_STOPs(:,1:(end-1))];
        
        Christmas=0;
        
        % SET UP THE COMPILANCE
        for Region = 2:11
            
            RUN_START=[82 RUN_STOPs(Region,1:(end-1))];
            
            m=find(RUN_START>=datenum(2021,2,1)+1-datenum(2020,1,1) & RUN_STOPs(Region,:)<=datenum(2021,3,12)+1-datenum(2020,1,1));
            LDComp=max(COMPLIANCE(m,Region));
            
            Comps(Region,1:Num_Comp)=COMPLIANCE(1:Num_Comp,Region);
            Comps(Region,Num_Comp:size(RUN_STOPs,2))=Comps(Region,Num_Comp);
            
            m=find(RUN_START>=STEP4_TIME); % Start STEP 4
            SeasonDate=RUN_START(RUN_START>datenum(DATE_STR,'dd_mm_yy')+1-datenum(2020,1,1));
            
            for i=1:length(m)
                V=STEP4_LEVEL*(1-(RUN_START(m(i))-RUN_START(m(1)))/(STEP4_UNTIL-STEP4_TIME));
                V(V<0)=0;
                Comps(Region,m(i))=V*Comps(Region,Num_Comp);
                Comps(Region,m(i))=Comps(Region,m(i)).*SeasonScaling(Region,RUN_START(m(1)));
            end
            
            
        end
        
        MAX_REGION=8;
        
        % Main Run !
        
        Christmas=0;
        
        if Remove_New_Variant
            IMPORT_DAY(:)=max(max(RUN_STOPs))+100;
        end
        
        %parfor Region=2:MAX_REGION
        for Region=2:MAX_REGION
            
            fprintf(1,'Simulating %s, PreCautionaryBehaviour %d, Vaccine Strategy %d\n',Names{Region}, PreCBeh, VACC_LOOP);
            
            [nDC, nDCD, nDCH, nHospital, inHospital, nICU, inICU, nDeaths, ~, All]=Simulate_One_Region(Region, TAU, ALPHA, INC_P, SCALING(Region), FACTOR(Region), H_FACTOR(Region+[0 11 22]), I_FACTOR(Region+[0 11 22]), D_FACTOR(Region+[0 11 22]),  H_STRETCH(Region+[0 11 22]), ...
                I_STRETCH(Region+[0 11 22]), LAG(Region), START_DATE(Region)+1, 1, Comps(Region,:), RUN_STOPs(Region,:), [NV_BETA(Region+[0]) NV_BETA(Region+[11])], [NV_SPEED(Region+[0 11])]', [IMPORT_DAY(Region+[0 11])], IMPORT_LEVEL(Region+[0 11]),...
                squeeze(V1(Region,:,:)), squeeze(V2(Region,:,:)), VET1(Region,:), VET2(Region,:), Transmission_Reduction, squeeze(RC(Region,:,:)), squeeze(RCH(Region,:,:)), squeeze(RCD(Region,:,:)),VoC_Sus, General_Reduction,...
                RatioPf, Detection, Susceptibility, GAMMA, [Seasonality SeasonDate], XTrans, XTransDate);
            
            RS=sum(inHospital,2)./sum(Region_PP(Region,:),2);  RS(RS>0.002)=0.002;
            
            nDeaths=nDeaths.* ((1+DH_SCALING(Region)*RS)*ones(1,size(nDeaths,2)));
            
            T=1:size(All,1);
            
            padding = maxtime-length(T);
            nALL(Region,:,:)=[All; zeros(padding,63)];
            nDEATHS(Region,:,:)=[nDeaths; zeros(padding,63)];
            nHOSP_OCC(Region,:,:)=[inHospital; zeros(padding,63)];
            nHOSP_AD(Region,:,:) = [nHospital; zeros(padding,63)];
            nICU_OCC(Region,:,:) = [inICU; zeros(padding,63)];
            nICU_AD(Region,:,:) = [nICU; zeros(padding,63)];
        end
        
        
        %% SAVE OUTPUT
        if Matlab_Save
            
            ALL(PreCBeh,:,:,:)=nALL(:,:,1:21)+nALL(:,:,22:42)+nALL(:,:,43:63);
            DEATHS(PreCBeh,:,:,:)=nDEATHS(:,:,1:21)+nDEATHS(:,:,22:42)+nDEATHS(:,:,43:63);
            HOSP_OCC(PreCBeh,:,:,:)=nHOSP_OCC(:,:,1:21)+nHOSP_OCC(:,:,22:42)+nHOSP_OCC(:,:,43:63);
            HOSP_AD(PreCBeh,:,:,:)=nHOSP_AD(:,:,1:21)+nHOSP_AD(:,:,22:42)+nHOSP_AD(:,:,43:63);
            ICU_OCC(PreCBeh,:,:,:)=nICU_OCC(:,:,1:21)+nICU_OCC(:,:,22:42)+nICU_OCC(:,:,43:63);
            ICU_AD(PreCBeh,:,:,:)=nICU_AD(:,:,1:21)+nICU_AD(:,:,22:42)+nICU_AD(:,:,43:63);
            
            save(['Warwick_Output_Loop' num2str(VACC_LOOP) '_' DATE_STR '.mat'], 'ALL', 'DEATHS', 'HOSP_OCC', 'HOSP_AD', 'ICU_OCC','ICU_AD');
            
        end
    end
end


