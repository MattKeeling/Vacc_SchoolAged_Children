clear
Generate_Output('06_11_21');

%%
Color=[0.5 0.5 0.5; 1 0 0; 0 0.9 0];

QQ=[0 3 5];

figure(1); clf;

for LOOP=0:2
     
    load(['Warwick_Output_Loop' num2str(LOOP) '_06_11_21.mat'])
    T=1:size(nALL,2); Weighting=1+0*T; Weighting(1:(datenum(2021,7,18)+1-datenum(2020,1,1)))=0;
    Weighting((datenum(2022,1,1)+1-datenum(2020,1,1)):(datenum(2022,12,31)+1-datenum(2020,1,1)))=1-[1:365]/365;
    Weighting((datenum(2023,1,1)+1-datenum(2020,1,1)):end)=0;
    

    Y=nALL; 
    subplot(4,3,1); y=Y(:,:,2)+0.4*Y(:,:,3); 
    fill(LOOP+0.3*[-1 1 1 -1],[0 0 1 1]*(sum(y(2:8,:),1)*Weighting'),'r','FaceColor',Color(LOOP+1,:)); hold on
    set(gca,'XTick',[0:2],'XTickLabel',{'>18 only','12-17','5-17'});
    ylabel('Aged 5-11'); title('Infections'); set(gca,'XLim',[-0.5 2.5]);
    
    subplot(4,3,4); y=0.6*Y(:,:,3)+0.6*Y(:,:,4); 
    fill(LOOP+0.3*[-1 1 1 -1],[0 0 1 1]*(sum(y(2:8,:),1)*Weighting'),'r','FaceColor',Color(LOOP+1,:)); hold on
    set(gca,'XTick',[0:2],'XTickLabel',{'>18 only','12-17','5-17'});
    ylabel('Aged 12-17'); set(gca,'XLim',[-0.5 2.5]);
    
    subplot(4,3,7); y=Y(:,:,2)+Y(:,:,3)+0.6*Y(:,:,4); 
    fill(LOOP+0.3*[-1 1 1 -1],[0 0 1 1]*(sum(y(2:8,:),1)*Weighting'),'r','FaceColor',Color(LOOP+1,:)); hold on
    set(gca,'XTick',[0:2],'XTickLabel',{'>18 only','12-17','5-17'});
    ylabel('Aged 5-17'); title('Infections'); set(gca,'XLim',[-0.5 2.5]);
    
    subplot(4,3,10); y=sum(Y(:,:,:),3); 
    fill(LOOP+0.3*[-1 1 1 -1],[0 0 1 1]*(sum(y(2:8,:),1)*Weighting'),'r','FaceColor',Color(LOOP+1,:)); hold on
    set(gca,'XTick',[0:2],'XTickLabel',{'>18 only','12-17','5-17'});
    ylabel('All Ages'); xlabel('Vaccinated ages');
    set(gca,'XLim',[-0.5 2.5]);

    Y=nHOSP_AD; 
    subplot(4,3,2); y=Y(:,:,2)+0.4*Y(:,:,3); 
    fill(LOOP+0.3*[-1 1 1 -1],[0 0 1 1]*(sum(y(2:8,:),1)*Weighting'),'r','FaceColor',Color(LOOP+1,:)); hold on
    set(gca,'XTick',[0:2],'XTickLabel',{'>18 only','12-17','5-17'});
    ylabel('Aged 5-11'); title('Hosptial Admissions'); set(gca,'XLim',[-0.5 2.5]);
    
    subplot(4,3,5); y=0.6*Y(:,:,3)+0.6*Y(:,:,4); 
    fill(LOOP+0.3*[-1 1 1 -1],[0 0 1 1]*(sum(y(2:8,:),1)*Weighting'),'r','FaceColor',Color(LOOP+1,:)); hold on
    set(gca,'XTick',[0:2],'XTickLabel',{'>18 only','12-17','5-17'});
    ylabel('Aged 12-17'); set(gca,'XLim',[-0.5 2.5]);
    
    subplot(4,3,8); y=Y(:,:,2)+Y(:,:,3)+0.6*Y(:,:,4); 
    fill(LOOP+0.3*[-1 1 1 -1],[0 0 1 1]*(sum(y(2:8,:),1)*Weighting'),'r','FaceColor',Color(LOOP+1,:)); hold on
    set(gca,'XTick',[0:2],'XTickLabel',{'>18 only','12-17','5-17'});
    ylabel('Aged 5-17'); title('Infections'); set(gca,'XLim',[-0.5 2.5]);
    
    subplot(4,3,11); y=sum(Y(:,:,:),3); 
    fill(LOOP+0.3*[-1 1 1 -1],[0 0 1 1]*(sum(y(2:8,:),1)*Weighting'),'r','FaceColor',Color(LOOP+1,:)); hold on
    set(gca,'XTick',[0:2],'XTickLabel',{'>18 only','12-17','5-17'});
    ylabel('All Ages'); xlabel('Vaccinated ages');
    set(gca,'XLim',[-0.5 2.5]);

    Y=nDEATHS;
    subplot(4,3,3); y=Y(:,:,2)+0.4*Y(:,:,3); 
    fill(LOOP+0.3*[-1 1 1 -1],[0 0 1 1]*(sum(y(2:8,:),1)*Weighting'),'r','FaceColor',Color(LOOP+1,:)); hold on
    set(gca,'XTick',[0:2],'XTickLabel',{'>18 only','12-17','5-17'});
    ylabel('Aged 5-11'); title('Deaths'); set(gca,'XLim',[-0.5 2.5]);
    
    subplot(4,3,6); y=0.6*Y(:,:,3)+0.6*Y(:,:,4); 
    fill(LOOP+0.3*[-1 1 1 -1],[0 0 1 1]*(sum(y(2:8,:),1)*Weighting'),'r','FaceColor',Color(LOOP+1,:)); hold on
    set(gca,'XTick',[0:2],'XTickLabel',{'>18 only','12-17','5-17'});
    ylabel('Aged 12-17'); set(gca,'XLim',[-0.5 2.5]);
    
    subplot(4,3,9); y=Y(:,:,2)+Y(:,:,3)+0.6*Y(:,:,4); 
    fill(LOOP+0.3*[-1 1 1 -1],[0 0 1 1]*(sum(y(2:8,:),1)*Weighting'),'r','FaceColor',Color(LOOP+1,:)); hold on
    set(gca,'XTick',[0:2],'XTickLabel',{'>18 only','12-17','5-17'});
    ylabel('Aged 5-17'); title('Infections'); set(gca,'XLim',[-0.5 2.5]);
    
    subplot(4,3,12); y=sum(Y(:,:,:),3); 
    fill(LOOP+0.3*[-1 1 1 -1],[0 0 1 1]*(sum(y(2:8,:),1)*Weighting'),'r','FaceColor',Color(LOOP+1,:)); hold on
    set(gca,'XTick',[0:2],'XTickLabel',{'>18 only','12-17','5-17'});
    ylabel('All Ages'); xlabel('Vaccinated ages');
    set(gca,'XLim',[-0.5 2.5]);
    drawnow;
end


