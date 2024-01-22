%% LOAD BEHAVIOR FILE
loadbehavior =1
if loadbehavior ==1
    server = 'D:\Dropbox'
    cd ([  server '\Work_related\Duration opto_2' ])
    clear all
    load DataB DataB
end

thismarkerdize =30
ErrobarlineWidth= 1
%%
% Select Animal actions (left 0 vs right 1 turn)
Animal_Action_Left_Right = DataB.Response.Action;
% [6 8 2 3 5 9 16 13]these rats have swithced ruls that is T2>T1 go left
% and the rest T2>T1 go right. 
Animal_Action_Left_Right(ismember(DataB.Info.Rat,[6 8 2 3 5 9 16 13]))= 1-Animal_Action_Left_Right(ismember(DataB.Info.Rat,[6 8 2 3 5 9 16 13]));
thisactionnew=Animal_Action_Left_Right;
NTD= (DataB.Stimulus.Duration{2} -DataB.Stimulus.Duration{1} )./(DataB.Stimulus.Duration{2} +DataB.Stimulus.Duration{1} );
NTD = (round(NTD*100)/100);
opto = DataB.Stimulus.opto;
Pscyometric_dataromance = DataB.Response.Correctness;




Bootstrap_percnet_choice =[];
Fit_Parameteres=[];
theseNTD = [-0.3485   -0.2400      -0.1170       0  0.1170     0.2400    0.3496 ];
Pscyometric_data=[];
Correct=[];
thyesedates = unique(DataB.Info.Date(ismember(DataB.Info.Rat ,9)));





%%


%% creat figure  sup 2d halo
Rats = [9 10];
thiscolor=[60 88 167;0 0 0; 93 186 70]./255;
bin=theseNTD;
figure
thesedays = unique(DataB.Info.Date );
thefirst_7= find(opto==7,1,'first');
thuissum=[]
perf=[];
% andthiscg =(DataB.Stimulus.Duration{1} <400&DataB.Stimulus.Duration{1} >200&ismember(DataB.Info.Rat ,Rat))&ismember(DataB.Info.Date ,[thesedays])&DataB.Stimulus.Nominal{2} ==DataB.Stimulus.Nominal{1} ;
for Ratid=1:numel(Rats)
    thesedays=unique(DataB.Info.Date(ismember(DataB.Info.Rat,Rats(Ratid))));
        if Ratid==5
        thesedays=thesedays(~ismember(thesedays,[  737103    737104  737107]));
        end
    andthiscg =DataB.Info.TrialNumber>0&DataB.Stimulus.Duration{1}<335&...
        DataB.Stimulus.Duration{1}>333&ismember(DataB.Info.Rat,Rats(Ratid))...
        &DataB.Stimulus.Nominal{2}>78&...
        DataB.Stimulus.Nominal{1}>78&ismember(DataB.Info.Date ,thesedays([ 1:end]))&DataB.Response.Perfromance>0.5;

% NTDcounter=0;

    for thiss=1:size(theseNTD,2)

        % NTDcounter=NTDcounter+1;
        thiscondition = NTD>theseNTD(thiss)-0.05&NTD<theseNTD(thiss)+0.05&andthiscg;
        thuissum{1}(Ratid,thiss) = sum((opto==7)&thiscondition);
        thuissum{3}(Ratid,thiss) = sum((opto==5)&thiscondition);
        thuissum{2}(Ratid,thiss) = sum((opto==3|opto==1)&thiscondition);

        thuissumaction{1}(Ratid,thiss) = sum(thisactionnew((opto==7)&thiscondition));
        thuissumaction{3}(Ratid,thiss) = sum(thisactionnew((opto==5)&thiscondition));
        thuissumaction{2}(Ratid,thiss) = sum(thisactionnew((opto==3|opto==1)&thiscondition));


        perf{1}(Ratid,thiss)=mean(thisactionnew((opto==7)&thiscondition));
        %     perf(2,thiss)=mean(thisactionnew(DataB.Stimulus.Nominal{2} >DataB.Stimulus.Nominal{1} &(opto==3|opto==1)&thiscondition));
        perf{2}(Ratid,thiss)=mean(thisactionnew(DataB.Stimulus.Nominal{1}<90 ...
            &DataB.Stimulus.Nominal{2}<90&DataB.Stimulus.Nominal{2} ==DataB.Stimulus.Nominal{1}...
            &(opto==3|opto==1)&thiscondition));
        %      perf(2,thiss)=mean(thisactionnew(DataB.Stimulus.Nominal{1} <90&DataB.Stimulus.Nominal{2} <90&(opto==3|opto==1)&thiscondition));
        perf{3}(Ratid,thiss)=mean(thisactionnew((opto==5)&thiscondition));
        perf{4}(Ratid,thiss)=mean(thisactionnew(DataB.Stimulus.Nominal{2} <DataB.Stimulus.Nominal{1} &(opto==3|opto==1)&thiscondition));
        perf{5}(Ratid,thiss)=mean(thisactionnew(DataB.Stimulus.Nominal{1} <90&DataB.Stimulus.Nominal{2} ==DataB.Stimulus.Nominal{1} &(opto==3|opto==1)&thiscondition));
        perf{6}(Ratid,thiss)=mean(thisactionnew(DataB.Stimulus.Nominal{1} >90&DataB.Stimulus.Nominal{2} ==DataB.Stimulus.Nominal{1} &(opto==3|opto==1)&thiscondition));
        perf{6}(Ratid,thiss)=mean(thisactionnew(DataB.Stimulus.Nominal{1} >90&DataB.Stimulus.Nominal{2} ==DataB.Stimulus.Nominal{1} &(opto==3|opto==1)&thiscondition));


    end

end




F=@(g,I,u,v,x) (g+(1-g-I)./(1+exp(-(x-u)/v)))*100;

UL=[0.5, 0.5, 10, 5];
SP=[ 0.4, 0.4, 5, 0.2];
LM=[0, 0, -10, 0.1];


UL=[0.5, 0.5, 10, 6];
SP=[ 0.4, 0.4, 5, 0.3];
LM=[0, 0, -10, 0.01];
plot([0 0], [ 0 100],'-','color',[0.5 0.5 0.5],'Linewidth',2)
hold on
plot([-0.36 0.36], [ 50 50],'-','color',[0.5 0.5 0.5],'Linewidth',2)


for i=[1  3  ]

    bootfit= (fit(bin',mean(perf{i})'*100,F,'StartPoint',SP,'Upper',UL,'Lower',LM));
    hold on
    hL1 = plot(bootfit,'-');
    hold on
    set([hL1 ],'color',thiscolor(i,:),'Linewidth',2)
    plot(bin, mean(perf{i})*100,'.','color',thiscolor(i,:),'Markersize',thismarkerdize)
    h=errorbar(bin, mean(perf{i})*100,std(perf{i}*100)./sqrt(size(perf{i},1)),...
        'color',thiscolor(i,:),'Markersize',30);
    h.LineWidth = ErrobarlineWidth;
    h.LineStyle = 'none';

end


title('External light duration rats (Control)')
set(gca,'xtick',[-0.3:0.1:0.3])
set(gca,'ytick',[0:20:100])
set(gca,'xtick',[-0.2 0 0.2])
set(gca,'xtickLabel',[-0.2 0 0.2])
ylabel('Choice T2>T1 (%)')
xlabel('\DeltaT')


set(gca,'FontSize',15)
hold on
axis square
set(gca,'TickDir', 'out' )
box off

xlim([-0.36 0.36])
ylim([0 100])
legend off
figurecorrectionOPTO(gca,[0.05 0.05],12)

%% suplimentary fig 3


%%
hold off
thuissumaction=[];
thuissum=[];
%   DataB.Stimulus.Nominal{2}
thisactionnew = DataB.Response.Action;
thisactionnew(ismember(DataB.Info.Rat,[6 8 2 3 5 9]))=1-thisactionnew(ismember(DataB.Info.Rat,[6 8 2 3 5 9]));
correct = DataB.Response.Correctness;
NTD = (DataB.Stimulus.Duration{2} -DataB.Stimulus.Duration{1} )./(DataB.Stimulus.Duration{2} +DataB.Stimulus.Duration{1} );



[thesepairs, ~] = unique([DataB.Stimulus.Duration{2} ;DataB.Stimulus.Duration{1} ]','rows');


perf = [];
number = [];
theseNTD = [-0.3485   -0.2400      -0.1170       0  0.1170     0.2400    0.3496 ];
NTDcounter = 0;
SDI = (DataB.Stimulus.Nominal{2} -DataB.Stimulus.Nominal{1} )./(DataB.Stimulus.Nominal{2} +DataB.Stimulus.Nominal{1} );


% close all
NTD = (DataB.Stimulus.Duration{2} -DataB.Stimulus.Duration{1} )./(DataB.Stimulus.Duration{2} +DataB.Stimulus.Duration{1} );

opto = DataB.Stimulus.opto;
%  Rat =[1:4 ];
Rats =unique(DataB.Info.Rat(DataB.Info.Optoinfo==4));
Rats = [   9 15]
Rats = [ 10    11        14 15]
Rats =unique(DataB.Info.Rat(DataB.Info.Optoinfo==1));
Rats = 1:4
Rats = [1:5];
% Rats =unique(DataB.Info.Rat(DataB.Info.Optoinfo==1));
%   Rats =[9 10]
% DataB.Info.Rat(DataB.Info.Rat==16)=15;

thesedays = unique(DataB.Info.Date );
thefirst_7= find(opto==7,1,'first');
thuissum = []
perf=[];
% andthiscg =(DataB.Stimulus.Duration{1} <400&DataB.Stimulus.Duration{1} >200&ismember(DataB.Info.Rat ,Rat))&ismember(DataB.Info.Date ,[thesedays])&DataB.Stimulus.Nominal{2} ==DataB.Stimulus.Nominal{1} ;
%  for Ratid=1:numel(Rats)
   for Ratid=1:5

    andthiscg =DataB.Info.TrialNumber>0&DataB.Stimulus.Duration{1}<335&...
        DataB.Stimulus.Duration{1}>333&ismember(DataB.Info.Rat ,Rats(Ratid))...
        &ismember(DataB.Info.Date ,thesedays([1:end]))&DataB.Response.Perfromance>0;
    for thiss=1:size(theseNTD,2)
        
        NTDcounter=NTDcounter+1;
%         thiscondition = (NTD>theseNTD(thiss)-0.05&NTD<theseNTD(thiss)+0.05&andthiscg)&(DataB.Stimulus.Nominal{1}<90&DataB.Stimulus.Nominal{2}<90)&(DataB.Stimulus.Nominal{1}>70&DataB.Stimulus.Nominal{2}>70);
        thiscondition = (NTD>theseNTD(thiss)-0.05&NTD<theseNTD(thiss)+0.05&andthiscg&DataB.Stimulus.Nominal{1}==80&DataB.Stimulus.Nominal{2}==80);
       
        thuissum{1}(Ratid,thiss) = sum((opto==7)&thiscondition);
        thuissum{3}(Ratid,thiss) = sum((opto==5)&thiscondition);
        thuissum{2}(Ratid,thiss) = sum((opto==3|opto==1)&thiscondition);
        
        thuissumaction{1}(Ratid,thiss) = sum(thisactionnew((opto==7)&thiscondition));
        thuissumaction{3}(Ratid,thiss) = sum(thisactionnew((opto==5)&thiscondition));
        thuissumaction{2}(Ratid,thiss) = sum(thisactionnew((opto==3|opto==1)&thiscondition));
        
        
        perf{1}(Ratid,thiss)=mean(thisactionnew((opto==7)&thiscondition));
        %     perf(2,thiss)=mean(thisactionnew(DataB.Stimulus.Nominal{2} >DataB.Stimulus.Nominal{1} &(opto==3|opto==1)&thiscondition));
        perf{2}(Ratid,thiss)=mean(thisactionnew(DataB.Stimulus.Nominal{1}<90 ...
            &DataB.Stimulus.Nominal{2}<90&DataB.Stimulus.Nominal{2} ==DataB.Stimulus.Nominal{1}...
            &(opto==3|opto==1)&thiscondition));
        %      perf(2,thiss)=mean(thisactionnew(DataB.Stimulus.Nominal{1} <90&DataB.Stimulus.Nominal{2} <90&(opto==3|opto==1)&thiscondition));
        perf{3}(Ratid,thiss)=mean(thisactionnew((opto==5)&thiscondition));
        perf{4}(Ratid,thiss)=mean(thisactionnew(DataB.Stimulus.Nominal{2} <DataB.Stimulus.Nominal{1} &(opto==3|opto==1)&thiscondition));
        perf{5}(Ratid,thiss)=mean(thisactionnew(DataB.Stimulus.Nominal{1} <90&DataB.Stimulus.Nominal{2} ==DataB.Stimulus.Nominal{1} &(opto==3|opto==1)&thiscondition));
        perf{6}(Ratid,thiss)=mean(thisactionnew(DataB.Stimulus.Nominal{1} >90&DataB.Stimulus.Nominal{2} ==DataB.Stimulus.Nominal{1} &(opto==3|opto==1)&thiscondition));
        perf{6}(Ratid,thiss)=mean(thisactionnew(DataB.Stimulus.Nominal{1} >90&DataB.Stimulus.Nominal{2} ==DataB.Stimulus.Nominal{1} &(opto==3|opto==1)&thiscondition));
        
        
    end
    
end






bin=theseNTD;
thiscolor=[0 0 1; 0 1 0; 1 0 0];
thiscolor=[0 0 0; 0 1 0.5; 0.8 0 0.6];
thiscolor=[0 0 0; 0.4 0.4 0.4; 0.7 0.7 0.7];
thiscolor=[0 0.4 1;0 0 0; 0 1 0; 0.5 0.5 0.5; 0.5 0.5 1; 0 0.5 0.5];

F=@(g,I,u,v,x) (g+(1-g-I)./(1+exp(-(x-u)/v)))*100;

UL=[0.5, 0.5, 10, 5];
SP=[ 0.4, 0.4, 5, 0.2];
LM=[0, 0, -10, 0.1];


UL=[0.5, 0.5, 10, 5];
SP=[ 0.4, 0.4, 5, 0.3];
LM=[0, 0, -10, 0.02];

for Ratid = 1:numel(Rats)
        hold off

    subplot(3,2,6-Ratid)
    plot([0 0], [ 0 100],'-','color',[0.5 0.5 0.5],'Linewidth',2)
    hold on
    plot([-0.36 0.36], [ 50 50],'-','color',[0.5 0.5 0.5],'Linewidth',2)
    psycho=thesepairs(:,2)==334;
    for i=[ 1 3 ]
        
        bootfit= (fit(bin',(perf{i}(Ratid,:))'*100,F,'StartPoint',SP,'Upper',UL,'Lower',LM));
        hold on
        hL1 = plot(bootfit,'-');
        hold on
        set([hL1 ],'color',thiscolor(i,:),'Linewidth',4)
        plot(bin, (perf{i}(Ratid,:))*100,'.','color',thiscolor(i,:),'Markersize',40)
        
        
    end
    
    
    % set(gca,'ytickLabel',[])
    set(gca,'xtick',[-0.3:0.1:0.3])
    set(gca,'ytick',[0:20:100])
    set(gca,'xtick',[-0.2 0 0.2])
    set(gca,'xtickLabel',[-0.2 0 0.2])
 ylabel('Choice T2>T1 (%)')
xlabel('\DeltaT')
    
    
    set(gca,'FontSize',15)
    hold on
    axis square
    set(gca,'TickDir', 'out' )
    box off
    
    xlim([-0.36 0.36])
    ylim([0 100])
    legend off
    figurecorrectionOPTO(gca,[0.05, 0.05],12)
end