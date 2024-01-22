   % load('D:\Dropbox\Work_related\Duration opto\datastr\DataB.mat')
 % yourfolder = 'D:\Dropbox\OptoDurationPaper\A_NATCOMSUBMISSION';
   load(fullfile('D:\Dropbox\OptoDurationPaper\A_NATCOMSUBMISSION','Data_figure1.mat'))
  DataB = Data_figure1;
DataB.Stimulus.TDI= (DataB.Stimulus.Duration{2}-DataB.Stimulus.Duration{1})./(DataB.Stimulus.Duration{2}+DataB.Stimulus.Duration{1});
% Ratswithswitchrule T2>T1 go left or I2 >I1 go left) instead of right
changerats = [{[60005 60006 61001  61003]},{[61002 60020 60032 60031  60016]}];
allrats= [{ [60001 60005 60006 61001 60004  61003 ]}...
    ,{[60002   60008 60009 60002 60016 61002  60020  60021 60032 60031]}];
% 60007
% 
 % indices=~ismember(DataB.Info.Rat,[allrats{1} allrats{2}] );
  % DataB = removeDaysfromDataSTRDuration(DataB,indices);
Minperf=0.65;
thesetrials =DataB.Response.Perfromance>Minperf;
   
        SpSide= sign(DataB.Stimulus.SDI);
        DSide= sign(DataB.Stimulus.TDI);
        IntensityColor= [93.7 18 12.9]/100;
        DurationColor = [47.1 47.1 47.1 ]./100;
        close all
        Minperf = 0.65;
        Lenghts=[0.1 0.1];
        Fontsize = 12;
        SDI = [-0.3:0.1:0.3];
        TDI = linspace(-0.35,0.35,7);
%% figure 1 a and b
figure
for ii=1:2




    Axes{1} = subplot(2,3,2);
    hold on
    Axes{3} = subplot(2,3,1);
    hold on
    Axes{2} = subplot(2,3,5);
    hold on
    Axes{4} = subplot(2,3,4);
    hold on

    Axes{5} = subplot(2,3,3);
    hold on
    Axes{6} = subplot(2,3,6);
    hold on
    norm=1;


    %
    % Thessname{60001}='bi1NTD';
    % Thessname{60002}='bi2NSD';
    % Thessname{60007}='fi7NSD';
    % Thessname{60008}='sht8NSD';
    % Thessname{60009}='sht9NSD';
    % Thessname{60002}='bi2NSD';
    % Thessname{60016}='ve16NSD';
    % Thessname{61001}='ss1NTD';
    % Thessname{61002}='ss2NSD';
    % Thessname{60003}='ss3NSD';
    % Thessname{60004}='ss4NTD';
    % Thessname{60005}='ss5NTD';
    % Thessname{60006}='ss6NTD';
    % Thessname{60020}='sr20NSD';
    % Thessname{60021}='sr21NSD';
    % Thessname{60032}='ar32NSD';
    % Thessname{60031}='ar31NSD';







    theserats=[  allrats{ii} ];
    % DataB.Stimulus.TDI= (DataB.Stimulus.Duration{2}-DataB.Stimulus.Duration{1})./(DataB.Stimulus.Duration{2}+DataB.Stimulus.Duration{1});
    subplot (Axes{ii+2})
    counter=0;
    SpSide= sign(DataB.Stimulus.SDI);
    perf=[];
    DSide= sign(DataB.Stimulus.TDI);
    % DSide(DSide==0)=(randi(2,1,sum(DSide==0))*2)-3;
    minDelta =0.01;
    for i=1:numel(theserats)



        counter=counter+1;

        generaldoncition=abs(DataB.Stimulus.SDI)>0.25&abs(DataB.Stimulus.TDI)>0.25&thesetrials&ismember(DataB.Info.Rat,theserats(i));
        generaldoncition=abs(DataB.Stimulus.SDI)>0&abs(DataB.Stimulus.TDI)>0.01;
        generaldoncition = abs(DataB.Stimulus.SDI) <1000;
        %         alltrials(counter)=sum(generaldoncition);
        % if ismember(theserats(i),[60001 60004 ])
        if ~ismember(theserats(i),[changerats{ii}])



            perf{1}(i,:)= nanmean( DataB.Response.Action(generaldoncition&DataB.Stimulus.Duration{1}>300&DataB.Stimulus.Duration{1}<400&DataB.Stimulus.TDI>=minDelta&ismember(DataB.Info.Rat,theserats(i))));
            perf{2}(i,:)= nanmean( DataB.Response.Action(generaldoncition&DataB.Stimulus.Duration{1}>300&DataB.Stimulus.Duration{1}<400&DataB.Stimulus.TDI<=-minDelta&ismember(DataB.Info.Rat,theserats(i))));
            perf{3}(i,:)= nanmean( DataB.Response.Action(generaldoncition&DataB.Stimulus.Nominal{1}>70&DataB.Stimulus.Nominal{1}<90&DataB.Stimulus.SDI>=minDelta&ismember(DataB.Info.Rat,theserats(i))));
            perf{4}(i,:)= nanmean( DataB.Response.Action(generaldoncition&DataB.Stimulus.Nominal{1}>70&DataB.Stimulus.Nominal{1} <90&DataB.Stimulus.SDI<=-minDelta&ismember(DataB.Info.Rat,theserats(i))));


        else

            perf{1}(i,:)= nanmean(1- DataB.Response.Action(generaldoncition&DataB.Stimulus.Duration{1}>300&DataB.Stimulus.Duration{1}<400&DataB.Stimulus.TDI>=minDelta&ismember(DataB.Info.Rat,theserats(i))));
            perf{2}(i,:)= nanmean( 1-DataB.Response.Action(generaldoncition&DataB.Stimulus.Duration{1}>300&DataB.Stimulus.Duration{1}<400&DataB.Stimulus.TDI<=-minDelta&ismember(DataB.Info.Rat,theserats(i))));
            perf{3}(i,:)= nanmean( 1-DataB.Response.Action(generaldoncition&DataB.Stimulus.Nominal{1}>70&DataB.Stimulus.Nominal{1}<90&DataB.Stimulus.SDI>=minDelta&ismember(DataB.Info.Rat,theserats(i))));
            perf{4}(i,:)= nanmean( 1-DataB.Response.Action(generaldoncition&DataB.Stimulus.Nominal{1}>70&DataB.Stimulus.Nominal{1}<90&DataB.Stimulus.SDI<=-minDelta&ismember(DataB.Info.Rat,theserats(i))));



        end

    end



    thisperf = 100*[[perf{1}+1-perf{2}]/2 [perf{3}+1-perf{4}]/2];


    hbar1 = bar ([1 ], nanmean(thisperf(:,1)));
    hold on
    hbar2 = bar ([ 3], nanmean(thisperf(:,2)));
    hbar2.FaceColor='none';
    hbar1.FaceColor='none';
    X = [ 1*ones(1,numel(thisperf(:,1)))+randn(1,numel(thisperf(:,1)))/5 ;...
        3*ones(1,numel(thisperf(:,2)))+randn(1,numel(thisperf(:,2)))/5 ];

    Y = [ thisperf(:,1)' ;...
        thisperf(:,2)' ];
    if ii==2
        plot (X ,Y,'color',[1 0.1 0.1])
        hold on
        plot (X ,Y,'o','color',[1 0.1 0.1])
        hold on
    else
        plot (X ,Y,'color',[0.4 0.4 0.4])
        hold on
        plot (X ,Y,'o','color',[0.4 0.4 0.4])
        hold on
    end
    h=errorbar([1 3], nanmean(thisperf),nanstd(thisperf)./sqrt(numel(theserats)*2));
    h.LineStyle='none';
    h.LineWidth=3;
    h.Color=[0 0 0];

    figurecorrectionOPTO(gca,[0.04 0.04],12)
    ylim([40 90])
    xlim([-0.5 4.5])
    set(gca, 'Xticklabel',['duration\newlinerule(T) ' ;'intensity\newlinerule(I)'])

    set(gca,'xtick',[1 3 ])
    ylabel('performance (%)')
end




perf1=[];
perf2 =[];
for iii=1:2
    theserats=[  allrats{iii} ];
    
    
    DataB.Response.Action(~ismember(DataB.Info.Rat,changerats{iii})) =...
        1-DataB.Response.Action(~ismember(DataB.Info.Rat,changerats{iii}));
    if iii==2
        DataB.Response.Action=1-DataB.Response.Action;
        
    end
    for ii=1:numel(theserats)
        
        
        DataB.Stimulus.TDI = (DataB.Stimulus.Duration{2}-DataB.Stimulus.Duration{1})./...
            (DataB.Stimulus.Duration{2}+DataB.Stimulus.Duration{1});
        
        
        SpSide= sign(DataB.Stimulus.SDI);
        DSide= sign(DataB.Stimulus.TDI);
        
        SDI = [-0.3:0.1:0.3];
        TDI = linspace(-0.35,0.35,7);
        
        for i =1:7
            
            perf1{iii}(ii,i)= nanmean( 1-DataB.Response.Action(DataB.Stimulus.Duration{1}>300&DataB.Stimulus.Duration{1}<500....
                &DataB.Stimulus.TDI>TDI(i)-0.05&DataB.Stimulus.TDI<TDI(i)+0.05&ismember(DataB.Info.Rat,theserats(ii))));
            
            perf2{iii}(ii,i)= nanmean(1- DataB.Response.Action(DataB.Stimulus.Nominal{1}>70&DataB.Stimulus.Nominal{1}<90&DataB.Stimulus.SDI>SDI(i)-0.05&...\
                DataB.Stimulus.SDI<SDI(i)+0.05&ismember(DataB.Info.Rat,theserats(ii))));
            
            
        end
    end
end

subplot(Axes{1})
F=@(g,I,u,v,x) (g+(1-g-I)./(1+exp(-(x-u)/v)))*100;
UL=[1, 1, 20, 10];
SP=[ 0.4, 0.4, 0, 0.2];
LM=[0, 0, -20, 0.05];
bin = TDI;
bootfit= (fit(bin',(nanmean(perf1{1})*100)',F,'StartPoint',SP,'Upper',UL,'Lower',LM));
hold on
hL1 = plot((Axes{1}) ,bin,bootfit(bin),'-');
hold on
set([hL1 ],'color',DurationColor,'Linewidth',4)


plot((Axes{1}) ,(bin), (nanmean(perf1{1})*100)','.','color',DurationColor,'Markersize',40)
bootfit= (fit(bin',(nanmean(perf1{2})*100)',F,'StartPoint',SP,'Upper',UL,'Lower',LM));
hold on
hL1 = plot((Axes{1}) ,bin,bootfit(bin),'-');
hold on
set([hL1 ],'color',IntensityColor,'Linewidth',4)
plot((Axes{1}) ,(bin),(nanmean(perf1{2})*100)','.','color',IntensityColor,'Markersize',40)
set(gca,'ytick',[0:20:100])
set(gca,'xtick',[-0.3 0 0.3])
% bin = round(155*(1+theseTDI)./(1-theseTDI));
set(gca,'xtickLabel',[-0.3 0 0.3])
ylabel('Choice T2>T1 (%)')
% xlabel('T2 ms (T1 = 155 ms)')
set(gca,'FontSize',15)
hold on
axis square
set(gca,'TickDir'     , 'out' )
xlim([-0.5 0.5])
ylim([0 100])
box off
legend off
figurecorrectionOPTO(gca,Lenghts,Fontsize)




subplot(Axes{2})


F=@(g,I,u,v,x) (g+(1-g-I)./(1+exp(-(x-u)/v)))*100;
UL=[0.5, 0.5, 10, 5];
SP=[ 0.4, 0.4, 5, 0.2];
LM=[0, 0, -10, 0.1];
bin = SDI;
bootfit= (fit(bin',(nanmean(perf2{1})*100)',F,'StartPoint',SP,'Upper',UL,'Lower',LM));
hold on
hL1 = plot((Axes{2}),bin,bootfit(bin),'-');
hold on
set([hL1 ],'color',DurationColor,'Linewidth',4)
plot(Axes{2},(bin), nanmean(perf2{1})*100,'.','color',DurationColor,'Markersize',40)
bootfit= (fit(bin',(nanmean(perf2{2})*100)',F,'StartPoint',SP,'Upper',UL,'Lower',LM));
hold on
hL1 = plot(Axes{2},bin,bootfit(bin),'-');
hold on
set([hL1 ],'color',IntensityColor,'Linewidth',4)
plot(Axes{2},(bin), nanmean(perf2{2})*100,'.','color',IntensityColor,'Markersize',40)
set(gca,'ytick',[0:20:100])
set(gca,'xtick',[-0.3 0 0.3])
% bin = round(155*(1+theseTDI)./(1-theseTDI));
set(gca,'xtickLabel',[-0.3 0 0.3])
ylabel('Choice I2>I1 (%)')
% xlabel('T2 ms (T1 = 155 ms)')
set(gca,'FontSize',15)
hold on
axis square
set(gca,'TickDir'     , 'out' )
xlim([-0.5 0.5])
ylim([0 100])
box off
legend off
figurecorrectionOPTO(gca,Lenghts,Fontsize)
%%  figure1 d

load(fullfile('D:\Dropbox\OptoDurationPaper\A_NATCOMSUBMISSION','Data_figure1.mat'))

DataB=Data_figure1;
perf1=[];
perf2 =[];
for iii=1:2
    theserats=[  allrats{iii} ];


    DataB.Response.Action(~ismember(DataB.Info.Rat,changerats{iii})) =...
        1-DataB.Response.Action(~ismember(DataB.Info.Rat,changerats{iii}));
    if iii==2
        DataB.Response.Action=1-DataB.Response.Action;

    end
    for ii=1:numel(theserats)


        DataB.Stimulus.TDI = (DataB.Stimulus.Duration{2}-DataB.Stimulus.Duration{1})./...
            (DataB.Stimulus.Duration{2}+DataB.Stimulus.Duration{1});


        SpSide= sign(DataB.Stimulus.SDI);
        DSide= sign(DataB.Stimulus.TDI);

        SDI = [-0.3:0.1:0.3];
        TDI = linspace(-0.35,0.35,7);

        for i =1:7

            perf1{iii}(ii,i)= nanmean( 1-DataB.Response.Action(DataB.Stimulus.Duration{1}>320&DataB.Stimulus.Duration{1}<340....
                &DataB.Stimulus.TDI>TDI(i)-0.05&DataB.Stimulus.TDI<TDI(i)+0.05&ismember(DataB.Info.Rat,theserats(ii))&abs(DataB.Stimulus.SDI)<0.2));

            perf2{iii}(ii,i)= nanmean(1- DataB.Response.Action(DataB.Stimulus.Nominal{1}>77&DataB.Stimulus.Nominal{1}<82&DataB.Stimulus.SDI>SDI(i)-0.02&...\
                DataB.Stimulus.SDI<SDI(i)+0.02&ismember(DataB.Info.Rat,theserats(ii))&abs(DataB.Stimulus.TDI)<0.2));



        end





    end
end

%%
% plot(nanmean(perf1{2}))
subplot(Axes{5})

for i=1:size(perf2{1},1)

    normperf(i,:) =(perf2{1}(i,:)-perf2{1}(i,4))*100;

end



% h1=plot(SDI,nanmean(normperf),'color',IntensityColor,'Markersize',40,'linewidth',4);
h1=plot(SDI,nanmean(normperf),'color',DurationColor,'Markersize',40,'linewidth',4);

hold on
% h1.FaceColor='none'

% plot((normperf'),'k');



h2=errorbar(SDI,nanmean(normperf),nanstd(normperf)./sqrt(size((normperf),1)))

h2.LineStyle='none';
h2.LineWidth=3;
h2.Color=DurationColor;
ylim([-13 13])
xlim([-0.4 0.4])

figurecorrectionOPTO(gca,Lenghts,Fontsize)

%%
normperf=[]
for i=1:size(perf1{1},1)

    %   normperf(i,:) =(perf1{2}(i,:)-perf1{2}(i,4))*100
    normperf(i,:) =(perf1{2}(i,:)-perf1{2}(i,4))*100

end



subplot(Axes{6})

h1=plot(TDI,nanmean(normperf),'color',IntensityColor,'linewidth',4);
hold on
% h1.FaceColor='none'

% plot((normperf'),'k');



h2=errorbar(TDI,nanmean(normperf),nanstd(normperf)./sqrt(size((normperf),1)));

h2.LineStyle='none';
h2.LineWidth=3;
h2.Color=IntensityColor;
ylim([-13 13])


figurecorrectionOPTO(gca,Lenghts,Fontsize)
%%
subplot(Axes{5})
xlabel('\Delta I')
ylabel('Bias in perceived duration')
legend ('Duration rats')

subplot(Axes{6})

ylabel('Bias in perceived intensity')
legend ('Intensity rats')
xlabel('\Delta T')
%%


subplot(Axes{2})
xlabel('\Delta I')
legend (['Duration rats '  ;'              ' ;'Intensity rats'] )

subplot(Axes{1})
xlabel('\Delta T')
legend (['Duration rats '  ;'              ' ;'Intensity rats'] )


