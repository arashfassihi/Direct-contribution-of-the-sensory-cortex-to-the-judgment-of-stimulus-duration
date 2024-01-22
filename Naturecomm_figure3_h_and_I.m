% clear all
dothereampling = 1;
if dothereampling==1
    subplotnumbers = [ 1 3 2 4];
    colorcode = {'b-'; 'k-' ; 'k:'; 'b:'};
    loadagain = 1;
    if loadagain==1
        %      cd 'D:\Dropbox\Work_related\Duration opto\Paper figures\figure3'
        cd 'D:\Dropbox\OptoDurationPaper\A_NATCOMSUBMISSION\related to figure3'

        load STIM1
        load STIM2
    end
    thissize = [2172        2216        2275        2345        2433        2556        2705];
    D = [    161   205   264   334   422   545   694  ];
    stimulus = zeros(1,2001);
    thisstim = randn(1,334);
    stimulus1 = stimulus;
    stimulus1 (1001:1000+numel(thisstim))=(thisstim);
    thisstim = randn(1,D(7));
    Smoothingsize =10;
    stimulus2 = stimulus;
    stimulus2 (1001:1000+numel(thisstim))=(thisstim);
    % plot(abs(smooth(stimulus1,10)),'k')
    thiscolormap = jet;
    thiscolormap(1,:)=[ 1 1 1];
    for istim=1:2

        for  i = 1:7
            counter = 0;
            a = NaN(1 ,thissize(i)+1);
            c = NaN(1 ,thissize(i)+1);

            figure(i)

            for j=1:size(STIM1,2)

                %% add the first line of each matrix to be the  neuron number
                if ~isempty(STIM1{i,j})&&~isempty(STIM2{i,j})

                    counter = counter+1;
                    b = STIM1{i,j};
                    %                  a= (b.PSTH{istim});

                    a =[a ; [b.PSTH{istim},ones(size(b.PSTH{istim},1),1)*j ]];
                    b = STIM2{i,j};
                    c =[c ; [b.PSTH{istim},ones(size(b.PSTH{istim},1),1)*j ]];
                    %     end
                end
            end
            Clim = [0 94.5]/10;
            Clim = [0 1500];
            beforestimperiod = [1:950 ];
            % beforestimperiod
            sortingperiod = 1000:1334;
            % subplot 211
            %         M1 = movmean(a',Smoothingsize)*1000;
            %         M2 = movmean(c',Smoothingsize)*1000;
            %         M1 = a;
            %         M2 = c;
            %          figure(3)
            %      plot(var(M1'))
            % hold on
            %           plot(var(M2'))
            %         Mean = mean([M1(beforestimperiod,:); M2(beforestimperiod,:)],1);
            %         STD = std([M1(beforestimperiod,:); M2(beforestimperiod,:)],[],1);
            %         Mz_exited = (M2-(repmat(Mean,size(M2,1),1)))./(repmat(STD,size(M2,1),1));
            %         Mz_Control =  (M1-(repmat(Mean,size(M1,1),1)))./(repmat(STD,size(M1,1),1));
            Mz_exited =c';
            Mz_Control=a';
            Period = 1:2000;
            %    [sortevalue ,sorted_Ind]= sort(mean(Mz_Control(sortingperiod,:)));
            subplot (2,2,subplotnumbers(istim*2-1))
            colormap(thiscolormap)
%             if istim ==1
%                 [~, sorted_Ind]= sort(nanmean(Mz_exited(sortingperiod,:)));
%             end
            [~, sorted_Ind]= sort(nanmean(Mz_exited(sortingperiod,:)));

            imagesc(Mz_exited(Period,sorted_Ind(numel(sorted_Ind):-1:1))',Clim)
            xlim([600 2000])

            subplot (2,2,subplotnumbers(istim*2))
            colormap(thiscolormap)

            [sortevalue ,sorted_Ind]= sort(mean(Mz_Control(sortingperiod,:)));
            if istim ==1

                thismean1{i}.mzcontrol=(Mz_Control);
                thismean1{i}.Mz_exited=(Mz_exited);

                thismean1raw{i}.mzcontrol=(Mz_Control);
                thismean1raw{i}.Mz_exited=(Mz_exited);
            else
                thismean2{i}.mzcontrol=(Mz_Control);
                thismean2{i}.Mz_exited=(Mz_exited);
                thismean2raw{i}.mzcontrol=(Mz_Control);
                thismean2raw{i}.Mz_exited=(Mz_exited);
            end


            imagesc(Mz_Control(Period,sorted_Ind(numel(sorted_Ind):-1:1))',Clim)

            xlim([600 2000])

            subplot (2,2,(istim*2-1))
            hold on
            ylim([-40 300])
            set (gca,'Ydir','reverse')
            plot(-20-abs(smooth(stimulus1,10))*15,'k','Linewidth',2)
            xlim([600 2000])
            ylim([-40 250])
            % axis off
            subplot (2,2,(istim*2))
            hold on
            ylim([-40 250])
            set (gca,'Ydir','reverse')



            plot(-20-abs(smooth(stimulus2,10))*15,'k','Linewidth',2)
            % ylim([-10 5])

            xlim([600 2000])
            ylim([-40 250])

            %         axis off
        end


    end
    thismean2raw{i}.mzcontrol(1,:)=[];
    thismean2raw{i}.Mz_exited(1,:)=[];
end
%%
% thismean2raw{i}.mzcontrol(1,:)=[];
% thismean2raw{i}.Mz_exited(1,:)=[];
figure
hold off
clear D2

Threshold=2.2;

D = [    161   205   264   334   422   545   694  ];
thiscolor=[0 0.4 1;0 0 0; 0 1 0; 0.5 0.5 0.5; 0.5 0.5 1; 0 0.5 0.5];
stimonset = 995;
fixdelay = 0;
minsize=5;
plotting= 0;
Smoothingsize = 15;


% here i  should write sampling from neurons id =thismean2raw{i}.mzcontrol(end,:) and then withing this subset of neurons then resample trials other wise I apply same threshold to differnet neuronal subset and compare them

plotting = 0;

for j=1:1000
    tic
    for i=1:7
%         if j==1
            %     subplot(7,1,i)
            M1 = movmean(thismean2raw{i}.mzcontrol(1:end-1,:),Smoothingsize)'*1000;
            M2 = movmean(thismean2raw{i}.Mz_exited(1:end-1,:),Smoothingsize)'*1000;
            M1all = [M1 nan(size(M1,1),3000-size(M1,2))];
            M2all = [M2 nan(size(M2,1),3000-size(M2,2))];
            subgroupofneuronsall = [thismean2raw{i}.mzcontrol(end,:)];
            if i<7
                for jjj = i+1:7
                    thisM1 = movmean(thismean2raw{jjj}.mzcontrol(1:end-1,:),Smoothingsize)'*1000;
                    thisM2 = movmean(thismean2raw{jjj}.Mz_exited(1:end-1,:),Smoothingsize)'*1000;
                    subgroupofneuronsall= [subgroupofneuronsall thismean2raw{jjj}.mzcontrol(end,:)];
                    thisM1 = [thisM1 nan(size(thisM1,1),3000-size(thisM1,2))];
                    thisM2 = [thisM2 nan(size(thisM2,1),3000-size(thisM2,2))];
                    M1all =[M1all ;thisM1];
                    M2all =[M2all ;thisM2];
                end
            end
%         end
        % [cnt_unique, unique_a] = hist(thismean2raw{i}.mzcontrol(:,end),unique(thismean2raw{i}.mzcontrol(:,end)))
        % thisunique =
        %  G = findgroups(thismean2raw{i}.mzcontrol(:,end)(ismember(thismean2raw{i}.mzcontrol(:,end),thisunique)));

        RandC = randi( size(M1,1),1, size(M1,1));
        subgroupofneurons = thismean2raw{i}.mzcontrol(end,:);
        [theseneurons,~] = unique(subgroupofneurons(RandC));
        theseneurons=theseneurons(~isnan(theseneurons));
        Neuroncount = hist(subgroupofneurons(RandC),theseneurons);
        AllrandE=[];


        for jj=1:numel(Neuroncount)

            RandEdummy = find(thismean2raw{i}.Mz_exited(end,:)==theseneurons(jj));
            thesetrials = randi(numel(RandEdummy),1,Neuroncount(jj));
            AllrandE =[AllrandE RandEdummy(thesetrials)];


        end
        % randomzing for spikecount population

        RandCSpike = randi( numel(subgroupofneuronsall),1, numel(subgroupofneuronsall));

        [theseneurons,~] = unique(subgroupofneuronsall(RandCSpike));
        theseneurons=theseneurons(~isnan(theseneurons));
        Neuroncount = hist(subgroupofneuronsall(RandCSpike),theseneurons);
        AllrandEspikes=[];

        for jj=1:numel(Neuroncount)
            RandEdummy=[];
            for jjj=i:7

                RandEdummy =[RandEdummy find(thismean2raw{jjj}.Mz_exited(end,:)==theseneurons(jj))];

            end


            thesetrials = randi(numel(RandEdummy),1,Neuroncount(jj));
            AllrandEspikes =[AllrandEspikes RandEdummy(thesetrials)];


        end



        X = nanmean(M1(RandC,:));
        X2 = nanmean(M1(RandC,:));
        Y = nanmean(M2(AllrandE,:));
        Y2 = nanmean(M2(AllrandE,:));


        RandE=AllrandE;



        X=X(500:1000+D(i)+200);
        Y=Y(500:1000+D(i)+200);
        if plotting==1
            figure(j)

            subplot(7,1,i)

            timevalues= (1:numel( nanmean(M1)))-stimonset;
            plot(timevalues, nanmean(M1(RandC,:)),'k','linewidth',2)
            hold on
            plot(timevalues, nanmean(M2(RandE,:)),'color',thiscolor(1,:),'linewidth',2)

        end
        thistrheshold(i,j)= ((mean(X2(500:980))+Threshold*std(X2(500:980))+mean(Y2(500:980))+Threshold*std(Y2(500:980))))./2;
        Actualthreshod = thistrheshold(i,j);

        [a,b]=find(X>Actualthreshod);
        b=findcontinuessegments(b,minsize);
        XX = zeros(size(X));
        XX(b)=1;
        b=find(smooth(XX,20)>0.9);
        if plotting==1
            figure(j)
            subplot(7,1,i)
            plot([1 2000]-stimonset,ones(1,2)*(Actualthreshod),'k--')

            plot(b+500-stimonset,X(b),'.g')
            plot([b(1) b(1)]+500-stimonset,[0 max(X)],'k--')
            plot([b(end) b(end)]+500-stimonset,[0 max(X)],'k--')
        end
        D2.C(i,j) = b(end)-b(1);
        D2.COnset(i,j) = b(1);
        D2.COffset(i,j) = b(end);
        %         ThisMat =M1(RandC,:);
        %         ThisMat= M1all;
        ThisMat= M1all(RandCSpike,:);

        D2.MeanC(i,j)=nanmean(nansum(ThisMat(:,b(1)+500-fixdelay:b(end)+500-fixdelay)'));
        D2.firingrateC(i,j)=nanmean(nansum(ThisMat(:,1000-fixdelay:1000+D(i)-fixdelay,:)'))./D(i);

        %     D2.MeanC(i)=nanmean(nansum(thismean2raw{i}.mzcontrol(1000-fixdelay:D(i)+1000-fixdelay,:)));

        [a,b]=find(Y>Actualthreshod);
        b=findcontinuessegments(b,minsize);
        XX = zeros(size(Y));
        XX(b)=1;
        b=find(smooth(XX,20)>0.9);
        if plotting==1
            figure(j)
            subplot(7,1,i)

            plot(b+500-stimonset,Y(b),'.r')
            plot([b(1) b(1)]+500-stimonset,[0 max(Y)],'color',thiscolor(1,:))
            plot([b(end) b(end)]+500-stimonset,[0 max(Y)],'color',thiscolor(1,:))
            xlim([-500 1000])
        end




        D2.E(i,j) = b(end)-b(1);
        %         D2.MeanE(i,j) = nanmean(nansum(thismean2raw{i}.Mz_exited(b(1)+500-fixdelay:b(end)+500-fixdelay,:)));
        %         D2.firingrateE(i,j)=nanmean(nansum(thismean2raw{i}.Mz_exited(1000-fixdelay:1000+D(i)-fixdelay,:)))./D(i);

        %         ThisMat =M2(AllrandE,:);
        ThisMat= M2all(AllrandEspikes,:);

        D2.MeanE(i,j)=nanmean(nansum(ThisMat(:,b(1)+500-fixdelay:b(end)+500-fixdelay)'));
        % here I should use the true onset ofset
        D2.firingrateE(i,j)=nanmean(nansum(ThisMat(:,1000-fixdelay:1000+D(i)-fixdelay,:)'))./D(i);

        D2.EOnset(i,j) = b(1);
        D2.EOffset(i,j) = b(end);
        %     D2.MeanE(i)=nanmean(nansum(thismean2raw{i}.Mz_exited(1000-fixdelay:D(i)+1000-fixdelay,:)));
    end
    if plotting==1

        saveas(figure(j),[num2str(j) '.jpeg']);
        saveas(figure(j),[num2str(j) '.fig']);
        close all
    end
    toc
end
%% make figure 3 h
% repmat(D,1,1000)
figure
thiscolor = jet(7);
for i=1:7

    thesevalidpoints_onset = abs(D2.EOnset(i,:)-490)<25&abs(D2.COnset(i,:)-490)<25;
    thesevalidpoints_offset = abs(D2.EOffset(i,:)-490-D(i))<25&abs(D2.COffset(i,:)-490-D(i))<25;

    %   thesevalidpoints_onset=thesevalidpoints_offset
    %   thesevalidpoints_onset = abs(D2.EOnset(i,:)-490)<40&abs(D2.COnset(i,:)-490)<40;
    %   thesevalidpoints_offset = abs(D2.EOffset-490)<40;

    thisdits = D2.E(i,:)-D2.C(i,:);
    thisdits1 = D2.E(i,:);
    thisdits2 = D2.C(i,:);
    TF = thesevalidpoints_offset&thesevalidpoints_onset &~isoutlier(thisdits,'quartiles')&~isoutlier(thisdits2,'quartiles')&~isoutlier(thisdits1,'quartiles');

    % TF = isoutlier(thisdits1,'quartiles')|isoutlier(thisdits2,'quartiles');
    %     plot(D2.C(i,~TF),D2.E(i,~TF),'.','markersize',10,'color',[0.3 0.3 0.3])
    %     plot(D2.C(i,~TF),D2.E(i,~TF),'.','markersize',10,'color',thiscolor(i,:))
    thisC=D2.COffset(i,TF)-D2.COnset(i,TF);
    thisE=D2.EOffset(i,TF)-D2.EOnset(i,TF);

    plot(thisC,thisE,'.','markersize',10,'color',thiscolor(i,:))


    Tvalues(i,:) = (D2.E(i,:)-D2.C(i,:));
    normTvalues(i,:) = (D2.E(i,:)-D2.C(i,:))./(D2.E(i,:)+D2.C(i,:));
    normTvalues(i,~TF)=NaN;
    Tvaluesall(i)=sum(Tvalues(~TF)>0)./sum(~TF);

    alltfs(i,:)=TF;
    hold on
    %     plot(mean(D2.C(i,~TF)),mean(D2.E(i,~TF)),'.','markersize',30,'color',thiscolor(i,:))

    %     plot([D(i) D(i)] ,[D(1) D(end)],'--','markersize',30,'color',thiscolor(i,:))
    %     plot([D(1) D(end)],[D(i) D(i)] ,'--','markersize',30,'color',thiscolor(i,:))

end
plot([120 750],[120 750],'k')
xlim([120 750])
ylim([120 750])
xlabel(['EST. duration' '\_' 'light off (ms)'])
ylabel(['EST. duration' '\_' 'light on (ms)'])
figurecorrectionOPTO(gca,[0.05, 0.05],12)


set(gca, 'YScale', 'log')
set(gca, 'xScale', 'log')
%%
figure

% TableDuration = struct2table(D2.EOnset');
% D2.COnset'
Dtable.Duration = D2.EOnset(:);
Dtable.Duration(alltfs(:)) =NaN;
Dtable.thisx= (repmat(1:7,1,1000))';
Dtable.Lighton= ones(size(D2.EOnset(:)));

Y =D2.COnset(:);
Y(alltfs(:)) = NaN;

Dtable.Duration = [Dtable.Duration;Y];
Dtable.thisx=[Dtable.thisx ;(repmat(1:7,1,1000))'];
Dtable.Lighton=[Dtable.Lighton;zeros(size(D2.COnset(:)))];


boxchart(Dtable.thisx,Dtable.Duration-490 ,'GroupByColor',Dtable.Lighton)
% boxchart(D2.COnset')
hold on
%%

plot([1 9],[-40 -40],'k--')
plot([1 9],[40 40],'k--')

ylim([-85 85])


%%
% D2.COnset'
Dtable.Duration = D2.EOffset(:);
Dtable.Duration(alltfs(:)) =NaN;
Dtable.thisx= (repmat(1:7,1,1000))';
Dtable.Lighton= ones(size(D2.EOffset(:)));

Y =D2.COffset(:);
Y(alltfs(:)) = NaN;

Dtable.Duration = [Dtable.Duration;Y];
Dtable.thisx=[Dtable.thisx ;(repmat(1:7,1,1000))'];
Dtable.Lighton=[Dtable.Lighton;zeros(size(D2.COffset(:)))];




boxchart(Dtable.thisx,Dtable.Duration-480 ,'GroupByColor',Dtable.Lighton)


%%
figure

% TableDuration = struct2table(D2.EOnset');
% D2.COnset'
Dtable.Duration = D2.EOffset(:);

% Dtable.Duration(alltfs(:)) =NaN;
Dtable.thisx= (repmat(1:7,1,1000))';
Dtable.thisD= (repmat(D,1,1000))';

Dtable.Lighton= ones(size(D2.EOffset(:)));

Y =D2.COffset(:);
Y(alltfs(:)) = NaN;

Dtable.Duration = [Dtable.Duration;Y];
Dtable.thisx=[Dtable.thisx ;(repmat(1:7,1,1000))'];
Dtable.Lighton=[Dtable.Lighton;zeros(size(D2.COffset(:)))];
Dtable.thisD= [Dtable.thisD ;(repmat(D,1,1000))'];


boxchart(Dtable.thisx,Dtable.Duration-Dtable.thisD-490 ,'GroupByColor',Dtable.Lighton)
% boxchart(D2.COnset')
hold on

%%

figure
thiscolor = jet(7);

for i=1:7

    thesevalidpoints_onset = abs(D2.EOnset(i,:)-490)<25&abs(D2.COnset(i,:)-490)<25;
    thesevalidpoints_offset = abs(D2.EOffset(i,:)-490-D(i))<25&abs(D2.COffset(i,:)-490-D(i))<25;
    thisdits = D2.MeanC(i,:)-D2.MeanE(i,:);

    TF = thesevalidpoints_offset&thesevalidpoints_onset &~isoutlier(thisdits,'quartiles');
    %     TF = isoutlier(thisdits);
    plot(D2.MeanC(i,TF)/1000,D2.MeanE(i,TF)/1000,'.','markersize',10,'color',thiscolor(i,:))
    hold on

    normTvaluesSpikecount(i,:) = (D2.MeanE(i,:)-D2.MeanC(i,:))./(D2.MeanE(i,:)+D2.MeanC(i,:));
    normTvaluesSpikecount(i,~TF)=NaN;
    TvaluesSpikecount(i,:) = (D2.MeanE(i,:)-D2.MeanC(i,:));
    TvaluesSpikecount(i,~TF)=NaN;

    %     plot(mean(D2.MeanC(i,~TF)),mean(D2.MeanE(i,~TF)),'.','markersize',30,'color',thiscolor(i,:))

end
plot([4 27],[4 27],'k')
xlabel('light off(spikes)')
ylabel('light on(spikes)')
set(gca, 'YScale', 'log')
set(gca, 'xScale', 'log')
figurecorrectionOPTO(gca,[0.05, 0.05],12)

%%
% plot([1 9],[-40 -40],'k--')
% plot([1 9],[40 40],'k--')

% ylim([-85 85])
%%
thisesge=-0.3:0.03:0.3;

[a,b] = histc(normTvalues',thisesge);
thisorder = [ 1 3 5 7 9 11 13 ]
thisorder = [ thisorder thisorder+1 ];

figure
for i=1:7
    subplot(7,2,thisorder(i))

    stairs(thisesge,a(:,i),'color',thiscolor(i,:))

    if i==7
        figurecorrectionOPTO(gca,[0.05, 0.05],12)

    else

        figurecorrectionOPTO(gca,[0.0, 0.0],1)

    end
    axis normal
    ylim([1 1000])
end


[a,b]=histc(normTvaluesSpikecount',thisesge);


for i=1:7
    subplot(7,2,thisorder(i+7))

    stairs(thisesge,a(:,i),'color',thiscolor(i,:))

    if i==7
        figurecorrectionOPTO(gca,[0.05, 0.05],12)

    else

        figurecorrectionOPTO(gca,[0.0, 0.0],1)

    end
    axis normal
    ylim([1 1000])
end


%%



figure
thiscolor = jet(7);

for i=1:7

    thesevalidpoints_onset = abs(D2.EOnset(i,:)-490)<25&abs(D2.COnset(i,:)-490)<25;
    thesevalidpoints_offset = abs(D2.EOffset(i,:)-490-D(i))<25&abs(D2.COffset(i,:)-490-D(i))<25;
    thisdits = D2.MeanC(i,:)-D2.MeanE(i,:);

    TF = thesevalidpoints_offset&thesevalidpoints_onset &~isoutlier(thisdits,'quartiles');
    %     TF = isoutlier(thisdits);
    plot(D2.MeanC(i,TF),D2.MeanE(i,TF),'.','markersize',10,'color',thiscolor(i,:))
    hold on
    normTvaluesSpikecount(i,:) = (D2.MeanE(i,:)-D2.MeanC(i,:))./(D2.MeanE(i,:)+D2.MeanC(i,:));
    normTvaluesSpikecount(i,~TF)=NaN;




    %     plot(mean(D2.MeanC(i,~TF)),mean(D2.MeanE(i,~TF)),'.','markersize',30,'color',thiscolor(i,:))

end
plot([4 27],[4 27],'k')
xlabel('light off(spikes)')
ylabel('light on(spikes)')

%%
plot([1 9],[-40 -40],'k--')
plot([1 9],[40 40],'k--')

ylim([-85 85])
%%
thisesge=-0.3:0.02:0.3;

[a,b]=histc(normTvalues',thisesge);
thisorder = [ 1 3 5 7 9 11 13 ]
thisorder = [ thisorder thisorder+1 ];

figure
for i=1:7
    subplot(7,1,(i))

    stairs(thisesge,a(:,i),'color',[ 0 0 0])

    if i==7
        figurecorrectionOPTO(gca,[0.05, 0.05],12)

    else

        figurecorrectionOPTO(gca,[0.0, 0.0],1)

    end
    axis normal
    ylim([1 1000])
end


[a,b]=histc(normTvaluesSpikecount',thisesge);


for i=1:7
    subplot(7,1,(i))

    h=stairs(thisesge,a(:,i),'color',thiscolor(i,:))
    h.LineStyle='-';
    if i==7
        figurecorrectionOPTO(gca,[0.05, 0.05],12)

    else

        figurecorrectionOPTO(gca,[0.0, 0.0],1)

    end
    axis normal
    ylim([1 1000])
end


%%
thisesge=-0.3:0.02:0.3;

[a,b] = histc(normTvalues',thisesge);
thisorder = [ 1 3 5 7 9 11 13 ]
thisorder = [ thisorder thisorder+1 ];

figure

stairs(thisesge,sum(a'),'color',[ 0 0 0])
figurecorrectionOPTO(gca,[0.0, 0.0],1)



[a,b]=histc(normTvaluesSpikecount',thisesge);
stairs(thisesge,sum(a'),'color',[ 0 0 1])



%%
figure
% TvaluesSpikecount
% Tvaluesall
thesedurations= [1:2 4:7]
Y =normTvaluesSpikecount (thesedurations,:);
Dtable.Duration =Y(:);
% Dtable.Duration(alltfs(:)) =NaN;
Dtable.thisx= (repmat(1:size(Y,1),1,1000))';

Dtable.Lighton= ones(size(Y(:)));

Y =normTvalues(thesedurations,:);

Dtable.Duration = [Dtable.Duration;Y(:)];
Dtable.thisx=[Dtable.thisx ;(repmat(1:size(Y,1),1,1000))'];
Dtable.Lighton=[Dtable.Lighton;zeros(size(Y(:)))];
Dtable.thisx= ones(size(Dtable.Duration))';



boxchart(Dtable.thisx,Dtable.Duration ,'GroupByColor',Dtable.Lighton)

%%
figure
errorbar(1,nanmean(normTvaluesSpikecount(:)),nanstd(normTvaluesSpikecount(:)))
hold on
h=errorbar(1-0.1,nanmean(normTvalues(:)),nanstd(normTvalues(:)))


% boxchart(D2.COnset')
%% p values

% each duration
1-(sum( normTvaluesSpikecount' >0)./sum(~isnan(normTvaluesSpikecount')))
% total
Y =normTvaluesSpikecount ([1:2 3 4:7],:);

text(1,0.19,['pvalue=' num2str(1-(sum(Y(:)>0)./sum(~isnan(Y(:)))))])
%%%%%%
%  for onset offset

1-(sum( normTvalues' >0)./sum(~isnan(normTvalues')))
% total

1-(sum( normTvalues(:) >0)./sum(~isnan(normTvalues(:))))


text(0.6,0.19,['pvalue=' num2str(1-(sum( normTvalues(:) >0)./sum(~isnan(normTvalues(:)))))])
%%

figure
thiscolor = jet(7);

for i=1:7

    thesevalidpoints_onset = abs(D2.EOnset(i,:)-490)<25&abs(D2.COnset(i,:)-490)<25;
    thesevalidpoints_offset = abs(D2.EOffset(i,:)-490-D(i))<25&abs(D2.COffset(i,:)-490-D(i))<25;
    thisdits = D2.firingrateC(i,:)-D2.firingrateE(i,:);
    thisdits1 = D2.C(i,:)-D2.E(i,:);

    TF = thesevalidpoints_offset&thesevalidpoints_onset &~isoutlier(thisdits,'quartiles')&~isoutlier(thisdits1,'quartiles');
    TF = thesevalidpoints_offset&thesevalidpoints_onset &~isoutlier(thisdits,'quartiles');

    %     TF = isoutlier(thisdits);

    thisC=D2.COffset(i,TF)-D2.COnset(i,TF);
    thisE=D2.EOffset(i,TF)-D2.EOnset(i,TF);

    plot(D2.firingrateC(i,TF),D2.firingrateE(i,TF),'.','markersize',10,'color',thiscolor(i,:))
    hold on
    %     normTvaluesSpikecount(i,:) = (D2.MeanE(i,:)-D2.MeanC(i,:))./(D2.MeanE(i,:)+D2.MeanC(i,:));


end
plot([0.031 0.044],[0.031 0.044],'k')

xlabel('light off(spike rate hz)')

ylabel('light on(spike rate hz)')