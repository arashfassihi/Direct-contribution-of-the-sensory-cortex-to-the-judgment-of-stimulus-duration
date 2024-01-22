subplotnumbers = [ 1 3 2 4];
colorcode = {'b-'; 'k-' ; 'k:'; 'b:'};
colormap hot
D = [161   205   264   334   422   545   694];
cd 'D:\Dropbox\OptoDurationPaper\A_NATCOMSUBMISSION\related to figure3'
load STIM1
load STIM2
stimulus = zeros(1,2334);
thisstim = randn(1,334);
stimperiod = (1001:1000+numel(thisstim))-20;
stimulus1 = stimulus;
stimulus1 (stimperiod)=(thisstim);
thisstim = randn(1,D(7));
stimshift = 20;
stimulus2 = stimulus;
stimperiod = (1001:1000+numel(thisstim))-20;
stimulus2 (stimperiod)=(thisstim);

% plot(abs(smooth(stimulus1,10)),'k')
for istim=[2 1]
    
    i = 7;
    counter = 0;
    a = [];
    c = [];
    for j=1:size(STIM1,2)
        
        %          for i=1:7
        %                       if istim==1
        %                           i =4;
        %                      end
        if ~isempty(STIM1{i,j})&&~isempty(STIM2{i,j})
            
            counter = counter+1;
            b = STIM1{i,j};
            a(counter,:) = mean(b.PSTH{istim});
            b = STIM2{i,j};
            c(counter,:) = mean(b.PSTH{istim});
            
        end
    end
    
    
    Clim = [0 4.5];
    beforestimperiod = [1:950 ];
    % beforestimperiod
    sortingperiod = 1000:1334;
    %     sortingperiod = 1000:1900;
    
    % subplot 211
    M1 = movmean(a',80)*1000;
    M2 = movmean(c',80)*1000;
      M1 = movmean(a',2)*1000;
    M2 = movmean(c',2)*1000;
    if istim ==2
        
        Mean = mean([M1(beforestimperiod,:); M2(beforestimperiod,:)],1);
        STD = std([M1(beforestimperiod,:); M2(beforestimperiod,:)],[],1);
    end
    Mz_exited = (M2-(repmat(Mean,size(M2,1),1)))./(repmat(STD,size(M2,1),1));
    Mz_Control =  (M1-(repmat(Mean,size(M1,1),1)))./(repmat(STD,size(M1,1),1));
    %             Mz_exited =M2;
    %            Mz_Control=M1;
    Period1 = 1:2334;
    Period2 = 1:2000;
    
    %[sortevalue ,sorted_Ind]= sort(mean(Mz_Control(sortingperiod,:)));
    
    figure(1)
    
    
    subplot (2,2,subplotnumbers(istim*2-1))
    colormap jet
    
    if istim ==2
        [~, sorted_Ind]= sort(mean(Mz_exited(sortingperiod,:)));
        %  [~, sorted_Ind]= sort(mean(Mz_exited(sortingperiod,:)+Mz_Control(sortingperiod,:)));
        
    end
    %     [~, sorted_Ind]= sort(mean(Mz_exited(sortingperiod,:)));
    if istim ==2
        
        Period = Period2;
    elseif istim ==1
        
        Period  = Period1;
    end
    XX{istim}=Mz_exited(Period,sorted_Ind(numel(sorted_Ind):-1:1))';
    XX2{istim} = Mz_Control(Period,sorted_Ind(numel(sorted_Ind):-1:1))';
    imagesc(Mz_exited(Period,sorted_Ind(numel(sorted_Ind):-1:1))',Clim)
    xlim([600 2705])
    subplot (2,2,subplotnumbers(istim*2))
    colormap jet
    
    [sortevalue ,sorted_Ind]= sort(mean(Mz_Control(sortingperiod,:)));
    imagesc(Mz_Control(Period,sorted_Ind(numel(sorted_Ind):-1:1))',Clim)
    
    xlim([600 2705])
    figure(2)
    plot(nanmean(Mz_Control(Period,sorted_Ind(numel(sorted_Ind):-1:1))'),colorcode{(istim*2-1)})
    hold on
    plot(nanmean(Mz_exited(Period,sorted_Ind(numel(sorted_Ind):-1:1))'),colorcode{(istim*2)})
    figure(1)
    subplot (2,2,(istim*2-1))
    hold on
    ylim([-40 300])
    set (gca,'Ydir','reverse')
    plot(-stimshift-abs(smooth(stimulus1,10))*15,'k','Linewidth',2)
    xlim([600 2705])
    ylim([-40 250])
    % axis off
    subplot (2,2,(istim*2))
    hold on
    ylim([-40 250])
    set (gca,'Ydir','reverse')
    
    plot(-stimshift-abs(smooth(stimulus2,10))*15,'k','Linewidth',2)
    % ylim([-10 5])
    xlim([600 2705])
    ylim([-40 250])
        
    
end

%%

figure
subplot (2,1,1)
Clim = [0 4];

imagesc([XX2{1} XX{2}],Clim)
colormap jet
xlim([600 4200])
ylim([-40 320])
subplot (2,1,2)
Clim = [0 6];
imagesc([XX{1} XX2{2}],Clim)
colormap jet
xlim([600 4200])
ylim([-40 320])