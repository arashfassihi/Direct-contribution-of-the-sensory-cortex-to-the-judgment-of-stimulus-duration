% Stimulusonset
loadthis=1
if loadthis==1
    cd 'D:\Dropbox\OptoDurationPaper\A_NATCOMSUBMISSION\related to figure3'
    load STIM1
    load STIM2
end
subplotnumbers = [ 1 3 2 4];
close all
colorcode = {'b-'; 'k-' ; 'k:'; 'b:'};
colormap hot
D = [    161   205   264   334   422   545   694  ];
X = -999:1000;
stimulus = zeros(1,2001);
thisstim = randn(1,334);
stimulus1 = stimulus;
stimulus1 (1001:1000+numel(thisstim))=(thisstim);
thisstim = randn(1,D(7));
colorcodeduration = jet(7);
stimulus2 = stimulus;
stimulus2 (1001:1000+numel(thisstim))=(thisstim);
% plot(abs(smooth(stimulus1,10)),'k')
%%



for istim=2
    for i = 1:7;
        counter = 0;
        a = [];
        c = [];
        spn=[];
        
        M1=[];
        M2=[];
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
                spn(counter) = mean(b.ratesp{istim});

            end
        end
        
        
        Clim = [0 4.5];
        beforestimperiod = [1:950 ];
        % beforestimperiod
        sortingperiod = 1000:1334;
        % subplot 211
        Thissigma=15;
        pd1 = makedist('HalfNormal','mu',0,'sigma',25);
        x = 1:200;
        pdf1 = pdf(pd1,x);
        % plot(pdf1)
        M2=[];
        M1=[];
        for xxx=1:size(a,1)
            Thiy = conv2(a(xxx,1:1000),pdf1(end:-1:1),'full')'*1000./sum(pdf1);
            M1.bef(:,xxx) = Thiy(size(pdf1,2):end-Thissigma*3);
            Thiy = conv2(c(xxx,1:1000),pdf1(end:-1:1),'full')'*1000./sum(pdf1);
            M2.bef(:,xxx) = Thiy(size(pdf1,2):end-Thissigma*3);
            
            %             Thiy = a(xxx,1:1000);
            %             M1.bef(:,xxx) =Thiy;
            %             Thiy = c(xxx,1:1000);
            %             M2.bef(:,xxx) = Thiy;
            %
            
            
            M1.after(:,xxx) = conv2(a(xxx,:),pdf1)'*1000./sum(pdf1);
            M2.after(:,xxx) = conv2(c(xxx,:),pdf1)'*1000./sum(pdf1);
            Thiy = conv2(a(xxx,900:end),pdf1,'full')'*1000./sum(pdf1);
            %
            M1.during(:,xxx) = Thiy(Thissigma*3+50:end-size(pdf1,2));
            Thiy = conv2(c(xxx,900:end),pdf1,'full')'*1000./sum(pdf1);
            M2.during(:,xxx) = Thiy(Thissigma*3+50:end-size(pdf1,2));
            
%             Thiy = a(xxx,900:end);
%             M1.during(:,xxx) = Thiy;
%             Thiy = c(xxx,900:end);
%             M2.during(:,xxx) = Thiy;
            
        end
        %    plot(nanmean(M2.during,2))
        
        
        Mz_exited =M2;
        Mz_Control=M1;
        
        
        
        
        
        
        %         [sortevalue ,sorted_Ind]= sort(mean(Mz_Control(sortingperiod,:)));
        sorted_Ind=find(  spn>0.15);
        
        subplot(1,4,1)
        X = (1:numel(nanmean(Mz_Control.bef(:,sorted_Ind),2)))-size(Mz_Control.bef,1);
        plot(X,nanmean(Mz_Control.bef(:,sorted_Ind),2),'color',colorcodeduration(i,:))
        hold on
        X = (1:numel(nanmean(Mz_Control.during(:,sorted_Ind),2)));
        plot(X,nanmean(Mz_Control.during(:,sorted_Ind),2),'color',colorcodeduration(i,:))
        subplot(1,4,3)
        X = (1:numel(nanmean(Mz_Control.bef(:,sorted_Ind),2)))-size(Mz_Control.bef,1);
        plot(X,nanmean(Mz_exited.bef(:,sorted_Ind),2),'color',colorcodeduration(i,:))
        hold on
        X = (1:numel(nanmean(Mz_Control.during(:,sorted_Ind),2)));
        plot(X,nanmean(Mz_exited.during(:,sorted_Ind),2),'color',colorcodeduration(i,:))
        hold on
        
        
    end
end
%%

for i=1:2
    subplot(1,4,(i*2)-1)
    box off
    
    figurecorrectionOPTO(gca,[0.01 0.01],12)
    xlim([-100 D(1)])
    ylim([22 55])
end

%% Stimulus offset
% Stimulusonset
subplotnumbers = [ 1 3 2 4];
colorcode = {'b-'; 'k-' ; 'k:'; 'b:'};
colormap hot
D = [    161   205   264   334   422   545   694  ];
X = -999:1000;
stimulus = zeros(1,2001);
thisstim = randn(1,334);
stimulus1 = stimulus;
stimulus1 (1001:1000+numel(thisstim))=(thisstim);
thisstim = randn(1,D(7));
colorcodeduration = jet(7);
stimulus2 = stimulus;
stimulus2 (1001:1000+numel(thisstim))=(thisstim);
% plot(abs(smooth(stimulus1,10)),'k')

for istim=2
    for i = 1:7;
        counter = 0;
        a = [];
        c = [];
        spn=[];
        
        M1=[];
        M2=[];
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
                spn(counter) = mean(b.ratesp{istim});
                
                %     end
            end
        end
        
        
        Clim = [0 4.5];
        beforestimperiod = [1:950 ];
        % beforestimperiod
        sortingperiod = 1000:1334;
        % subplot 211
        pd1 = makedist('HalfNormal','mu',0,'sigma',30);
        x = 1:600;
        pdf1 = pdf(pd1,x);
        % plot(pdf1)
        M2=[];
        M1=[];
        for xxx=1:size(a,1)
            
            Thiy = conv2(a(xxx,1000+D(i)-800:1000+D(i)+400),pdf1(end:-1:1))'*1000./sum(pdf1);
            M1.during(:,xxx) = Thiy(size(pdf1,2)+50:end-400-Thissigma*3);
            Thiy = conv2(c(xxx,1000+D(i)-800:1000+D(i)+400),pdf1(end:-1:1))'*1000./sum(pdf1);
            M2.during(:,xxx) = Thiy(size(pdf1,2)+50:end-400-Thissigma*3);
            %             M1.after(:,xxx) = conv2(a(xxx,:),pdf1)'*1000./sum(pdf1);
            %             M2.after(:,xxx) = conv2(c(xxx,:),pdf1)'*1000./sum(pdf1);
            Thiy = conv2(a(xxx,1000+D(i):end),pdf1,'full')'*1000./sum(pdf1);
            %             M1.after(:,xxx) = Thiy(100:end-size(pdf1,2));
            M1.after(:,xxx) = Thiy(Thissigma*3+50:end-size(pdf1,2));
            
            Thiy = conv2(c(xxx,1000+D(i):end),pdf1,'full')'*1000./sum(pdf1);
            %             M2.after(:,xxx) = Thiy(100:end-size(pdf1,2));
            M2.after(:,xxx) = Thiy(Thissigma*3+50:end-size(pdf1,2));
            
            %             er(:,xxx) = Thiy(90:end-size(pdf1,2));
        end
        
        
        Mz_exited =M2;
        Mz_Control=M1;
        
        
        
        
        
        
        %         [sortevalue ,sorted_Ind]= sort(mean(Mz_Control(sortingperiod,:)));
        sorted_Ind=find(  spn>0.15);
        
        subplot(1,4,2)
        X = (1:numel(nanmean(Mz_Control.during(:,sorted_Ind),2)))-size(Mz_Control.during,1);
        plot(X,nanmean(Mz_Control.during(:,sorted_Ind),2),'color',colorcodeduration(i,:))
        hold on
        X = (1:numel(nanmean(Mz_Control.after(:,sorted_Ind),2)));
        plot(X,nanmean(Mz_Control.after(:,sorted_Ind),2),'color',colorcodeduration(i,:))
        subplot(1,4,4)
        X = (1:numel(nanmean(Mz_exited.during(:,sorted_Ind),2)))-size(Mz_exited.during,1);
        plot(X,nanmean(Mz_exited.during(:,sorted_Ind),2),'color',colorcodeduration(i,:))
        hold on
        X = (1:numel(nanmean(Mz_exited.after(:,sorted_Ind),2)));
        plot(X,nanmean(Mz_exited.after(:,sorted_Ind),2),'color',colorcodeduration(i,:))
        hold on
        
        
    end
end
%%

for i=1:2
    subplot(1,4,i*2)
    box off
    
    figurecorrectionOPTO(gca,[0.01 0.01],12)
    xlim([-D(i)+100 100])
    ylim([22 55])
end
%%



    subplot(1,4,1)
    ylabel('Response rate')

    xlabel('time from stim onset ')

    
       subplot(1,4,2)

    xlabel('time from stim offset ')
    
        subplot(1,4,1)
        title('Light off')
        
        
            subplot(1,4,3)
        title('Light on ')
    

