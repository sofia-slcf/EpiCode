function proba_density(config,rongeur,per_rat)

S = {'CC','HPC', 'NC', 'PTA', 'S1', 'TH'};

if ispc
    savedir = '\\lexport\iss01.charpier\analyses\wod\Sofia\images';
elseif isunix
    savedir = '/network/lustre/iss01/charpier/analyses/wod/Sofia/images';
end

%% pour chaque structure
for icat=1:6
    wod_ini=[];
    for irat=1:length(config)
        if irat==4
            continue
        end
        
        for itrial=1%:size(rongeur(irat).catego(icat).wod_time_ini,2)
            
            wod_ini(irat)=rongeur(irat).catego(icat).wod_time_ini(itrial);
            
        end
        wod_ini(wod_ini==0) = nan ;
        
    end
    
    
    
    %% plot histogram
    fig1=figure(1);
    hist = histogram(wod_ini,100:20:250); %histogram(wod_ini,6) binage pas terrible
    legend(S);
    xlabel ('time (s)','Fontsize',15);
    ylabel ('anaimal proportion','Fontsize',15);
    title({'Distribution of Intitating WoD accross brain structure (6 bins)'} ,'Fontsize',12)
    hold on
    
    %% plot probability density % pas l'idéal : trop peu precis
    fig2=figure(2);
    [f,xi] = ksdensity(wod_ini); 
    plot(xi,f,'Linewidth',1.5)
    legend(S);
    xlabel ('time (s)','Fontsize',15);
    ylabel ('Density','Fontsize',15);
    title({'Density Probabilty of Intitating WoD accross brain structure '} ,'Fontsize',12)
    hold on
    
    %% plot histogram as a line, smoothed % beaucoup plus precis
    fig3=figure(3);
    x = hist.BinEdges(1:end-1) + diff(hist.BinEdges)/2; %center of each histogram bar
    xinterp = linspace(x(1), x(end), 1000);
    yinterp = pchip(x, hist.Values, xinterp);
    plot(xinterp, yinterp, 'linewidth', 2);
    legend(S);
    xlabel ('time (s)','Fontsize',15);
    ylabel ('anaimal proportion','Fontsize',15);
    title({'Smoothed distribution of Intitating WoD accross brain structure (6 bins)'} ,'Fontsize',12)
    set(gca, 'tickdir', 'out', 'fontsize', 15);
    hold on
    
    %% plot histogram as a line % raw representation
    fig4=figure(4);
    x = hist.BinEdges(1:end-1) + diff(hist.BinEdges)/2; %center of each histogram bar
    plot(x, hist.Values, 'linewidth', 2);
    legend(S);
    xlabel ('time (s)','Fontsize',15);
    ylabel ('anaimal proportion','Fontsize',15);
    title({'Raw distribution of Intitating WoD accross brain structure '} ,'Fontsize',12)
    set(gca, 'tickdir', 'out', 'fontsize', 15);
    hold on
    
end
hold off;hold off;hold off;hold off;

if ~exist(fullfile(savedir,'density'))
mkdir(fullfile(savedir,'density'));
end
saveas(fig1,fullfile(savedir,'density', ['histogram_initiating.jpg']));
saveas(fig2,fullfile(savedir,'density', ['KS_Density_initiating.jpg']));
saveas(fig3,fullfile(savedir,'density', ['Density_initiating.jpg']));
saveas(fig4,fullfile(savedir,'density', ['Raw_initiating.jpg'])); 
close(fig1); close(fig2);close(fig3); close(fig4);

%%
[per_rat_normwodmin,per_rat_normdelta]=wod_time_norm(config,rongeur,per_rat);

for icat=1:6
    wod_ini=[];
    for irat=1:length(config)
        if irat==4
            continue
        end
        
        for itrial=1%:size(rongeur(irat).catego(icat).wod_time_ini,2)
            
            wod_ini(irat)=per_rat_normwodmin(irat).trial(itrial).struct(icat).wod_time(itrial);
            
        end
        wod_ini(wod_ini==0) = nan ;
        
    end
  
    %% plot histogram
    fig5=figure(5);
    hist = histogram(wod_ini,0:8:100); %histogram(wod_ini,6) binage pas terrible
    legend(S);
    xlabel ('time (s)','Fontsize',15);
    ylabel ('animal proportion','Fontsize',15);
    title({'Distribution of normalized intitating WoD accross brain structure '} ,'Fontsize',12)
    hold on
    
    %% plot probability density % pas l'idéal : trop peu precis
    fig6=figure(6);
    [f,xi] = ksdensity(wod_ini); 
    plot(xi,f,'Linewidth',1.5)
    legend(S);
    xlabel ('time (s)','Fontsize',15);
    ylabel ('Density','Fontsize',15);
    title({'Density Probabilty of normalized intitating WoD accross brain structure '} ,'Fontsize',12)
    hold on
    
    %% plot histogram as a line, smoothed % beaucoup plus precis
    fig7=figure(7);
    x = hist.BinEdges(1:end-1) + diff(hist.BinEdges)/2; %center of each histogram bar
    xinterp = linspace(x(1), x(end), 1000);
    yinterp = pchip(x, hist.Values, xinterp);
    plot(xinterp, yinterp, 'linewidth', 2);
    legend(S);
    xlabel ('time (s)','Fontsize',15);
    ylabel ('animal proportion','Fontsize',15);
    title({'Smoothed distribution of Normalized intitating WoD accross brain structure '} ,'Fontsize',12)
    set(gca, 'tickdir', 'out', 'fontsize', 15);
    hold on
    
    %% plot histogram as a line % raw representation
    fig8=figure(8);
    x = hist.BinEdges(1:end-1) + diff(hist.BinEdges)/2; %center of each histogram bar
    plot(x, hist.Values, 'linewidth', 2);
    legend(S);
    xlabel ('time (s)','Fontsize',15);
    ylabel ('animal proportion','Fontsize',15);
    title({'Raw distribution of normalized Intitating WoD accross brain structure '} ,'Fontsize',12)
    set(gca, 'tickdir', 'out', 'fontsize', 15);
    hold on
    
end
hold off;hold off;hold off;hold off;

if ~exist(fullfile(savedir,'density'))
mkdir(fullfile(savedir,'density'));
end
saveas(fig5,fullfile(savedir,'density', ['histogram_initiating_norm.jpg']));
saveas(fig6,fullfile(savedir,'density', ['KS_Density_initiating_norm.jpg']));
saveas(fig7,fullfile(savedir,'density', ['Density_initiating_norm.jpg']));
saveas(fig8,fullfile(savedir,'density', ['Raw_initiating_norm.jpg'])); 
close(fig5); close(fig6);close(fig7); close(fig8);




end