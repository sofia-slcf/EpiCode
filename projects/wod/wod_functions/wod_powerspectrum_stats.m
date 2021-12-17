function wod_powerspectrum_stats(config, power_detection, analysis_list)

%power_detection : 
% - pour chaque structure, temps où pic max est détecté, par bande de
% fréquence





analysis_names = {'timefreq_wod', 'timefreq_wod_timenorm', 'timefreq_baseline','timefreq_wod_blcorrected', 'timefreq_wod_timenorm_blcorrected', 'timefreq_baseline_blcorrected','log_timefreq_wod', 'log_timefreq_wod_timenorm', 'log_timefreq_baseline','log_timefreq_wod_blcorrected', 'log_timefreq_wod_timenorm_blcorrected','log_timefreq_baseline_blcorrected'};
brain_struct = {'CC','HPC', 'NC', 'PTA', 'S1', 'TH'};
brain_struct_color = linspecer(length(brain_struct));

% if ispc
%     savedir = '\\lexport\iss01.charpier\analyses\wod\Sofia\images';
% elseif isunix
%     savedir = '/network/lustre/iss01/charpier/analyses/wod/Sofia/images';
% end

do_plot=true; %false

for analysis = analysis_list
    
%% organiser les donner pour un plot statistic
    for icat=1:6%size(config.classe.structure,2)
        nLF(icat)=length(catego2stat(icat).peak_LF);
        nTF(icat)=length(catego2stat(icat).peak_TF);
        nHF(icat)=length(catego2stat(icat).peak_HF);
        nVHF(icat)=length(catego2stat(icat).peak_VHF)
    end
    N=max([max(nLF),max(nTF),max(nHF),max(nVHF)]);
    data=NaN(N,6) ;
    jcat=1;
    for icat=1:6%size(config.classe.structure,2)
        data(1:length(catego2stat(icat).peak_LF),jcat)=catego2stat(icat).peak_LF';
        data(1:length(catego2stat(icat).peak_TF),jcat+1)=catego2stat(icat).peak_TF';
        data(1:length(catego2stat(icat).peak_HF),jcat+2)=catego2stat(icat).peak_HF';
        data(1:length(catego2stat(icat).peak_VHF),jcat+3)=catego2stat(icat).peak_VHF';
        jcat=jcat+4;
    end
    
    
    
    %% plot les peak de range du power spectrum across brain regions
    
    sd=nanstd(data);
    m=nanmean(data);
    
    fig5=figure(5)
    x=data;
    n=length(data); xx=([1:24])';
    r=repmat(xx,1,n)';
    g=r(:)';
    positions = 1:24;%[1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24];
    h=boxplot(x,g, 'positions', positions);
    set(h,'linewidth',2)
    set(gca,'xtick',[mean(positions(1:4)) mean(positions(5:8)) mean(positions(9:12)) mean(positions(13:16)) mean(positions(17:20)) mean(positions(21:24))])
    set(gca,'xticklabel',brain_struct,'Fontsize',20)
    color = repmat({'c'; 'y' ;'r';'b'} ,6,1)'
    h = findobj(gca,'Tag','Box');
    for j=1:length(h)
        patch(get(h(j),'XData'),get(h(j),'YData'),cell2mat(color(j)),'FaceAlpha',0.5);
    end
    yt = get(gca, 'YTick');
    axis([xlim    0  ceil(max(yt)*1.2)])
    xt=positions
    text(xt(1)-0.3,m(1)+3500,{[num2str(m(1))],['\pm', num2str(sd(1))]})
    for p=2:length(xt)
        text(xt(p)-0.3,m(p)+2000,{[num2str(m(p))],['\pm', num2str(sd(p))]})
    end
    hold on
    lgd=legend({'LF','TF','HF','VHF'},'Fontsize',16)
    xlabel ('Structures')
    ylabel ('Peak delay from Vent-off','Fontsize',15)
    title({['Comparaison accross brain regions of the '], ['Power spectrum peak delay from Vent Off of 3 frenquency ranges N=22']},'Fontsize',15)
    figname =fullfile(config{1}.imagesavedir,'Statistics',['peak_power_spectrum']);
    dtx_savefigure(fig5, figname, 'png', 'pdf', 'jpeg', 'close');
    hold off
    
    %% statistics ANOVAN des variable frequce /structure
    
    groupeHF = repmat({'HF'} ,size(data,1),1)'
    groupeLF = repmat({'LF'} ,size(data,1),1)'
    groupeTF = repmat({'TF'} ,size(data,1),1)'
    groupeTF = repmat({'VHF'} ,size(data,1),1)'
   
    groupe_cat= repmat(brain_struct ,size(data,1)*4,1)
    groupe= repmat([groupeHF; groupeTF; groupeLF; groupeVHF] ,size(brain_struct,2),1)'
    
    Groupe_cat= reshape(groupe_cat,[],1)';
    Groupe=reshape(groupe,[],1)';
    Y=reshape(data,1,[]);
    
    [p,tbl,stats,terms] = anovan(Y,{Groupe,Groupe_cat},'model','interaction',...
        'varnames',{'Groupe','Groupe_cat'},'display','on')
    results = multcompare(stats,'Dimension',[1 2])
    
    %% Anova 1 facteur comparaison d'une frequence entre les structures
    groupe_cat= repmat(brain_struct ,size(data,1),1);
    Groupe=reshape(groupe_cat,[],1)';
    dataHF_stat =[];
    dataTF_stat =[];
    dataLF_stat =[];
    dataVHF_stat =[];
    i=1;j=2;k=3;
    while i<=size(data,2)
        dataHF_stat  = [dataHF_stat data(:,i)'];
        dataTF_stat  = [dataTF_stat data(:,j)'];
        dataLF_stat  = [dataLF_stat data(:,k)'];
        dataVHF_stat = [dataVHF_stat data(:,m)'];
        i=i+3;
        j=j+3;
        k=k+3;
        m=m+4;
    end
    
    [p_HF,tbl_HF,statsHF,terms] = anovan(dataHF_stat,{Groupe},'model','interaction',...
        'varnames',{'Groupe'},'display','on');
    resultsHF = multcompare(statsHF,'Dimension',[1]);
    [p_LF,tbl_LF,statsLF,terms] = anovan(dataLF_stat,{Groupe},'model','interaction',...
        'varnames',{'Groupe'},'display','on');
    resultsLF = multcompare(statsLF,'Dimension',[1]);
    [p_TF,tbl_TF,statsTF,terms] = anovan(dataTF_stat,{Groupe},'model','interaction',...
        'varnames',{'Groupe'},'display','on');
    resultsTF = multcompare(statsTF,'Dimension',[1]);
    [p_VHF,tbl_VHF,statsVHF,terms] = anovan(dataVHF_stat,{Groupe},'model','interaction',...
        'varnames',{'Groupe'},'display','on');
    resultsVHF = multcompare(statsVHF,'Dimension',[1]);
    
    
end
    