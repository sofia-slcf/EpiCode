function freq_data_stats(freq_data,cfg)


%% ANOVA for all comparison factors

%make table with all comparison variables
depth_start = 10;
depth_end = 2210; %�m
depth_step = 100;

statstable = table.empty;

analysis_names={'timefreq_wod','timefreq_wod_timenorm','timefreq_wod_blcorrected','timefreq_wod_timenorm_blcorrected'};
freq_names=["HF","LF"];


for inorm=1:size(analysis_names,2)
    
    %% Make table with all data
    irow = 0;
    
    
    for irat = 1:size(freq_data,2)
        
        if isempty(cfg{irat})
            continue
        end
        
        for freqband = ["HF" "LF"]
            
            if freqband=="HF"
                freq_idx=1;
            elseif freqband=="LF"
                freq_idx=2;
            end
            
            
            for idata= ["peak_time" "peak_value"]
                for itrial = 1:size(freq_data{irat}.(analysis_names{inorm}).(idata).(freqband),2)
                    %for ichan = 1:size(cfg{irat}.LFP.chan_depth,2) %channame
                    for idepth = depth_start:depth_step:depth_end
                        sel = abs(freq_data{irat}.Depth(:, itrial) - idepth) < depth_step/2;
                        if sum(sel) == 1
                            irow = irow + 1;
                            statstable.peaktime(irow) = freq_data{irat}.(analysis_names{inorm}).peak_time.(freqband)(sel,itrial);
                            statstable.peakval(irow) = freq_data{irat}.(analysis_names{inorm}).peak_value.(freqband)(sel,itrial);
                            statstable.chandepth(irow) = idepth;%cfg{irat}.LFP.chan_depth{ichan};
                            statstable.freq(irow) = freq_idx;
                            statstable.trial(irow) = itrial;
                            statstable.rat(irow) = irat;
                        elseif sum(sel) == 0
                            continue
                            
                        elseif sum(sel) > 1
                            error('it should have only one electrode for one deepness');
                        end
                        
                        
                    end %ichan
                end %itrial
            end %idata
        end %freqband
    end %irat
    
    %delete NaNs
    A=find(isnan(statstable.peaktime));
    statstable([A],:)=[];
    
    
    
    %% Make ANOVA for time and value of peak
    
    anovasavedir= fullfile(cfg{5}.statsavedir,'freq_data','anova');
    
    if ~isfolder(anovasavedir)
        mkdir(anovasavedir);
    end
    
    fname_peaktime=fullfile(anovasavedir,sprintf('%s_peaktime',analysis_names{inorm}));
    fname_peakval=fullfile(anovasavedir,sprintf('%s_peakval',analysis_names{inorm}));
    
    
    mdl_peaktime              = fitlm(statstable, 'peaktime ~ chandepth + freq + trial + chandepth:freq + chandepth:trial + freq:trial'); %group:time (equivalent à time:group) : effet du groupe auquel l'effet du temps a été retiré
    stats_peaktime            = anova(mdl_peaktime,'component');
    writetable(stats_peaktime,fname_peaktime,'WriteRowNames',true,'FileType','spreadsheet');
    
    mdl_peakval               = fitlm(statstable, 'peakval ~ chandepth  + freq + trial + chandepth:freq + chandepth:trial + freq:trial'); %group:time (equivalent à time:group) : effet du groupe auquel l'effet du temps a été retiré
    stats_peakval             = anova(mdl_peakval,'component');
    writetable(stats_peakval,fname_peakval,'WriteRowNames',true,'FileType','spreadsheet');
    
    
    %% Compare values 2 by 2
    
    pval_depthtable=table.empty;
    %Compare peak values according to depth for each frequency band
    irow=0;
    for idepth = depth_start:depth_step:depth_end
        irow=irow+1;
        for ifreq= 1:4
            
            if idepth==depth_end
                continue
            end
            
            sel= statstable.chandepth==idepth & statstable.freq==ifreq;
            datachan1= statstable(sel,:);
            
            
            sel_2=statstable.chandepth==idepth+depth_step & statstable.freq==ifreq;
            datachan2= statstable(sel_2,:);
            
            if isempty(datachan1) | isempty(datachan2)
                continue
            end
            
            %compare peakval between depth
            p(irow,ifreq)= ranksum(datachan1.peakval,datachan2.peakval);
            
            %Make table with p-values
            pval_depthtable.depth(irow)=idepth+depth_step;
            pval_depthtable.(freq_names{ifreq})(irow)= p(irow,ifreq);
        end %ifreq
    end %idepth
    
    %Correcting p-values
    [~, ~, ~, adj_p]=fdr_bh(p);
    adj_pval_peakval.(analysis_names{inorm})=adj_p;
    clear p adj_p
    
    %Make table with corrected p-values
    adj_pval_depthtable=table.empty;
    adj_pval_depthtable.depth=pval_depthtable.depth;
    for ifreq=1:2
        adj_pval_depthtable.(freq_names{ifreq})=adj_pval_peakval.(analysis_names{inorm})(:,ifreq);
    end %ifreq
    
    %Save tables
    fname_peakvaldepth= fullfile(anovasavedir,'multiple_comp',sprintf('peak_values_%s',analysis_names{inorm}));
    fname_adjpeakvaldepth= fullfile(anovasavedir,'multiple_comp',sprintf('adjusted_peak_values_%s',analysis_names{inorm}));
    
    if ~isfolder(fullfile(anovasavedir,'multiple_comp'))
        mkdir(fullfile(anovasavedir,'multiple_comp'));
    end
    writetable(pval_depthtable,fname_peakvaldepth,'FileType','spreadsheet');
    writetable(adj_pval_depthtable,fname_adjpeakvaldepth,'FileType','spreadsheet');
    
    
    
    %Compare peak timings according to frequency bands
    pval_freqtable=table.empty;
    irow=0;
    for ifreq= 1:4
        irow=irow+1;
        
        if ifreq==4
            continue
        end
        
        sel= statstable.freq==ifreq;
        datafreq1= statstable.peaktime(sel);
        sel_2=statstable.freq==ifreq+1;
        datafreq2= statstable.peaktime(sel_2);
        
        if isempty(datafreq1) | isempty(datafreq2)
            continue
        end
        %store p-values
        p(irow,ifreq)= ranksum(datafreq1,datafreq2);
        
        %Make table
        pval_freqtable.(freq_names{ifreq})(irow)= p(irow,ifreq);
    end %ifreq
    
    [~, ~, ~, adj_p]=fdr_bh(p);
    adj_pval_peaktime.(analysis_names{inorm})=adj_p;
    clear p adj_p
    
    adj_pvalpeaktime_table=table.empty;
    
    for ifreq=1:2
        adj_pvalpeaktime_table.(freq_names{ifreq})=adj_pval_peaktime.(analysis_names{inorm});
    end
    
    fname_pvaltime= fullfile(anovasavedir,'multiple_comp',sprintf('peak_timings_freq_%s',analysis_names{inorm}));
    fname_adjpvaltime= fullfile(anovasavedir,'multiple_comp',sprintf('adj_peak_timings_freq_%s',analysis_names{inorm}));
    
    writetable(pval_freqtable,fname_pvaltime,'FileType','spreadsheet');
    writetable(adj_pvalpeaktime_table,fname_adjpvaltime,'FileType','spreadsheet');
    
    % compare peak timings according to trials
    pval_trialtable=table.empty;
    itrial=1;
    for ifreq=1:4
        
        sel= statstable.trial== itrial & statstable.freq==ifreq;
        datatrial1= statstable.peaktime(sel);
        sel_2=statstable.trial== itrial+1 & statstable.freq==ifreq;
        datatrial2= statstable.peaktime(sel_2);
        
        if isempty(datafreq1) | isempty(datafreq2)
            continue
        end
        %store p-values
        p(ifreq)= ranksum(datafreq1,datafreq2);
        
        %Make table
        pval_trialtable.(freq_names{ifreq})(irow)= p(ifreq);
    end %ifreq
    
    fname_pvaltrial=fullfile(anovasavedir,'multiple_comp',sprintf('peak_timings_%s',analysis_names{inorm}));
    writetable(pval_trialtable,fname_pvaltrial,'FileType','spreadsheet');
    
    
    %% Make boxplots of comparisons
    
    boxplotpath=fullfile(cfg{5}.imagesavedir,'freq_band_peaks','boxplots');
    
    if ~isfolder(boxplotpath)
        mkdir(boxplotpath);
    end
    
    %plot peak timings according to frequency band
    boxplot(statstable.peaktime,statstable.freq,'Colors','k','Symbol','o','Labels',freq_names);
    
    freqTimes=gcf;
    fname_freqtimes=fullfile(boxplotpath,sprintf('peaktime_frequency_%s',analysis_names{inorm}));
    dtx_savefigure(freqTimes,fname_freqtimes,'pdf','png','close');
    
    %plot peak timings per frequency band and trial
    
    %preparing data
    sel_ftrial=statstable.trial==1;
    peaktime_ftrial=statstable(sel_ftrial,:);
    
    sel_strial=statstable.trial==2;
    peaktime_strial=statstable(sel_strial,:);
    
    for ifreq=1:2
        sel=peaktime_ftrial.freq==ifreq;
        sel_2=peaktime_strial.freq==ifreq;
        %first trial
        ftrialpeaktime_perfreq(:,ifreq)=peaktime_ftrial.peaktime(sel);
        %second trial
        strialpeaktime_perfreq(:,ifreq)=peaktime_strial.peaktime(sel_2);
        
    end
    
    % put data into cell arrays
    data=cell(2,2);
    for ii=1:size(data,1)
        ftrialpeaktime_perfreq_c{ii}=ftrialpeaktime_perfreq(:,ii);
        strialpeaktime_perfreq_c{ii}=strialpeaktime_perfreq(:,ii);
    end
    data=vertcat(ftrialpeaktime_perfreq_c,strialpeaktime_perfreq_c);
    
    xlab={'HF','LF'};
    %determine colors
    col=[255,64,64, 200;
        255,127,36, 200];
    col=col/255;
    
    multiple_boxplot(data',xlab,{'first', 'second'},col');
    title('Power peak by frequency and trials')
    ylabel('Time from Vent. Off (s)','FontWeight','bold')
    fig_trials=gcf;
    
    fname_figtrials=fullfile(boxplotpath,sprintf('peaktime_freqtrials_%s',analysis_names{inorm}));
    dtx_savefigure(fig_trials,fname_figtrials,'pdf','png','close');
    
    
end %inorm



end %function