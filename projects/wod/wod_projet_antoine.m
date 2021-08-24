function wod_project_antoine(slurm_task_id)

% Ce script projet sert Ã  calculer les donnÃ©es de chaque rat
% 'irat' est en input car il prendra la valeur du array slurm sur le cluster
% Rassembler les donnÃ©es de tous les rats pour faire des moyennes ou stats
% sur tous les rats : dans un autre script (si possible organisÃ© comme
% celui-ci)

% a faire pour rÃ©organiser les fonctions WOD :
% - mettre en input : cfg, MuseStruct, LFP, et toutes les structures qui sont calculÃ©es
% dans une autre fonction
% - retirer les addpath, les config = wod_setparams, les boucles irat
% - remplacer tous les config{irat} par cfg
% - si possible, sauvegarder les donnÃ©es en fin de script, et crÃ©er la
% possibilitÃ© de les charger avec l'argument 'force'

%% set parameters
try %en local
    scriptpath = matlab.desktop.editor.getActiveFilename;
catch %cluster
    scriptpath = mfilename('fullpath');
end

epicodepath = [fileparts(fileparts(fileparts(scriptpath))), filesep];

addpath (genpath([epicodepath,'development']))
addpath (genpath([epicodepath,'shared']))
addpath (genpath([epicodepath,'external']))
addpath (genpath([epicodepath,'templates']))
addpath (genpath([epicodepath,'projects', filesep, 'wod']))
addpath (genpath([epicodepath,'projects', filesep, 'dtx']))
addpath (genpath([epicodepath,'projects', filesep, 'wod',filesep,'wod_functions']))

if ispc
    addpath \\lexport\iss01.charpier\analyses\wod\fieldtrip-20200607
elseif isunix
    addpath /network/lustre/iss01/charpier/analyses/wod/fieldtrip-20200607
end

ft_defaults
%remove fieldtrip's output
ft_warning off
ft_notice off
ft_info off
ft_debug off
global ft_default
ft_default.checkconfig = 'silent';
ft_default.checkpath = 'once';
ft_default.showcallinfo = 'no';
ft_default.trackcallinfo = 'no';
ft_default.tracktimeinfo = 'no';

config16 = wod_setparams;
config32 = wod_setparams_32chan;
config = config16;
cfgorig = config;

ipart= 1;

anovastatpath=fullfile(config16{4}.statsavedir,'Waves_detection','anova');
tablesavedir=fullfile(fileparts(config16{4}.datasavedir),'tables','matlab');


if ~isfolder(anovastatpath)
    mkdir(anovastatpath);
end

if ~isfolder(tablesavedir)
    mkdir(tablesavedir)
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if slurm_task_id >0
    %% analysis by rat
    
    for irat= slurm_task_id
        
        if isempty(config{irat})
            continue
        end
        
        %find concatenated LFP (see wod_concatenateLFP.m)
        [~,dir_name]                       = fileparts(cfgorig{irat}.rawdir);
        config{irat}.rawdir                = fullfile(config{irat}.concatdata_path);
        %when config 32
        %config{irat}.directorylist{ipart}  = {config{irat}.prefix};
        
        %when config16
        config{irat}.directorylist{ipart}(1,:)=[];
        config{irat}.directorylist{ipart}{ipart}=dir_name;
        %read Muse markers
        MuseStruct = readMuseMarkers(config{irat}, true);
        %save(fullfile(config{irat}.datasavedir,sprintf('%s-MuseStruct.mat',config{irat}.prefix)),'MuseStruct');
        
        %read LFP, append electrodes, and cut into trials according to Muse Markers
        LFP = readLFP(config{irat}, MuseStruct, true);
        LFP = LFP{1}.(config{irat}.LFP.name{1}); %remove this 'epicode' organisation for now.
        %end
        
        %rename channels according to depth
        for ichan = 1:size(config{irat}.LFP.channel, 2)
            idx = strcmp(config{irat}.LFP.channel{ichan}, LFP.label);
            label_renamed{idx} = config{irat}.LFP.rename{ichan};
        end
        LFP.label = label_renamed';
        clear label_renamed
        
        %vÃ©rifier qu'il y a bien autant de trials que de marqueurs Vent_Off
        startmarker = config{irat}.muse.startmarker.(config{irat}.LFP.name{1});
        if size(LFP.trial,2) ~= size(MuseStruct{1}{1}.markers.(startmarker).synctime,2)
            error('Not the same number of trials that of marker start for %s. \nCheck that begin/end of each trial is not before start of file or after end of file', config{irat}.prefix(1:end-1));
        end
        
        
        % Compute TFR for each rat
        %wod_tfr_compute(config{irat}, MuseStruct,LFP);
        
        %Plot TFR data for each rat
        %wod_tfr_plotrat(config{irat});
        
        
        
        
    end %irat
end %if slurm_task_id

if slurm_task_id==0
    
    %% Analysis for pooled rats
    %% Waves delay data and separated
    %Detect waves, extract timings, origin depth and propagation speed
    stats_all=wod_wavedetection([config16 config32],false);
    
    %make stats on calculated data
    wod_calculated_data_stats([config16,config32],stats_all);
    %% Pool data according to depth
    config=[config16 config32];
    
    %order protocols starting at 1000 µm according to channel depth
    depth_start = 10;
    depth_end = 2210; %µm
    depth_step = 100;
    
    icol = 0;
    for idepth = depth_start:depth_step:depth_end %step de 250 car les électrodes sont espacées de 250µm
        icol = icol+1;
        irow = 0;
        for irat = 1:size(stats_all, 2)
            if isempty(config{irat})
                continue
            end
            
            for itrial = 1:size(stats_all{irat}.WoD.peak_time, 2)
                if stats_all{irat}.oridepthclass(itrial)==1000
                    irow = irow+1;
                    sel = abs(stats_all{irat}.Depth(:, itrial) - idepth) < depth_step/2;
                    for iwave=["WoD" "WoR"]
                        
                        allfields=fieldnames(stats_all{irat}.(iwave))';
                        excludefields = allfields(~ismember(allfields, {'speed_up', 'speed_dn'}));  %remove fields 'E' and 'G'
                        
                        for ifield = string(excludefields)
                            
                            if sum(sel) == 1
                                ordered_data_1000.(iwave).(ifield)(irow, icol) = stats_all{irat}.(iwave).(ifield)(sel, itrial);
                                ordered_data_1000.Depth(irow,icol)= stats_all{irat}.Depth(sel,itrial);
                                ordered_data_1000.plateau_duration(irow,icol)= stats_all{irat}.plateau_duration(sel,itrial);
                            elseif sum(sel) == 0
                                ordered_data_1000.(iwave).(ifield)(irow, icol) = nan;
                                ordered_data_1000.Depth(irow,icol)= nan;
                                ordered_data_1000.plateau_duration(irow,icol)=nan;
                            elseif sum(sel) > 1
                                error('it should have only one electrode for one deepness');
                            end
                            
                        end %ifield
                    end %iwave
                elseif stats_all{irat}.oridepthclass(itrial)==1400
                    irow = irow+1;
                    sel = abs(stats_all{irat}.Depth(:, itrial) - idepth) < depth_step/2;
                    for iwave=["WoD" "WoR"]
                        
                        allfields=fieldnames(stats_all{irat}.(iwave))';
                        excludefields = allfields(~ismember(allfields, {'speed_up', 'speed_dn'}));  %remove fields 'E' and 'G'
                        
                        for ifield = string(excludefields)
                            
                            if sum(sel) == 1
                                ordered_data_1400.(iwave).(ifield)(irow, icol) = stats_all{irat}.(iwave).(ifield)(sel, itrial);
                                ordered_data_1400.Depth(irow,icol)= stats_all{irat}.Depth(sel,itrial);
                                ordered_data_1400.plateau_duration(irow,icol)= stats_all{irat}.plateau_duration(sel,itrial);
                            elseif sum(sel) == 0
                                ordered_data_1400.(iwave).(ifield)(irow, icol) = nan;
                                ordered_data_1400.Depth(irow,icol)= nan;
                                ordered_data_1400.plateau_duration(irow,icol)=nan;
                                
                            elseif sum(sel) > 1
                                error('it should have only one electrode for one deepness');
                            end
                            
                        end %ifield
                    end %iwave
                end
            end %itrial
        end
    end
    
    %remove empty rows from arrays
    for iwave= ["WoD" "WoR"]
        for ifield= string(fieldnames(ordered_data_1000.(iwave))')
            toremove_1=ordered_data_1000.(iwave).(ifield)==0;
            toremove_1=all(toremove_1, 2);
            ordered_data_1000.(iwave).(ifield)=ordered_data_1000.(iwave).(ifield)(~toremove_1, :);
            toremove_2=ordered_data_1400.(iwave).(ifield)==0;
            toremove_2=all(toremove_2, 2);
            ordered_data_1400.(iwave).(ifield)=ordered_data_1400.(iwave).(ifield)(~toremove_2, :);
        end %ifield
    end %iwave
    ordered_data_1000.Depth=ordered_data_1000.Depth(~toremove_1, :);
    ordered_data_1400.Depth=ordered_data_1400.Depth(~toremove_2, :);
    ordered_data_1000.plateau_duration=ordered_data_1000.plateau_duration(~toremove_1,:);
    ordered_data_1400.plateau_duration=ordered_data_1400.plateau_duration(~toremove_2,:);
    
    %plot delays pooled 1000µm origin
    wod_plot_delays(ordered_data_1000,config);
    %% Plot plateau duration according to depth
    
    %For superficial origins
    for irat = 1:size(ordered_data_1000.plateau_duration, 1)
        startinterp = find(~isnan(ordered_data_1000.plateau_duration(irat,:)), 1, 'first');
        endinterp = find(~isnan(ordered_data_1000.plateau_duration(irat,:)), 1, 'last');
        ordered_data_interp1000.plateau_duration(irat,startinterp:endinterp) = fillmissing(ordered_data_1000.plateau_duration(irat,startinterp:endinterp), 'linear'); %ou pchip, ou spline
        ordered_data_interp1000.Depth(irat,startinterp:endinterp)= fillmissing(ordered_data_1000.Depth(irat,startinterp:endinterp), 'linear');
    end
    
    %For deep origins
    for irat = 1:size(ordered_data_1400.plateau_duration, 1)
        startinterp = find(~isnan(ordered_data_1400.plateau_duration(irat,:)), 1, 'first');
        endinterp = find(~isnan(ordered_data_1400.plateau_duration(irat,:)), 1, 'last');
        ordered_data_interp1400.plateau_duration(irat,startinterp:endinterp) = fillmissing(ordered_data_1400.plateau_duration(irat,startinterp:endinterp), 'linear'); %ou pchip, ou spline
        ordered_data_interp1400.Depth(irat,startinterp:endinterp)= fillmissing(ordered_data_1400.Depth(irat,startinterp:endinterp), 'linear');
    end
    
    %replace zeros by NaNs
    idx_zeros_1000= find(ordered_data_interp1000.plateau_duration(:,:)==0);
    idx_zeros_1400= find(ordered_data_interp1400.plateau_duration(:,:)==0);
    
    idx_zeros_depth_1000=find(ordered_data_interp1000.Depth(:,:)==0);
    idx_zeros_depth_1400=find(ordered_data_interp1400.Depth(:,:)==0);
    
    ordered_data_interp1400.plateau_duration(idx_zeros_1400)=nan;
    ordered_data_interp1000.plateau_duration(idx_zeros_1000)=nan;
    ordered_data_interp1400.Depth(idx_zeros_depth_1400)=nan;
    ordered_data_interp1000.Depth(idx_zeros_depth_1000)=nan;
    
    clear idx_zeros_1000 idx_zeros_1400 idx_zeros_depth_1000 idx_zeros_depth_1400
    
    %calculate median and standard deviation for superficial origin
    for ichan= 1:22%size(ordered_data.(iwave).(itime),2)
        data_med_1000(1,ichan)= nanmedian(ordered_data_interp1000.plateau_duration(:,ichan));
        data_std_1000(1,ichan)= std(ordered_data_interp1000.plateau_duration(:,ichan),'omitnan');
        data_med_depth_1000(1,ichan)=nanmedian(ordered_data_interp1000.Depth(:,ichan));
        data_mean_1000(1,ichan)=nanmean(ordered_data_interp1000.plateau_duration(:,ichan));
        data_sem_1000=data_std_1000 / (sqrt(length(ordered_data_interp1000.plateau_duration)));
        data_mean_depth_1000(1,ichan)=nanmean(ordered_data_interp1000.Depth(:,ichan));
    end %ichan
    
    %calculate median and standard deviation for deep origin
    for ichan= 1:22%size(ordered_data.(iwave).(itime),2)
        data_med_1400(1,ichan)= nanmedian(ordered_data_interp1400.plateau_duration(:,ichan));
        data_std_1400(1,ichan)= std(ordered_data_interp1400.plateau_duration(:,ichan),'omitnan');
        data_med_depth_1400(1,ichan)=nanmedian(ordered_data_interp1400.Depth(:,ichan));
        data_mean_1400(1,ichan)=nanmean(ordered_data_interp1400.plateau_duration(:,ichan));
        data_sem_1400=data_std_1400 / (sqrt(length(ordered_data_interp1400.plateau_duration)));
        data_mean_depth_1400(1,ichan)=nanmean(ordered_data_interp1400.Depth(:,ichan));
    end %ichan
    clear data_mean_1000 data_mean_1400 data_mean_depth_1000 data_mean_depth_1400 data_med_1000 data_med_1400 data_med_depth_1000 data_med_depth_1400
    clear data_sem_1000 data_sem_1400 data_std_1000 data_std_1400
    %% Plot mean +/- SEM filled
    
    fig_meanraw=figure;hold
    
    C_prot={[0.5 0.5 0.5],[0.85 0.325 0.098]};
    C_med= {'k',[0.635 0.078 0.184]};
    C_std= {[0.5 0.5 0.5]; [1  0.549 0.412]};
    
    plot(data_mean_1000,data_mean_depth_1000,'Color',C_med{1},'LineWidth',1.5);
    patch([data_mean_1000- data_sem_1000, data_mean_1000(end:-1:1)+ data_sem_1000(end:-1:1)], [data_mean_depth_1000, data_mean_depth_1000(end:-1:1)], C_std{1}, 'facealpha', 0.3, 'edgecolor', 'none');
    
    plot(data_mean_1400,data_mean_depth_1400,'Color',C_med{2},'LineWidth',1.5);
    patch([data_mean_1400- data_sem_1400, data_mean_1400(end:-1:1)+ data_sem_1400(end:-1:1)], [data_mean_depth_1400, data_mean_depth_1400(end:-1:1)],  C_std{2}, 'facealpha', 0.3, 'edgecolor', 'none');
    
    
    ylim([0 2200]);
    xlim([0 80]);
    
    fname_meanraw=fullfile(config{4}.imagesavedir,'delays','plateau_mean_sem');
    
    if ~isfolder(fullfile(config{4}.imagesavedir,'delays'))
        mkdir(fullfile(config{4}.imagesavedir,'delays'));
    end
    %save figure
    dtx_savefigure(fig_meanraw,fname_meanraw,'png','pdf','close');
    clear C_prot C_med C_std
    %% Frequency band data
    
    %extract peak time and value for power bands
    freq_data=wod_tfr_extractdata([config16 config32],false);
    
    %gather 16 and 32 chans
    ordered_freqdata = wod_fusion_freqdata(freq_data,[config16 config32],stats_all,false);
    
    %plot freq data
    wod_plot_freqdata(ordered_freqdata,[config16 config32]);
    
    %freq data stats
    freq_data_stats(freq_data,[config16 config32]);
    
    
    %Plot average TFR
    %wod_tfr_grandaverage(config);
    %% Create table with every delay
    
    %make table with all comparison variables
    config_all=[config16 config32];
    delaytable = table.empty;
    irow=0;
    for irat=1:size(config_all,2)
        if isempty(config_all{irat})
            continue
        end
        
        for itrial=1:size(stats_all{irat}.WoD.peak_time,2)
            irow= irow+1;
            sel=stats_all{irat}.Depth(:,itrial)==stats_all{irat}.wod_origin_depth(itrial);
            
            delaytable.rat(irow)=irat;
            delaytable.trial(irow)=itrial;
            delaytable.surge_delay(irow)=freq_data{irat}.timefreq_wod_blcorrected.peak_time.HF(sel,itrial);
            delaytable.low_delay(irow)=freq_data{irat}.timefreq_wod_blcorrected.peak_time.LF(sel,itrial);
            delaytable.Iso_delay(irow)=stats_all{irat}.ISO(sel,itrial);
            delaytable.Wod_delay(irow)=stats_all{irat}.WoD.peak_time(sel,itrial);
            delaytable.plateau(irow)=stats_all{irat}.plateau_duration(sel,itrial);
            delaytable.recov(irow)=config_all{irat}.LFP.recov{itrial};
            delaytable.depthclass(irow)=stats_all{irat}.oridepthclass(itrial);
            
        end %itrial
    end %irat
    
    fname_delay=fullfile(tablesavedir,'table_delays');
    writetable(delaytable,fname_delay,'FileType','Spreadsheet')
    %% Make ANOVA for delays
    
    %HF surge
    mdl_surge= fitlm(delaytable, 'surge_delay ~ recov + depthclass + trial+  depthclass:recov + depthclass:trial + trial:recov');
    stats_surge           = anova(mdl_surge,'component');
    fname_surge=fullfile(anovastatpath,'anova_surge_time');
    
    %LF surge
    mdl_LF= fitlm(delaytable, 'low_delay ~ recov + depthclass + trial+  depthclass:recov + depthclass:trial + trial:recov');
    stats_LF           = anova(mdl_LF,'component');
    fname_LF=fullfile(anovastatpath,'anova_LF_time');
    
    %ISO delay
    mdl_ISO= fitlm(delaytable, 'Iso_delay ~ recov + depthclass + trial+  depthclass:recov + depthclass:trial + trial:recov');
    stats_ISO           = anova(mdl_ISO,'component');
    fname_ISO=fullfile(anovastatpath,'anova_ISO_time');
    
    %WoD delay
    mdl_WOD= fitlm(delaytable, 'Wod_delay ~ recov + depthclass + trial+  depthclass:recov + depthclass:trial + trial:recov');
    stats_WOD           = anova(mdl_WOD,'component');
    fname_WOD=fullfile(anovastatpath,'anova_WOD_time');
    
    %Plateau duration
    mdl_plateau= fitlm(delaytable, 'plateau ~ recov + depthclass + trial+  depthclass:recov + depthclass:trial + trial:recov');
    stats_plateau           = anova(mdl_plateau,'component');
    fname_plateau=fullfile(anovastatpath,'anova_plateau_time');
    
    %save tables
    writetable(stats_surge,fname_surge,'WriteRowNames',true,'FileType','spreadsheet');
    writetable(stats_LF,fname_LF,'WriteRowNames',true,'FileType','spreadsheet');
    writetable(stats_ISO,fname_ISO,'WriteRowNames',true,'FileType','spreadsheet');
    writetable(stats_WOD,fname_WOD,'WriteRowNames',true,'FileType','spreadsheet');
    writetable(stats_plateau,fname_plateau,'WriteRowNames',true,'FileType','spreadsheet');
    %% Make 2 by 2 comparisons by depth origin
    
    %compare by origin depth
    pval_delaytable=table.empty;
    %Compare peak values according to depth for each frequency band
    sel= delaytable.depthclass==1000;
    datachan1= delaytable(sel,:);
    
    
    sel_2=delaytable.depthclass==1400;
    datachan2= delaytable(sel_2,:);
    
    irow=0;
    for iana=["surge_delay" "low_delay" "Iso_delay" "Wod_delay" "plateau"]
        irow=irow+1;
        %compare peakval between depth
        p(irow)= ranksum(datachan1.(iana),datachan2.(iana));
        
        pval_delaytable.(iana)=p(irow);
    end %iana
    
    %correct all p-values
    [~, ~, ~, adj_p]=fdr_bh(p);
    
    fname_pvaldelay=fullfile(config16{4}.statsavedir,'Waves_detection','pval_delays_depthclass');
    writetable(pval_delaytable,fname_pvaldelay,'FileType','Spreadsheet')
    
    clear pval_delaytable p adj_p
    %% Make 2 by 2 comparisons by trial
    %compare by trials
    pval_delaytable=table.empty;
    %Compare peak values according to depth for each frequency band
    sel= delaytable.trial==1;
    datachan1= delaytable(sel,:);
    
    
    sel_2=delaytable.trial==2;
    datachan2= delaytable(sel_2,:);
    
    irow=0;
    for iana=["surge_delay" "low_delay" "Iso_delay" "Wod_delay" "plateau"]
        irow=irow+1;
        %compare peakval between depth
        p(irow)= ranksum(datachan1.(iana),datachan2.(iana));
        
        pval_delaytable.(iana)=p(irow);
    end %iana
    
    %correct all p-values
    [~, ~, ~, adj_p]=fdr_bh(p);
    
    fname_pvaldelay=fullfile(config16{4}.statsavedir,'Waves_detection','pval_delays_trials');
    writetable(pval_delaytable,fname_pvaldelay,'FileType','Spreadsheet')
    
    clear pval_delaytable p adj_p
    %% Make scatter plot for all depths
    
    
    figure_depths=figure;
    % Plot 1000
    irow=0;
    for irat=1:size(stats_all,2)
        if isempty(stats_all{irat})
            continue
        end
        for itrial=1:size(stats_all{irat}.oridepthclass,2)
            if stats_all{irat}.oridepthclass(itrial)==1000
                irow=irow+1;
                Value_1000(irow) = stats_all{irat}.wod_origin_depth(itrial);
            end
        end
    end
    
    x = ones(1, length(Value_1000));
    scatter(x, Value_1000,70, 'ko',  'LineWidth', 1,'MarkerFaceColor','w','jitter', 'on', 'jitterAmount', 0.05);
    hold on;
    % Plot 14000
    irow=0;
    for irat=1:size(stats_all,2)
        if isempty(stats_all{irat})
            continue
        end
        for itrial=1:size(stats_all{irat}.oridepthclass,2)
            if stats_all{irat}.oridepthclass(itrial)==1400
                irow=irow+1;
                Value_1400(irow) = stats_all{irat}.wod_origin_depth(itrial);
            end
        end
    end
    x = 1.2 * ones(1, length(Value_1400));
    scatter(x, Value_1400,70, 'ko',  'LineWidth', 1,'MarkerFaceColor','w','jitter', 'on', 'jitterAmount', 0.05);
    % Set up axes.
    xlim([0.5, 1.5]);
    ylim([700, 2000]);
    ax = gca;
    ax.XTick = [1, 1.2];
    ax.XTickLabels = {'1000','1400'};
    
    fname=fullfile(config16{4}.imagesavedir,'delays','boxplots','scatter_alldepths');
    dtx_savefigure(figure_depths,fname,'pdf','png','close');
    clear fname Value_1000 Value_1400 x
    %% GET DATA for plotting and filtering
    config=[config16 config32];
    
    for irat= 1:size(config,2)
        
        
        if isempty(config{irat})
            continue
        end
        
        
        temp=load(fullfile(config{irat}.datasavedir,sprintf('%sLFP_WoD',config{irat}.prefix)));
        LFP{irat}=temp.LFP{1}.WoD;
        clear temp
        
        MuseStruct{irat}=readMuseMarkers(config{irat},false);
        
        %rename channels according to depth
        for ichan = 1:size(config{irat}.LFP.channel, 2)
            idx = strcmp(config{irat}.LFP.channel{ichan}, LFP{irat}.label);
            label_renamed{idx} = config{irat}.LFP.rename{ichan};
        end
        LFP{irat}.label = label_renamed';
        clear label_renamed
        
        %remove breathing and ekg channel
        cfgtemp         = [];
        cfgtemp.channel = {'all', '-E0', '-Respi', '-ECG','-Puff'};
        LFP{irat}             = ft_selectdata(cfgtemp, LFP{irat});
        LFP_cleaned{irat}     = LFP{irat}; %save for later removing of artefacts
        
        %low-pass filter data
        cfgtemp=[];
        cfgtemp.lpfilter='yes';
        cfgtemp.lpfreq= config{irat}.LFP.lpfilter_wod_detection;
        cfgtemp.lpfilttype='fir';
        LFP_lpfilt{irat}= ft_preprocessing(cfgtemp,LFP_cleaned{irat});
        
        
    end %irat
    %% Overdraws
    for irat=1:size(config,2)
        if isempty(config{irat})
            continue
        end
        chan_list=LFP_lpfilt{irat}.label;
        
        for itrial=1:size(LFP{irat}.trial,2)
            %% Overdraws channels to origin with colors
            
            %find initiation channel index
            idx_origin= find(stats_all{irat}.Depth(:,itrial)==stats_all{irat}.wod_origin_depth(itrial));
            
            starttrial              = LFP_lpfilt{irat}.trialinfo.begsample(itrial) / LFP_lpfilt{irat}.fsample;
            endtrial                = LFP_lpfilt{irat}.trialinfo.endsample(itrial) / LFP_lpfilt{irat}.fsample;
            offsettrial             = LFP_lpfilt{irat}.trialinfo.offset(itrial) / LFP_lpfilt{irat}.fsample;
            
            v_off= MuseStruct{irat}{1}{1}.markers.Vent_Off.synctime(itrial)-starttrial+offsettrial;
            v_on= MuseStruct{irat}{1}{1}.markers.Vent_On.synctime(itrial)-starttrial+offsettrial;
            
            C_1= winter(idx_origin);
            C_2= winter(size(chan_list,1)-idx_origin+1);
            C_2(end,:)=[];
            C_all=vertcat(C_1,flip(C_2));
            
            %FIXME FAIRE NORMALISATION AMPLITUDE DES WAVESs
            
            fig_overdraw_wod=figure;
            sgtitle(sprintf('Rat_%i_WoD_%i_on_%i_overdraw',irat,itrial,size(LFP{irat}.trial,2)), 'interpreter', 'none')
            for ichan=1:size(chan_list,1)
                LFP_color_wod{irat}=LFP_lpfilt{irat};
                LFP_color_wod{irat}.trial{itrial}(ichan,:)= movmean(LFP_color_wod{irat}.trial{itrial}(ichan,:),1000);
                %baseline correction and normalisation
                idx_wod=LFP_color_wod{irat}.time{itrial}== floor(MuseStruct{irat}{1}{1}.markers.WOD.synctime(itrial)-starttrial+offsettrial-10);
                baseline(ichan)= LFP_color_wod{irat}.trial{itrial}(ichan,idx_wod);
                LFP_color_wod{irat}.trial{itrial}(ichan,:)=((LFP_color_wod{irat}.trial{itrial}(ichan,:)-baseline(ichan))./ -stats_all{irat}.WoD.peak_value(ichan,itrial));
                
                
                plot(LFP_color_wod{irat}.time{itrial},LFP_color_wod{irat}.trial{itrial}(ichan,:),'Color',C_all(ichan,:))
                hold on
                
                xlim([min(stats_all{irat}.WoD.peak_time(:,itrial))-10 max(stats_all{irat}.WoD.peak_time(:,itrial))+20])
            end %ichan
            fname_overdraw_wod=fullfile(config{irat}.imagesavedir,'overdraw_traces','color_from_origin',sprintf('%s_WoD_%i_on_%i',config{irat}.prefix,itrial,size(LFP{irat}.trial,2)));
            
            if ~isfolder(fullfile(config{irat}.imagesavedir,'overdraw_traces','color_from_origin'))
                mkdir(fullfile(config{irat}.imagesavedir,'overdraw_traces','color_from_origin'));
            end
            
            dtx_savefigure(fig_overdraw_wod,fname_overdraw_wod,'pdf','png','close');
            
            
            if config{irat}.LFP.recov{itrial}==1
                fig_overdraw_wor=figure;
                sgtitle(sprintf('Rat_%i_WoD_%i_on_%i_overdraw',irat,itrial,size(LFP{irat}.trial,2)))
                for ichan=1:size(chan_list,1)
                    LFP_color_wor{irat}=LFP_lpfilt{irat};
                    LFP_color_wor{irat}.trial{itrial}(ichan,:)= movmean(LFP_color_wor{irat}.trial{itrial}(ichan,:),1000);
                    %baseline correction and normalisation
                    idx_wor=LFP_color_wor{irat}.time{itrial}== ~(MuseStruct{irat}{1}{1}.markers.WOR.synctime(itrial)-starttrial+offsettrial);
                    baseline= LFP_color_wor{irat}.trial{itrial}(ichan,idx_wor);
                    LFP_color_wor{irat}.trial{itrial}=(LFP_color_wor{irat}.trial{itrial})./stats_all{irat}.WoR.peak_value(ichan,itrial);
                    
                    plot(LFP_color_wor{irat}.time{itrial},LFP_color_wor{irat}.trial{itrial}(ichan,:),'Color',C_all(ichan,:))
                    hold on
                    
                    xlim([min(stats_all{irat}.WoR.peak_time(:,itrial))-10 max(stats_all{irat}.WoR.peak_time(:,itrial))+20])
                end %ichan
                fname_overdraw_wor=fullfile(config{irat}.imagesavedir,'overdraw_traces','color_from_origin',sprintf('%s_WoR_%i_on_%i',config{irat}.prefix,itrial,size(LFP{irat}.trial,2)));
                
                if ~isfolder(fullfile(config{irat}.imagesavedir,'overdraw_traces','color_from_origin'))
                    mkdir(fullfile(config{irat}.imagesavedir,'overdraw_traces','color_from_origin'));
                end
                xlim([v_on v_on+80]);
                
                dtx_savefigure(fig_overdraw_wor,fname_overdraw_wor,'pdf','png','close');
            end
            
            
            clear C_1 C_2 C_all
            
            
            %% Overdraw all traces
            
            %WOD
            fig_multitraces_wod = figure;
            sgtitle(sprintf('Rat_%i_WoD_%i_on_%i_multitraces',irat,itrial,size(LFP{irat}.trial,2)))
            d=6000;
            a=0;
            for ichan=1:size(chan_list,1)
                
                plot(LFP_lpfilt{irat}.time{itrial},LFP_lpfilt{irat}.trial{itrial}(ichan,:)-(a*d),'Color','k','LineWidth',2);
                hold on
                ylabel(sprintf('%i',config{irat}.LFP.chan_depth{ichan}));
                
                a= a+1;
                
            end %ichan
            xlim([v_off+60 v_on+10]);
            xline(v_off);
            xline(v_on);
            
            fname_multitraces=fullfile(config{irat}.imagesavedir,'overdraw_traces','multitraces',sprintf('%s_WoD_%i_on_%i',config{irat}.prefix,itrial,size(LFP{irat}.trial,2)));
            
            if ~isfolder(fullfile(config{irat}.imagesavedir,'overdraw_traces','multitraces'))
                mkdir(fullfile(config{irat}.imagesavedir,'overdraw_traces','multitraces'));
            end
            
            dtx_savefigure(fig_multitraces_wod,fname_multitraces,'pdf','png','close');
            
            % WOR
            fig_multitraces_wor = figure;
            sgtitle(sprintf('Rat_%i_WoD_%i_on_%i_multitraces',irat,itrial,size(LFP{irat}.trial,2)))
            d=3000;
            a=0;
            for ichan=1:size(chan_list,1)
                
                plot(LFP_lpfilt{irat}.time{itrial},LFP_lpfilt{irat}.trial{itrial}(ichan,:)-(a*d),'Color','k','LineWidth',2);
                hold on
                ylabel(sprintf('%i',config{irat}.LFP.chan_depth{ichan}));
                
                a= a+1;
                
            end %ichan
            xlim([v_on v_on+60]);
            xline(v_on);
            
            fname_multitraces=fullfile(config{irat}.imagesavedir,'overdraw_traces','multitraces',sprintf('%s_WoR_%i_on_%i',config{irat}.prefix,itrial,size(LFP{irat}.trial,2)));
            
            if ~isfolder(fullfile(config{irat}.imagesavedir,'overdraw_traces','multitraces'))
                mkdir(fullfile(config{irat}.imagesavedir,'overdraw_traces','multitraces'));
            end
            
            dtx_savefigure(fig_multitraces_wor,fname_multitraces,'pdf','png','close');
            
            
            
        end %itrial
    end %irat
    %% Overdraw origin channels WoD separated by class
    
    fig_origin_deep=figure;hold
    for irat=1:size(config,2)
        if isempty(config{irat})
            continue
        end
        
        Color_ori={'k',[0.8500 0.3250 0.0980]};
        data_plot=LFP_lpfilt{irat};
        
        for itrial= 1:size(LFP_lpfilt{irat}.trial,2)
            idx_ori= stats_all{irat}.Depth(:,itrial)==stats_all{irat}.wod_origin_depth(itrial);
            if stats_all{irat}.oridepthclass(itrial)== 1400
                plot(data_plot.time{itrial}-stats_all{irat}.WoD.peak_time(idx_ori,itrial),data_plot.trial{itrial}(idx_ori,:),'Color',[0.8500 0.3250 0.0980 0.5],'LineWidth',1);
                xlim([-10;10]);
            end
        end %itrial
    end %irat
    overdraw_oripath= fullfile(config{5}.imagesavedir,'overdraw_traces','origin_channel');
    
    if ~isfolder(overdraw_oripath)
        mkdir(overdraw_oripath);
    end
    
    fname_deep=fullfile(overdraw_oripath,'deep_origin');
    dtx_savefigure(fig_origin_deep,fname_deep,'pdf','png','close');
    
    
    fig_origin_sup=figure;hold
    for irat=1:size(config,2)
        if isempty(config{irat})
            continue
        end
        data_plot=LFP_lpfilt{irat};
        
        for itrial= 1:size(LFP_lpfilt{irat}.trial,2)
            idx_ori= stats_all{irat}.Depth(:,itrial)==stats_all{irat}.wod_origin_depth(itrial);
            
            if stats_all{irat}.oridepthclass(itrial)==1000
                plot(data_plot.time{itrial}-stats_all{irat}.WoD.peak_time(idx_ori,itrial),data_plot.trial{itrial}(idx_ori,:),'Color',[0 0 0 0.5],'LineWidth',1);
                xlim([-10;10]);
            end
        end %itrial
    end %irat
    
    fname_sup=fullfile(overdraw_oripath,'superficial_origin');
    dtx_savefigure(fig_origin_sup,fname_sup,'pdf','png','close');
    %% Distribution of origin depth
    cfg=[config16 config32];
    
    irow=0;
    for irat=1:size(cfg,2)
        
        if isempty(cfg{irat})
            continue
        end
        
        for itrial=1:size(stats_all{irat}.wod_origin_depth,2)
            irow=irow+1;
            if stats_all{irat}.oridepthclass(itrial)==1000
                origin_depth{1}(irow)=stats_all{irat}.wod_origin_depth(itrial);
            elseif stats_all{irat}.oridepthclass(itrial)==1400
                origin_depth{2}(irow)=stats_all{irat}.wod_origin_depth(itrial);
            end
            
        end %itrial
    end %irat
    
    %delete zeros by nans
    for i= 1:size(origin_depth,2)
        ZerosIdx=origin_depth{i}==0;
        origin_depth{i}(ZerosIdx)=[];
    end
    
    Color_hist={'k',[0.8500 0.3250 0.0980]};
    
    fig=figure;hold
    
    for i=1:size(origin_depth,2)
        h=histfit(origin_depth{i}',33,'wbl');
        set(h(2),'color',Color_hist{i});
    end
    xlim([600 2000]);
    fname_h=fullfile(cfg{irat}.imagesavedir,'delays','distrib','oridepth');
    
    if ~isfolder(fullfile(cfg{4}.imagesavedir,'delays','distrib'))
        mkdir(fullfile(cfg{4}.imagesavedir,'delays','distrib'));
    end
    
    dtx_savefigure(fig,fname_h,'png','pdf','close');
    fitmethis(origin_depth{2})
end %rat_list

end %wod_project