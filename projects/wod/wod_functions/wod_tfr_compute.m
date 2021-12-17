function wod_tfr_compute(cfg,MuseStruct,LFP)

%

if isempty(cfg)
    return
end

MuseStruct = concatenateMuseMarkers(cfg,MuseStruct,false);

%remove physio constant channels
cfgtemp         = [];
cfgtemp.channel = {'all', '-E0', '-Respi', '-ECG'};
LFP             = ft_selectdata(cfgtemp, LFP);
LFP_cleaned     = LFP; %save for later removing of artefacts


%hp filter lfp to exclude WoD and WoR and Notch 50Hz
cfgtemp             = [];
cfgtemp.hpfilter    = 'yes';
cfgtemp.hpfilttype  = 'fir';
cfgtemp.hpfreq      = cfg.LFP.hpfilter_wod_exclusion;
cfgtemp.bsfilter     = 'yes';
cfgtemp.bsfilttype     = 'fir';
cfgtemp.bsfreq          = [49 51];
LFP_cleaned          = ft_preprocessing(cfgtemp, LFP_cleaned);


%filter lfp to better recognize WOD/WOR peak
cfgtemp             = [];
cfgtemp.lpfilter    = 'yes';
cfgtemp.lpfilttype  = 'fir';

cfgtemp.lpfreq      = cfg.LFP.lpfilter_wod_detection;
LFP_lpfilt          = ft_preprocessing(cfgtemp, LFP_cleaned);

for itrial= 1:size(LFP.trial,2)

 %select one trial
    cfgtemp         = [];
    cfgtemp.trials  = itrial;
    LFP_trial       = ft_selectdata(cfgtemp, LFP_cleaned);

%FIXME
%correction approté par sofia grace à la correction préccédente de Paul le 01/12/2021    
%     %recover trial real timings to use it with muse markers
%     starttrial              = LFP_trial.trialinfo.begsample / LFP_trial.fsample;%temps en sec
%     endtrial                = LFP_trial.trialinfo.endsample / LFP_trial.fsample;
%     offsettrial             = LFP_trial.trialinfo.offset / LFP_trial.fsample;
   if size(MuseStruct{1}.starttime,2) == 1 || itrial == 1
        
        starttrial( itrial)              = LFP_lpfilt.trialinfo.begsample(itrial) / LFP_lpfilt.fsample;
        endtrial(itrial)                = LFP_lpfilt.trialinfo.endsample(itrial) / LFP_lpfilt.fsample;
        offsettrial( itrial)             = LFP_lpfilt.trialinfo.offset(itrial) / LFP_lpfilt.fsample;
        
    else %un cas particulier si plusieurs dossiers (au lieu de données concaténées)
        
        idir                            = LFP_lpfilt.trialinfo.idir(itrial);
        length_previous_dir = 0;
        for i = 2:idir
            length_previous_dir = length_previous_dir + seconds(MuseStruct{1}.endtime{idir-1} - MuseStruct{1}.starttime{idir-1});
        end
        starttrial(itrial)           = LFP_lpfilt.trialinfo.begsample(itrial) / LFP_lpfilt.fsample + length_previous_dir;
        endtrial(itrial)             = LFP_lpfilt.trialinfo.endsample(itrial) / LFP_lpfilt.fsample + length_previous_dir;
        offsettrial(itrial)         = LFP_lpfilt.trialinfo.offset(itrial) / LFP_lpfilt.fsample;
        
    end





    
    
        
        %do time frequency analysis
        cfgtemp                         = [];
        cfgtemp.channel                 = 'all';
        cfgtemp.method                  = 'mtmconvol';
        cfgtemp.output                  = 'pow';
        cfgtemp.taper                   = 'dpss'; %default = dpss
        cfgtemp.pad                     = 'nextpow2';
        cfgtemp.keeptrials              = 'no';
        cfgtemp.tapsmofrq               = cfg.timefreq.tapsmofrq;
        cfgtemp.foi                     = cfg.timefreq.foi;
        cfgtemp.t_ftimwin               = ones(size(cfgtemp.foi))*cfg.timefreq.t_ftimwin;
        cfgtemp.toi                     = cfg.timefreq.toi(1) : cfg.timefreq.timestep : cfg.timefreq.toi(end);
        timefreq_alldata{itrial}        = ft_freqanalysis(cfgtemp,LFP_trial);
        
        %replace artifacts by nans
        %need to remove artefacts after time freq analysis, because
        %any nan in the lfp data creates a time freq with only nan values
        if isfield(MuseStruct{1}.markers, 'BAD__START__')
            if isfield(MuseStruct{1}.markers.BAD__START__, 'synctime')
                %get bad timings
                bad_start                   = MuseStruct{1}.markers.BAD__START__.synctime;
                bad_end                     = MuseStruct{1}.markers.BAD__END__.synctime;
                if length(bad_start) ~= length(bad_end)
                    error('not the same numbers of bad start and end markers');
                end
                %                         t_tfr                   = timefreq_alldata{itrial}.time;
                t_lfp                   = LFP_trial.time{1};
                t_tfr                   = timefreq_alldata{itrial}.time;
                bad_sel                 = find(bad_start >= starttrial(itrial) & bad_start <= endtrial(itrial));
                %bad_sel2                 = find(bad_start >= starttrial & bad_start <= endtrial);
                %go through each bad timing
                
                %FIXME
                %correction apporté par Sofia le 7/12/2021
%                 for ibad = bad_sel
%                     %remove lfp artefacts
%                     bad_period_lfp = t_lfp >= (bad_start(ibad)- starttrial + offsettrial) & t_lfp <= (bad_end(ibad)- starttrial + offsettrial);
%                     LFP_cleaned.trial{itrial}(:,bad_period_lfp) = nan(size(LFP_cleaned.trial{itrial},1),sum(bad_period_lfp));
%                     %remove tfr artefacts: all window with at least one artefacted sample
%                     bad_period_tfr = t_tfr >= (bad_start(ibad)- starttrial + offsettrial) - cfg.timefreq.t_ftimwin/2 & t_tfr <= (bad_end(ibad)- starttrial + offsettrial) + cfg.timefreq.t_ftimwin/2;
%                     timefreq_alldata{itrial}.powspctrm(:,:,bad_period_tfr) = nan(size(timefreq_alldata{itrial}.powspctrm,1),size(timefreq_alldata{itrial}.powspctrm,2),sum(bad_period_tfr));
%                 end
               
                  %voir avec Paul si c'est la bonne correction
                  for bad_number=1:length(bad_sel)
                      for ibad = bad_sel(bad_number)
                          %remove lfp artefacts
                          bad_period_lfp = t_lfp >= (bad_start(ibad)- starttrial(itrial) + offsettrial(itrial)) & t_lfp <= (bad_end(ibad)- starttrial(itrial) + offsettrial(itrial));
                          LFP_cleaned.trial{itrial}(:,bad_period_lfp) = nan(size(LFP_cleaned.trial{itrial},1),sum(bad_period_lfp));
                          %remove tfr artefacts: all window with at least one artefacted sample
                          bad_period_tfr = t_tfr >= (bad_start(ibad)- starttrial(itrial) + offsettrial(itrial)) - cfg.timefreq.t_ftimwin/2 & t_tfr <= (bad_end(ibad)- starttrial(itrial) + offsettrial(itrial)) + cfg.timefreq.t_ftimwin/2;
                          timefreq_alldata{itrial}.powspctrm(:,:,bad_period_tfr) = nan(size(timefreq_alldata{itrial}.powspctrm,1),size(timefreq_alldata{itrial}.powspctrm,2),sum(bad_period_tfr));
                      end
                  end

                clear bad_sel t_tfr t_lfp
            end
        end
        
        for ichan = 1:size(timefreq_alldata{itrial}.label,1)

        %% select channel for short param
        ichan_name              = timefreq_alldata{itrial}.label{ichan};
         %ichan_name= LFP_lpfilt.label{ichan};
        cfgtemp                 = [];
        cfgtemp.channel         = ichan_name;
        timefreq_ichan_temp   	= ft_selectdata(cfgtemp,timefreq_alldata{itrial});
        %timefreq_ichan_temp   	= ft_selectdata(cfgtemp,LFP_lpfilt);

        %% WOD DATA : find WOD peak per channel, and normalize time per channel
        %use filtered data to find wod
        
        timefreq_wod{itrial}.(ichan_name)            = timefreq_ichan_temp;
        
        %select lfp channel with the same name as tfr channel (in
        %case channel numbers were schuffled by fieldtrip)
        chan_idx    = strcmp(LFP_lpfilt.label, ichan_name);
        
%         %get hand-annotated wod timing
%         %wod_marker = MuseStruct{1}.markers.WOD.synctime(itrial);%sofia changed :  MuseStruct{1} to  MuseStruct{1}{itrial} & synctime(itrial)to synctime
%         %%%sofia a rajouter la converstion suivante le 21/10/2021
%         wod_marker = MuseStruct{1}.markers.WOD.synctime(itrial);
%         temp_wod_marker = [];
%         if itrial>1
%         temp_wod_marker =((LFP_trial.sampleinfo(itrial)/MuseStruct{1, 1}.markers.Stopping_Recording.synctime(itrial))* wod_marker)/LFP_trial.fsample
%         wod_marker = temp_wod_marker
%         fname_out = fullfile(cfg.datasavedir,'sampleinfo', sprintf([cfg.prefix,'trial' ,num2str(itrial), '.mat']))
%         sampleinfo = LFP_trial.sampleinfo(itrial)
%         save(fname_out,'sampleinfo')
%         end
%         %%% fin de la modification 
        
       %correction apporté par sofia le 07/12/2021
        wod_marker= MuseStruct{1}.markers.WOD.synctime(itrial) ;

        %select times where to search WOD peak
        t      = LFP_lpfilt.time{itrial};%{1}; sofia changed time{itrial} time{1}
        t_1    = t > (wod_marker + cfg.LFP.wod_toisearch(1) - starttrial(itrial) + offsettrial(itrial));
        t_2    = t < (wod_marker + cfg.LFP.wod_toisearch(2) - starttrial(itrial) + offsettrial(itrial)) ;
        t_sel  = t_1 & t_2;
        %Search LFP maximum peak in this selected window. '-'LFP because wod is negative
        [v_peak_wod, t_peak_wod] = findpeaks(-LFP_lpfilt.trial{itrial}(chan_idx,t_sel),t(t_sel),'NPeaks',1,'SortStr','descend','WidthReference','Halfheight');
        clear t t_1 t_2 t_sel
        
        %keep only data between 0 and wod
        cfgtemp                                   = [];
        cfgtemp.latency                           = [0 t_peak_wod];
        timefreq_wod{itrial}.(ichan_name)  = ft_selectdata(cfgtemp,timefreq_wod{itrial}.(ichan_name));
        
        %normalize time
        timefreq_wod_timenorm{itrial}.(ichan_name)        = timefreq_wod{itrial}.(ichan_name);
        timefreq_wod_timenorm{itrial}.(ichan_name).time   = timefreq_wod{itrial}.(ichan_name).time ./ t_peak_wod;
        %check the location of the peak detection
        %                 fig= figure;hold
        %                 plot(LFP_trial_lpfilt.time{1}, LFP_trial_lpfilt.trial{1}(ichan,:));
        %                 scatter(t_peak_wod, -v_peak_wod,50,'xr');
        %                 xlim([t_peak_wod-10 t_peak_wod+10]);
        %                 print(fig, '-dpng', fullfile(cfg.imagesavedir,'Detection',sprintf('Rat%g_WOD%g_%s.png',irat,itrial,ichan_name)),'-r600');
        %                 close all
        
        for ifreq = 1:size(timefreq_wod_timenorm{itrial}.(ichan_name).powspctrm,2)
            %resample data to have the same number of data points, for
            %time-normalized data
            t_old                                                = timefreq_wod_timenorm{itrial}.(ichan_name).time;
            t_new                                                = linspace(0,1,1000);
            %powspctrm_new(1,1,:)                                = pchip(t_old,squeeze(timefreq_wod_timenorm{itrial}.(ichan_name).powspctrm(1,1,:)),t_new);
            powspctrm_new(1,ifreq,:)                             = interp1(t_old,squeeze(timefreq_wod_timenorm{itrial}.(ichan_name).powspctrm(1,ifreq,:)),t_new);
        end
        timefreq_wod_timenorm{itrial}.(ichan_name).time         = t_new;
        timefreq_wod_timenorm{itrial}.(ichan_name).powspctrm    = powspctrm_new;
        clear powspctrm_new
        %plot(timefreq_wod{itrial}.time{1},timefreq_wod{itrial}.trial{1}); xlim([0 0.95])
        %% baseline short
        
        cfgtemp = [];
        cfgtemp.latency = [-300 0];
        timefreq_baseline{itrial}.(ichan_name) = ft_selectdata(cfgtemp,timefreq_ichan_temp);
        %% MAKE BASELINE CORRECTION FOR SHORT PERIODS
        
        timefreq_baseline_blcorrected{itrial}.(ichan_name)= timefreq_baseline{itrial}.(ichan_name);
        timefreq_wod_blcorrected{itrial}.(ichan_name)= timefreq_wod{itrial}.(ichan_name);
        timefreq_wod_timenorm_blcorrected{itrial}.(ichan_name) = timefreq_wod_timenorm{itrial}.(ichan_name);
        
        for ifreq = 1:size(timefreq_wod{itrial}.(ichan_name).freq,2) %baseline
            baseline_ifreq = nanmean(squeeze(timefreq_baseline{itrial}.(ichan_name).powspctrm(1,ifreq,:))); %baseline
            timefreq_baseline_blcorrected{itrial}.(ichan_name).powspctrm(1,ifreq,:) = timefreq_baseline{itrial}.(ichan_name).powspctrm(1,ifreq,:) ./ baseline_ifreq; %baseline
            timefreq_wod_blcorrected{itrial}.(ichan_name).powspctrm(1,ifreq,:) = timefreq_wod{itrial}.(ichan_name).powspctrm(1,ifreq,:) ./ baseline_ifreq;
            timefreq_wod_timenorm_blcorrected{itrial}.(ichan_name).powspctrm(1,ifreq,:) = timefreq_wod_timenorm{itrial}.(ichan_name).powspctrm(1,ifreq,:) ./ baseline_ifreq;
        end %baseline
        
        clear timefreq_ichan_temp
        %% Make Logarythm for short periods
        
        %duplicate fieldtrip structure
        log_timefreq_baseline_blcorrected{itrial}.(ichan_name)= timefreq_baseline_blcorrected{itrial}.(ichan_name);
        log_timefreq_baseline{itrial}.(ichan_name)=timefreq_baseline{itrial}.(ichan_name);
        log_timefreq_wod_blcorrected{itrial}.(ichan_name)=timefreq_wod_blcorrected{itrial}.(ichan_name);
        log_timefreq_wod{itrial}.(ichan_name)=timefreq_wod{itrial}.(ichan_name);
        log_timefreq_wod_timenorm{itrial}.(ichan_name)=timefreq_wod_timenorm{itrial}.(ichan_name);
        log_timefreq_wod_timenorm_blcorrected{itrial}.(ichan_name)=timefreq_wod_timenorm_blcorrected{itrial}.(ichan_name);
        
        %apply logarythm to powspctrm
        log_timefreq_baseline_blcorrected{itrial}.(ichan_name).powspctrm= log10(timefreq_baseline_blcorrected{itrial}.(ichan_name).powspctrm);
        log_timefreq_baseline{itrial}.(ichan_name).powspctrm= log10(timefreq_baseline{itrial}.(ichan_name).powspctrm);
        log_timefreq_wod_blcorrected{itrial}.(ichan_name).powspctrm= log10(timefreq_wod_blcorrected{itrial}.(ichan_name).powspctrm);
        log_timefreq_wod{itrial}.(ichan_name).powspctrm= log10(timefreq_wod{itrial}.(ichan_name).powspctrm);
        log_timefreq_wod_timenorm{itrial}.(ichan_name).powspctrm= log10(timefreq_wod_timenorm{itrial}.(ichan_name).powspctrm);
        log_timefreq_wod_timenorm_blcorrected{itrial}.(ichan_name).powspctrm= log10(timefreq_wod_timenorm_blcorrected{itrial}.(ichan_name).powspctrm);
        %% SHORT remove cfg fields as it is what takes the most of place on disk, whereas we do not use it later
        timefreq_wod{itrial}.(ichan_name)            = rmfield(timefreq_wod{itrial}.(ichan_name),{'cfg'});
        timefreq_wod_blcorrected{itrial}.(ichan_name)            = rmfield(timefreq_wod_blcorrected{itrial}.(ichan_name),{'cfg'});
        timefreq_wod_timenorm{itrial}.(ichan_name) 	 = rmfield(timefreq_wod_timenorm{itrial}.(ichan_name),{'cfg'});
        timefreq_wod_timenorm_blcorrected{itrial}.(ichan_name) 	 = rmfield(timefreq_wod_timenorm_blcorrected{itrial}.(ichan_name),{'cfg'});
        timefreq_baseline{itrial}.(ichan_name)                  = rmfield(timefreq_baseline{itrial}.(ichan_name),{'cfg'});
        timefreq_baseline_blcorrected{itrial}.(ichan_name)       = rmfield(timefreq_baseline_blcorrected{itrial}.(ichan_name),{'cfg'});
        log_timefreq_baseline_blcorrected{itrial}.(ichan_name)    = rmfield(log_timefreq_baseline_blcorrected{itrial}.(ichan_name),{'cfg'});
        log_timefreq_baseline{itrial}.(ichan_name)                = rmfield(log_timefreq_baseline{itrial}.(ichan_name),{'cfg'});
        log_timefreq_wod_blcorrected{itrial}.(ichan_name)         = rmfield(log_timefreq_wod_blcorrected{itrial}.(ichan_name),{'cfg'});
        log_timefreq_wod{itrial}.(ichan_name)                     = rmfield(log_timefreq_wod{itrial}.(ichan_name),{'cfg'});
        log_timefreq_wod_timenorm{itrial}.(ichan_name)            = rmfield(log_timefreq_wod_timenorm{itrial}.(ichan_name),{'cfg'});
        log_timefreq_wod_timenorm_blcorrected{itrial}.(ichan_name)= rmfield(log_timefreq_wod_timenorm_blcorrected{itrial}.(ichan_name),{'cfg'});
end %ichan
end %itrial

% add empty missing channels channels to have the same channels between rats
for itrial = 1:size( timefreq_wod,2)
    for chan_name = string(cfg.LFP.allchannel)
        %chan_renamed = sprintf('E%d',str2num(regexp(chan_name,'\d*','Match')));
        if ~isfield(timefreq_wod{itrial}, chan_name)
           
            timefreq_wod{itrial}.(chan_name)                 = [];
            timefreq_wod_blcorrected{itrial}.(chan_name)                 = [];
            
            timefreq_wod_timenorm{itrial}.(chan_name)        = [];
            timefreq_wod_timenorm_blcorrected{itrial}.(chan_name)                 = [];
         
            timefreq_baseline{itrial}.(chan_name)            = [];
            timefreq_baseline_blcorrected{itrial}.(chan_name)            = [];
 
            log_timefreq_wod{itrial}.(chan_name)                 = [];
            log_timefreq_wod_blcorrected{itrial}.(chan_name)                 = [];
            
            log_timefreq_wod_timenorm{itrial}.(chan_name)        = [];
            log_timefreq_wod_timenorm_blcorrected{itrial}.(chan_name)                 = [];
          
            log_timefreq_baseline{itrial}.(chan_name)            = [];
            log_timefreq_baseline_blcorrected{itrial}.(chan_name)            = [];
            
        end
    end
end

%save time freq data to disk :

%WOD data (entre vent_off et pic de wod), power normalized with baseline period :
save(fullfile(cfg.datasavedir,[cfg.prefix, 'timefreq_wod.mat']),'timefreq_wod','-v7.3');
save(fullfile(cfg.datasavedir,[cfg.prefix, 'timefreq_wod_blcorrected.mat']),'timefreq_wod_blcorrected','-v7.3');
save(fullfile(cfg.datasavedir,[cfg.prefix, 'log_timefreq_wod.mat']),'log_timefreq_wod','-v7.3');
save(fullfile(cfg.datasavedir,[cfg.prefix, 'log_timefreq_wod_blcorrected.mat']),'log_timefreq_wod_blcorrected','-v7.3');

%WOD data : time normalized between vent_off and wod peak, per channel. power normalized with baseline period :
save(fullfile(cfg.datasavedir,[cfg.prefix, 'timefreq_wod_timenorm.mat']),'timefreq_wod_timenorm','-v7.3');
save(fullfile(cfg.datasavedir,[cfg.prefix, 'timefreq_wod_timenorm_blcorrected.mat']),'timefreq_wod_timenorm_blcorrected','-v7.3');
save(fullfile(cfg.datasavedir,[cfg.prefix, 'log_timefreq_wod_timenorm.mat']),'log_timefreq_wod_timenorm','-v7.3');
save(fullfile(cfg.datasavedir,[cfg.prefix, 'log_timefreq_wod_timenorm_blcorrected.mat']),'log_timefreq_wod_timenorm_blcorrected','-v7.3');

%Baseline data: t0 at Vent_Off analysis made before Vent_Off
save(fullfile(cfg.datasavedir,[cfg.prefix, 'timefreq_baseline.mat']),'timefreq_baseline','-v7.3');
save(fullfile(cfg.datasavedir,[cfg.prefix, 'timefreq_baseline_blcorrected.mat']),'timefreq_baseline_blcorrected','-v7.3');
save(fullfile(cfg.datasavedir,[cfg.prefix, 'log_timefreq_baseline.mat']),'log_timefreq_baseline','-v7.3');
save(fullfile(cfg.datasavedir,[cfg.prefix, 'log_timefreq_baseline_blcorrected.mat']),'log_timefreq_baseline_blcorrected','-v7.3');

end % wod_tfr_compute