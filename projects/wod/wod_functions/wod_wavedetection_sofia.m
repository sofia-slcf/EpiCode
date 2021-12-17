function  wod_wavedetection_sofia(irat, cfg, data)%(cfg, force)

% sofia a suprimé la boucle irat puis a rajouté irat comme variable de la fonction
% cfg{irat} est remplacé par cfg tout court

detectsavedir=fullfile(cfg.imagesavedir,'detection');

%Load LFP and Muse markers
temp= load(fullfile(cfg.datasavedir,sprintf('%s%s_%s.mat',cfg.prefix,'LFP',cfg.name{1})));
LFP = temp.LFP{1}.(cfg.LFP.name{1});
clear temp
MuseStruct               = readMuseMarkers(cfg,false);

MuseStruct = concatenateMuseMarkers(cfg,MuseStruct,true); %added by sofia

%determine a criteria to accept or reject a detection
WOD_threshold =  300;

%vérifier qu'il y a bien autant de trials que de marqueurs Vent_Off
%     startmarker = cfg.muse.startmarker.(cfg.LFP.name{1});
%     if size(LFP.trial,2) ~= size(MuseStruct{1}.markers.(startmarker).synctime,2)
%         error('Not the same number of trials that of marker start for %s. \nCheck that begin/end of each trial is not before start of file or after end of file', cfg.prefix(1:end-1));
%     end
%
%rename channels according to depth

for ichan = 1:size(cfg.LFP.channel, 2)
    idx = strcmp(cfg.LFP.channel{ichan}, LFP.label);
    label_renamed{idx} = cfg.LFP.rename{ichan};
end
LFP.label = label_renamed';
clear label_renamed

%remove breathing and ekg channel
cfgtemp         = []; %% ici
cfgtemp.channel = {'all'};%{'all', '-E0', '-Respi', '-ECG'};%%ici
% LFP             = ft_selectdata(cfgtemp, LFP);%modifié par sofia
%    LFP_cleaned     = LFP; %save for later removing of artefacts
LFP_cleaned     = data;
%remove 50Hz and interpolate with 49 and 51 Hz
%   LFP_cleaned= ft_preproc_dftfilter(LFP_cleaned,LFP_cleaned.fsample,50,'Flreplace','neighbour');


%filter lfp to better recognize WOD/WOR peak
cfgtemp             = [];
cfgtemp.lpfilter    = 'yes';
cfgtemp.lpfilttype  = 'fir';
cfgtemp.lpfreq      = cfg.LFP.lpfilter_wod_detection;
LFP_lpfilt          = ft_preprocessing(cfgtemp, LFP_cleaned);

%filter lfp to better recognize silence delay
cfgtemp             = [];
cfgtemp.bpfilter    = 'yes';
cfgtemp.bpfreq      = [1 30];
%     cfgtemp.lpfilter    = 'yes';
%     cfgtemp.lpfilttype  = 'fir';
%     cfgtemp.lpfreq      = 30;
%     cfgtemp.hpfilter    = 'yes';
%     cfgtemp.hpfilttype  = 'fir';
%     cfgtemp.hpfreq      = 1;
%     cfgtemp.bsfilter    = 'yes';
%     cfgtemp.bsfreq      = [49 51];
LFP_hpfilt          = ft_preprocessing(cfgtemp, LFP_cleaned);

for itrial = 1:size(LFP.trial,2)
    
    %recover trial real timings to use it with muse markers
    if size(MuseStruct{1}.starttime,2) == 1 || itrial == 1
        
        starttrial              = LFP_lpfilt.trialinfo.begsample / LFP_lpfilt.fsample;
        endtrial                = LFP_lpfilt.trialinfo.endsample / LFP_lpfilt.fsample;
        offsettrial             = LFP_lpfilt.trialinfo.offset / LFP_lpfilt.fsample;
        
    else %un cas particulier si plusieurs dossiers (au lieu de données concaténées)
        
        idir                            = LFP_lpfilt.trialinfo.idir(itrial);
        length_previous_dir = 0;
        for i = 2:idir
            length_previous_dir = length_previous_dir + seconds(MuseStruct{1}.endtime{idir-1} - MuseStruct{1}.starttime{idir-1});
        end
        starttrial(itrial)              = LFP_lpfilt.trialinfo.begsample(itrial) / LFP_lpfilt.fsample + length_previous_dir;
        endtrial(itrial)                = LFP_lpfilt.trialinfo.endsample(itrial) / LFP_lpfilt.fsample + length_previous_dir;
        offsettrial(itrial)             = LFP_lpfilt.trialinfo.offset(itrial) / LFP_lpfilt.fsample;
        
    end
    
    for ichan=1:size(LFP.label,1)
        
        ichan_name              = LFP_lpfilt.label{ichan};
        
        fprintf('Launching peak detections for %s trial %i and channel %s\n', cfg.prefix,itrial,ichan_name);
        
        
        
        %smoothing of signal with movmean
        LFP_lpfilt.trial{itrial}(ichan,:)= movmean(LFP_lpfilt.trial{itrial}(ichan,:),1000);
        
        
        
        
        
        %% WOD and WOR peak detection
        
        %             WOD detection
        %             select lfp channel (in
        %             case channel numbers were schuffled by fieldtrip)
        chan_idx    = strcmp(LFP_lpfilt.label, ichan_name);
        
        wod_marker = MuseStruct{1}.markers.WOD.synctime(itrial);
        
        %             %%%sofia a rajouter la converstion suivante le 21/10/2021
        %             temp_wod_marker = [];
        %             if itrial>1
        %                 LFP_sample=load(fullfile(cfg.datasavedir,'sampleinfo', sprintf([cfg.prefix,'trial' ,num2str(itrial), '.mat'])));
        %                 temp_wod_marker =((LFP_sample.sampleinfo/MuseStruct{1, 1}.markers.Stopping_Recording.synctime(itrial))* wod_marker)/LFP.fsample;
        %                 wod_marker = temp_wod_marker;
        %             end
        %             %%% fin de la modification
        
        wod_markertime = wod_marker - starttrial(itrial) + offsettrial(itrial);
        
        %select times where to search WOD peak
        t = LFP_lpfilt.time{itrial};
        t_1 = t > (wod_marker + cfg.LFP.wod_toisearch(1) - starttrial(itrial) + offsettrial(itrial));
        t_2 = t < (wod_marker + cfg.LFP.wod_toisearch(2) - starttrial(itrial) + offsettrial(itrial));
        t_sel = t_1 & t_2;
        
        [v_peak_wod, t_peak_wod] = findpeaks(-LFP_lpfilt.trial{itrial}(chan_idx,t_sel),t(t_sel),'NPeaks',1,'SortStr','descend','WidthReference','Halfheight');
        clear t t_1 t_2 t_sel
        
        
        %WOR detection
        wor_marker = MuseStruct{1}.markers.WOR.synctime(itrial);
        %             %modification apportée par sofia: correction synctiime wor
        %             temp_wor_marker = []
        %             if itrial>1
        %             temp_wor_marker =((LFP_sample.sampleinfo/MuseStruct{1, 1}.markers.Stopping_Recording.synctime(itrial))* wor_marker)/LFP.fsample
        %             wor_marker = temp_wor_marker
        %             end
        %             % fin de la modification
        
        
        %select times where to search WOR peak
        
        t = LFP_lpfilt.time{itrial};
        
        t_1 = t > (wor_marker + cfg.LFP.wor_toisearch(1) - starttrial(itrial) + offsettrial(itrial));
        t_2 = t < (wor_marker + cfg.LFP.wor_toisearch(2) - starttrial(itrial) + offsettrial(itrial));
        t_sel = t_1 & t_2;
        
        [v_peak_wor, t_peak_wor] = findpeaks(LFP_lpfilt.trial{itrial}(chan_idx,t_sel),t(t_sel),'NPeaks',1,'SortStr','descend','WidthReference','Halfheight');
        clear t t_1 t_2 t_sel
        
        
        %             store rat name
        stats_all.rat_name= string(cfg.prefix);
        %             store peak timings per channel in structure
        stats_all.WoD.peak_time(ichan,itrial)= t_peak_wod;
        stats_all.WoD.peak_value(ichan,itrial)= -v_peak_wod;
        
        % Calculate Plateau duration
        if cfg.LFP.recov{itrial}==0
            stats_all.plateau_duration(ichan,itrial)=nan;
        else
            stats_all.plateau_duration(ichan,itrial)=t_peak_wor-t_peak_wod;
        end
        
        %             express wor data compared to Vent_On
        t_VentOn= MuseStruct{1}.markers.Vent_On.synctime(itrial)-starttrial(itrial) +offsettrial(itrial);
        
        %             %%%sofia a rajouter la converstion suivante le 21/10/2021
        %             temp_t_VentOn = []
        %             if itrial>1
        %             temp_t_VentOn =((LFP_sample.sampleinfo/MuseStruct{1, 1}.markers.Stopping_Recording.synctime(itrial))* t_VentOn)/LFP.fsample
        %             t_VentOn = temp_t_VentOn
        %             end
        %             %%% fin de la modification
        
        if  cfg.LFP.recov{itrial}==0
            stats_all.WoR.peak_time(ichan,itrial)=nan;
            stats_all.WoR.peak_value(ichan,itrial)= nan;
            
        else
            stats_all.WoR.peak_time(ichan,itrial)= t_peak_wor-t_VentOn;
            stats_all.WoR.peak_value(ichan,itrial)= v_peak_wor;
        end
        
        
        fig_wod=figure;
        plot(LFP_lpfilt.time{itrial},LFP_lpfilt.trial{itrial}(ichan,:));
        title(sprintf('WOD peak rat %s trial %d',cfg.prefix,itrial),'interpreter', 'none')
        hold on
        scatter(t_peak_wod,-v_peak_wod,'rx');
        xlim([t_peak_wod-30 t_peak_wod+30]);
        
        fig_wor=figure;
        plot(LFP_lpfilt.time{itrial},LFP_lpfilt.trial{itrial}(ichan,:));
        title(sprintf('WOR peak rat %s trial %d',cfg.prefix,itrial),'interpreter', 'none')
        hold on
        scatter(t_peak_wor,v_peak_wor,'rx');
        xlim([t_peak_wor-30 t_peak_wor+30]);
        
        %%%% modification par sofia 30/11/2021
        %             trier les plot selon la condition que la detection est bonne ou non "accepted" rejected
        
        %             detect_wod=fullfile(detectsavedir,'WoD','peak',sprintf('%s',cfg.prefix));
        %             detect_wor=fullfile(detectsavedir,'WoR','peak',sprintf('%s',cfg.prefix));
        %fixme
        %if abs(stats_all.WoD.peak_value(ichan,itrial))> WOD_threshold && abs(stats_all.WoR.peak_value(ichan,itrial))>300
        if abs(stats_all.WoD.peak_value(ichan,itrial))> WOD_threshold
            detect_wod=fullfile(detectsavedir,'WoD','peak','accepted',sprintf('%s',cfg.prefix));
            detect_wor=fullfile(detectsavedir,'WoR','peak','accepted',sprintf('%s',cfg.prefix));
        else
            detect_wod=fullfile(detectsavedir,'WoD','peak','rejected',sprintf('%s',cfg.prefix));
            detect_wor=fullfile(detectsavedir,'WoR','peak','rejected',sprintf('%s',cfg.prefix));
        end
        
        
        
        
        
        %%%% fin de la modification
        
        
        %             if ~isfolder(detectsavedir)
        %                 mkdir(detectsavedir);
        %             end
        %
        %             if ~isfolder(detect_wod)
        %                 mkdir(detect_wod);
        %             end
        %
        %             if ~isfolder(detect_wor)
        %                 mkdir(detect_wor);
        %             end
        
        
        fname_wod=fullfile(detect_wod,sprintf('%s_WoD%i_of_%i',ichan_name,itrial,size(LFP_lpfilt.trial,2)));
        fname_wor=fullfile(detect_wor,sprintf('%s_WoR%i_of_%i',ichan_name,itrial,size(LFP_lpfilt.trial,2)));
        
        
        dtx_savefigure(fig_wod,fname_wod,'png','pdf','close');
        dtx_savefigure(fig_wor,fname_wor,'png','pdf','close');
        
        clear fname_wod fname_wod detect_wod detect_wor
        
        
        
        
        %% Determine minimum and maximum slopes and extract timings and values
        
        fprintf('Launching max slope detections for %s trial %i and channel %s\n', cfg.prefix,itrial,ichan_name);
        
        
        
        %             WOD window selection
        t1= t_peak_wod-30;
        t2= t_peak_wod+30;
        t_sel= [t1 t2];
        
        %             cut data to keep only WOD
        cfgtemp=[];
        cfgtemp.latency= t_sel;
        WOD_cut= ft_selectdata(cfgtemp,LFP_lpfilt);
        clear t1 t2 t_sel
        
        %             Transform into slope
        WOD_cut_slope=WOD_cut;
        WOD_cut_slope.trial{itrial}= ft_preproc_derivative(WOD_cut.trial{itrial});
        %             smooth slope
        %             WOD_cut_slope.trial{itrial}= movmean(WOD_cut_slope.trial{itrial},100,2);
        
        %             Search for peaks in slope data
        %             Determine time window to search
        t = WOD_cut.time{itrial};
        t_1 = t > (t_peak_wod - 10);
        t_2 = t < (t_peak_wod + 10);
        t_sel = t_1 & t_2;
        
        [v_peak_wodslope, t_peak_wodslope] = findpeaks(-WOD_cut_slope.trial{itrial}(chan_idx,t_sel),t(t_sel),'NPeaks',1,'SortStr','descend','WidthReference','Halfheight');
        clear t t_1 t_2 t_sel
        
        %             save values
        
        stats_all.WoD.min_slope_time(ichan,itrial)=  t_peak_wodslope;
        stats_all.WoD.min_slope_value(ichan,itrial)=   -v_peak_wodslope;
        
        %             WOR threshold
        
        %value= 10 et 10 par defaut dans setparam
        chosen_value(1)= cfg.chosen_value{itrial}(1)
        chosen_value(2)= cfg.chosen_value{itrial}(2)
        t1= t_peak_wor- chosen_value(1);        
        t2=  t_peak_wor+chosen_value(2);
        t_sel= [t1 t2];
        
%FIXME % test sofia a tenté une boucle car out of range 
%         t = LFP_lpfilt.time{itrial}
%         t2= t(t < t_peak_wor+20);
%         t2=t2(end)
%         t1= t(t >t_peak_wor-80);
%         t1=t1(1)
        
 
        % cut data to keep only WOR
        cfgtemp=[];
        cfgtemp.latency= t_sel;
        WOR_cut = ft_selectdata(cfgtemp,LFP_lpfilt);
        clear t1 t2 t_sel
        
        %             Transform into slope
        WOR_cut_slope=WOR_cut;
        WOR_cut_slope.trial{itrial}= ft_preproc_derivative(WOR_cut.trial{itrial});
        %             smooth slope
        WOR_cut_slope.trial{itrial}= movmean(WOR_cut_slope.trial{itrial},1000,2);
        
        %             Search for peaks in slope data
        %             Determine time window to search
        t = WOR_cut.time{itrial};
        t_1 = t > (t_peak_wor - 3);
        t_2 = t < (t_peak_wor + 3);
        t_sel = t_1 & t_2;
        
        [v_peak_worslope, t_peak_worslope] = findpeaks(WOR_cut_slope.trial{itrial}(chan_idx,:),t,'NPeaks',1,'SortStr','descend','WidthReference','Halfheight');
        clear t t_1 t_2 t_sel
        
        %             save values
        %             express timings of WoR according to Vent On
        %mis commentaire par sofia
        %t_VentOn= MuseStruct{1}.markers.Vent_On.synctime(itrial)-starttrial(itrial) +offsettrial(itrial);
        
        real_timeslope_wor= t_peak_worslope - t_VentOn;
        stats_all.WoR.min_slope_time(ichan,itrial)=  real_timeslope_wor;
        stats_all.WoR.min_slope_value(ichan,itrial)=   v_peak_worslope;
        
        if  cfg.LFP.recov{itrial}==0
            stats_all.WoR.min_slope_time(ichan,itrial)=nan;
            stats_all.WoR.min_slope_value(ichan,itrial)= nan;
        end
        %             plot for visual control
        
        fig_wodslope=figure;
        plot(WOD_cut_slope.time{itrial},WOD_cut_slope.trial{itrial}(ichan,:));
        title(sprintf('WOD slope rat %s trial %d',cfg.prefix,itrial),'interpreter', 'none')
        hold on
        scatter(t_peak_wodslope,-v_peak_wodslope,'x');
        xlim([t_peak_wod-10 t_peak_wod+10]);
        
        fig_worslope=figure;
        plot(WOR_cut_slope.time{itrial},WOR_cut_slope.trial{itrial}(ichan,:));
        title(sprintf('WOR slope rat %s trial %d',cfg.prefix,itrial),'interpreter', 'none')
        hold on
        scatter(t_peak_worslope,v_peak_worslope,'x');
        xlim([t_peak_wor-10 t_peak_wor+10]);
        
        
        %%%% modification par sofia 30/11/2021
        %             trier les plot selon la condition que la detection est bonne ou non "accepted" rejected
        
        %             detectslope_wod=fullfile(detectsavedir,'WoD','slope',sprintf('%s',cfg.prefix));
        %             detectslope_wor=fullfile(detectsavedir,'WoR','slope',sprintf('%s',cfg.prefix));
        
        if abs(stats_all.WoD.peak_value(ichan,itrial))> WOD_threshold
            detectslope_wod=fullfile(detectsavedir,'WoD','slope','accepted',sprintf('%s',cfg.prefix));
            detectslope_wor=fullfile(detectsavedir,'WoR','slope','accepted',sprintf('%s',cfg.prefix));
        else
            detectslope_wod=fullfile(detectsavedir,'WoD','slope','rejected',sprintf('%s',cfg.prefix));
            detectslope_wor=fullfile(detectsavedir,'WoR','slope','rejected',sprintf('%s',cfg.prefix));
        end
        
        %%%% fin de la modification
        
        
        %             if ~isfolder(detectsavedir)
        %                 mkdir(detectsavedir);
        %             end
        %
        %             if ~isfolder(detectslope_wod)
        %                 mkdir(detectslope_wod);
        %             end
        %
        %             if ~isfolder(detectslope_wor)
        %                 mkdir(detectslope_wor);
        %             end
        
        fname_wodslope=fullfile(detectslope_wod,sprintf('%s_WoD%i_of_%i',ichan_name,itrial,size(LFP_lpfilt.trial,2)));
        fname_worslope=fullfile(detectslope_wor,sprintf('%s_WoR%i_of_%i',ichan_name,itrial,size(LFP_lpfilt.trial,2)));
        
        dtx_savefigure(fig_wodslope,fname_wodslope,'png','pdf','close');
        dtx_savefigure(fig_worslope,fname_worslope,'png','pdf','close');
        
        clear fname_wodslope fname_worslope detectslope_wod detectslope_wor
        
        % Determine threshold and crossing points
        
        %             calculate threshold 20% of max slope
        wod_thr= -0.2*v_peak_wodslope;
        wor_thr= 0.2*v_peak_worslope;
        
        %             Determine crossing point
        %             define time window
        t = WOD_cut.time{itrial};
        t_1 = t > (t_peak_wodslope - 5);
        t_2 = t < (t_peak_wodslope + 5);
        t_sel = t_1 & t_2;
        
        %             Create curve and horizontal line
        x1 = WOD_cut_slope.time{itrial}(1,t_sel);
        y1 = WOD_cut_slope.trial{itrial}(ichan,t_sel);
        x2 = x1;
        y2 = ones(size(y1)) * wod_thr;
        %             Find values of intersection of 2 curves
        [x_wodintersect, y_wodintersect] = intersections(x1, y1, x2, y2);
        
        time_start_wod= x_wodintersect(1);
        value_start_wod= y_wodintersect(1);
        
        clear t t_1 t_2 t_sel
        
        t = WOR_cut.time{itrial};
        t_1 = t > (t_peak_worslope - 10);
        t_2 = t < (t_peak_worslope + 10);
        t_sel = t_1 & t_2;
        
        %             Create curve and horizontal line
        x1 = WOR_cut_slope.time{itrial}(1,t_sel);
        y1 = WOR_cut_slope.trial{itrial}(ichan,t_sel);
        x2 = x1;
        y2 = ones(size(y1)) * wor_thr;
        %             Find values of intersection of 2 curves
        [x_worintersect, y_worintersect] = intersections(x1, y1, x2, y2);
        time_start_wor= x_worintersect(1);
        value_start_wor= y_worintersect(1);
        
        clear t t_1 t_2 t_sel
        
        %             store values
        real_time_wor=time_start_wor-t_VentOn;
        
        
        stats_all.WoD.start_time(ichan,itrial)=time_start_wod;
        stats_all.WoD.start_slope_value(ichan,itrial)=value_start_wod;
        
        if cfg.LFP.recov{itrial}==0
            real_time_wor   = nan;
            value_start_wor = nan;
        end
        
        
        stats_all.WoR.start_time(ichan,itrial)=  real_time_wor;
        stats_all.WoR.start_slope_value(ichan,itrial)=value_start_wor;
        
        %plot for visual control
        
        fig_wodthr= figure;
        plot(WOD_cut.time{itrial},WOD_cut.trial{itrial}(ichan,:));
        title(sprintf('WOD start rat %s trial %d',cfg.prefix,itrial),'interpreter', 'none')
        xline(time_start_wod);
        
        fig_worthr=figure;
        plot(WOR_cut.time{itrial},WOR_cut.trial{itrial}(ichan,:));
        title(sprintf('WOR start rat %s trial %d',cfg.prefix,itrial),'interpreter', 'none')
        xline(time_start_wor);
        
        
        
        %%%% modification par sofia 30/11/2021
        %             trier les plot selon la condition que la detection est bonne ou non "accepted" rejected
        
        %             detectstart_wod=fullfile(detectsavedir,'WoD','start',sprintf('%s',cfg.prefix));
        %             detectstart_wor=fullfile(detectsavedir,'WoR','start',sprintf('%s',cfg.prefix));
        
        if abs(stats_all.WoD.peak_value(ichan,itrial))> WOD_threshold
            detectstart_wod=fullfile(detectsavedir,'WoD','start','accepted',sprintf('%s',cfg.prefix));
            detectstart_wor=fullfile(detectsavedir,'WoR','start','accepted',sprintf('%s',cfg.prefix));
        else
            detectstart_wod=fullfile(detectsavedir,'WoD','start','rejected',sprintf('%s',cfg.prefix));
            detectstart_wor=fullfile(detectsavedir,'WoR','start','rejected',sprintf('%s',cfg.prefix));
        end
        
        %%%% fin de la modification
        
        
        %             if ~isfolder(detectstart_wod)
        %                 mkdir(detectstart_wod);
        %             end
        %
        %             if ~isfolder(detectstart_wor)
        %                 mkdir(detectstart_wor);
        %             end
        
        fname_wodstart=fullfile(detectstart_wod,sprintf('%s_WoD%i_of_%i',ichan_name,itrial,size(LFP_lpfilt.trial,2)));
        fname_worstart=fullfile(detectstart_wor,sprintf('%s_WoR%i_of_%i',ichan_name,itrial,size(LFP_lpfilt.trial,2)));
        
        dtx_savefigure(fig_wodthr,fname_wodstart,'png','pdf','close');
        dtx_savefigure(fig_worthr,fname_worstart,'png','pdf','close');
        
        
        clear x_wodintersect x_worintersect detectstart_wod detectstart_wor fname_wodstart fname_worstart
        
        %% Determine Half-width of waves
        
        fprintf('Launching half-width of WoD for %s trial %i and channel %s\n', cfg.prefix,itrial,ichan_name);
        
        
        %Calculate half amplitude of waves
        wod_amp= -v_peak_wod ;
        half_wod= wod_amp/2;
        
        %Determine time window to search
        %WOD
        
        t = WOD_cut.time{itrial};
        t_1 = t > (t_peak_wod - 10);
        t_2 = t < (t_peak_wod + 10);
        t_sel = t_1 & t_2;
        
        x1 = WOD_cut.time{itrial}(1,t_sel);
        y1 = WOD_cut.trial{itrial}(ichan,t_sel);
        x2 = x1;
        y2 = ones(size(y1)) * half_wod;
        [x_wodintersect, y_wodintersect] = intersections(x1, y1, x2, y2);
        
        if isempty(x_wodintersect) || length(x_wodintersect)==1
            WOD_halfwi = nan;
        else
            
            WOD_halfwi= x_wodintersect(2)- x_wodintersect(1);
        end
        clear x1 y1 x2 y2
        
        %Store data
        
        stats_all.WoD.half_width(ichan,itrial)=WOD_halfwi;
        
        %plot for visual control
        
        fig_wodhalf=figure;
        plot(WOD_cut.time{itrial},WOD_cut.trial{itrial}(ichan,:));
        title(sprintf('WOD halfwidth rat %s trial %d',cfg.prefix,itrial),'interpreter', 'none')
        hold on
        scatter(x_wodintersect,y_wodintersect,'rx')
        yline(half_wod);
        
        
        %%%% modification par sofia 30/11/2021
        %             trier les plot selon la condition que la detection est bonne ou non "accepted" rejected
        
        %              detecthalf_wod=fullfile(detectsavedir,'WoD','half-width',sprintf('%s',cfg.prefix));
        
        if abs(stats_all.WoD.peak_value(ichan,itrial))> WOD_threshold
            detecthalf_wod=fullfile(detectsavedir,'WoD','half-width','accepted',sprintf('%s',cfg.prefix));
        else
            detecthalf_wod=fullfile(detectsavedir,'WoD','half-width','rejected',sprintf('%s',cfg.prefix));
        end
        
        %%%% fin de la modification
        
        
        %
        %             if ~isfolder(detecthalf_wod)
        %                 mkdir(detecthalf_wod);
        %             end
        
        fname_wodhalf=fullfile(detecthalf_wod,sprintf('%s_WoD%i_of_%i',ichan_name,itrial,size(LFP_lpfilt.trial,2)));
        
        dtx_savefigure(fig_wodhalf,fname_wodhalf,'png','pdf','close');
        
        
        clear x_wodintersect  y_wodintersect
        %% Create structure with electrode depths
        
        %stats_all.Depth(ichan,itrial)=cfg.LFP.chan_depth{ichan};
        stats_all.Depth(ichan,itrial)=cfg.LFP.chan_depth(ichan); % remplacé par sofia
        
        %% Find ISO onset (new method to test)
                
        threshold_1 = 1.5; %1er seuil : 1.5 fois l'amplitude du signal isoélectrique
        threshold_2 = 2; %2eme seuil : 2 fois le premier seuil, pour détecter des pics de grande amplitude s'il y en a après la première détection
        
        cfgtemp=[];
        cfgtemp.latency=[0 wod_markertime];
        cfgtemp.channel = ichan_name;
        Activity=ft_selectdata(cfgtemp,LFP_hpfilt);
           
        %go through each window
        winsize = 1;
        iso_baseline_size = 40; %seconds : taille de la fin de l'iso qui sert de baseline
        
        winstart = 0;
        winend = winstart + winsize;
        data_detec.x = [];
        data_detec.y = [];
        while winend < Activity.time{itrial}(end)
            
            t_sel = Activity.time{itrial} > winstart & Activity.time{itrial} < winend;
            
            data_detec.x(end+1) = winstart + (winend-winstart)/2; %milieu de la fenetre
            %data_detec.y(end+1) = mean(data_abs(t_sel), 'omitnan');
            up = mean(findpeaks(Activity.trial{itrial}(t_sel), 'NPeaks', 10, 'SortStr', 'descend'));
            lo = mean(findpeaks(-Activity.trial{itrial}(t_sel), 'NPeaks', 10, 'SortStr', 'descend'));
            data_detec.y(end+1) = up+lo;
            
            
            winstart = winstart + winsize;
            winend = winend + winsize;
        end
        
        t_sel = data_detec.x > data_detec.x(end) - iso_baseline_size;
        threshold = threshold_1 * mean(data_detec.y(t_sel)); 
        
        %find first threshold crossing
        ISO_idx = find(data_detec.y < threshold, 1, 'first');
        
        %search on the rest of the data if there is any big peak and if so,
        %move the detection after it
        new_detec_idx = [];
        for i = ISO_idx:size(data_detec.y, 2)
            if data_detec.y(i) > threshold_2 * threshold
                new_detec_idx = i;%moment où croise le 2ème seuil 
                %ISO_idx = i + 1;
            end
        end
        if ~isempty(new_detec_idx) %si le signal dépasse 2 * la baseline après la première détection, refaire une détection à partir de ce nouvel endroit
            ISO_idx = find(data_detec.x > data_detec.x(new_detec_idx) & data_detec.y < threshold, 1, 'first');
        end
              
        ISO_time = data_detec.x(ISO_idx);
        %added by sofia
        if isempty(ISO_time)
            ISO_time= nan;
        end
        %end added by sofia
        
        %             plot(Activity.time{itrial},amp_signal);
        %             yline(thr);
        
        stats_all.ISO(ichan,itrial)=ISO_time;
        
       
        
        fig=figure; hold on;
        plot(Activity.time{itrial},Activity.trial{itrial}, 'k');
        plot(data_detec.x, data_detec.y, 'r', 'linewidth', 2);
        yline(threshold, 'b');
        if ~isnan(ISO_time)
            xline(ISO_time, 'r');
        end
        title(sprintf('ISO detection rat %s trial %d',cfg.prefix,itrial),'interpreter', 'none')
        
        if abs(stats_all.WoD.peak_value(ichan,itrial))> WOD_threshold
            detect_iso=fullfile(detectsavedir,'WoD','ISO','accepted',sprintf('%s',cfg.prefix));
        else
            detect_iso=fullfile(detectsavedir,'WoD','ISO','rejected',sprintf('%s',cfg.prefix));
        end
        
        fname_iso=fullfile(detect_iso,sprintf('%s_WoD%i_of_%i',ichan_name,itrial,size(LFP_lpfilt.trial,2)));
        
        dtx_savefigure(fig,fname_iso,'png','pdf','close');
        
        clear ISO_time thr_cross idx_signal

        %% Determine Iso delay with signal amplitude
        %travailler sur un LFP filtré >1Hz pour retirer le bruit lent
        
        
        %select baseline
        %FIXME : pas besoin de retrouver v_off car le t0 du LFP
        %correspond déjà à v_off
%         voff_marker = MuseStruct{1}.markers.Vent_Off.synctime(itrial);% == 0
        
        %             %correction apporté par sofia
        %             temp_voff_marker = []
        %             if itrial>1
        %             temp_voff_marker =((LFP_sample.sampleinfo/MuseStruct{1, 1}.markers.Stopping_Recording.synctime(itrial))* voff_marker)/LFP.fsample
        %             voff_marker = temp_voff_marker
        %             end
        %             %fin de la correction
%         voff_time= voff_marker - starttrial(itrial) + offsettrial(itrial);
        
        
%         %select baseline
%         %             t_1=voff_time-60;
%         %             t_2=voff_time;
%         t_sel=[-60 0];
%         
%         cfgtemp=[];
%         cfgtemp.latency=t_sel;
%         Baseline=ft_selectdata(cfgtemp,LFP_hpfilt);
%         clear t_1 t_2 t_sel
%         
%         %selectuntil WoD
%         t_1=0;
%         t_2=wod_markertime;
%         t_sel= [t_1 t_2];
%         
%         cfgtemp=[];
%         cfgtemp.latency=t_sel;
%         Activity=ft_selectdata(cfgtemp,LFP_hpfilt);
%         clear t_1 t_2 t_sel
%         
%         % voir si retirer ci dessous :
%         [up,lo]=envelope(Baseline.trial{itrial}(ichan,:),1000,'peak');
%         baseline_amp=mean(up-lo);
%         clear up lo
%         
%         [up,lo]=envelope(Activity.trial{itrial}(ichan,:),1000,'peak');
%         amp_signal=up-lo; % cumule
%         clear up lo
%         thr=0.2*baseline_amp;
%         
%         %             % Test Paul : ne pas utiliser
%         %             baseline_amp = mean(abs(Baseline.trial{itrial}(ichan,:)));
%         %             thr=baseline_amp;
%         %             amp_signal   = abs(Activity.trial{itrial}(ichan,:));
%         
%         
%         idx_signal= amp_signal>thr(1);
%         diff_signal=diff(idx_signal);
%         thr_cross=find(diff_signal); %index du(des) samples qui dépssent le seuil
%         
%         if  ~isempty(thr_cross)
%         thr_cross = thr_cross(end);
%         else % si ne detecte pas l'iso augmenter le seuil
%         clear thr idx_signal diff_signal thr_cross
%         thr=0.5*baseline_amp;
%         idx_signal= amp_signal>thr(1);
%         diff_signal=diff(idx_signal);
%         thr_cross=find(diff_signal); 
%         thr_cross = thr_cross(1);%thr_cross(end)
%         end
%         
%         ISO_time = Activity.time{itrial}(thr_cross); %temps en secondes après V_off
%         
%         %             if isempty(thr_cross)
%         %                 thr_iso=0.5*baseline_amp;
%         %                 idx_signal= Activity.trial{itrial}>thr_iso;
%         %                 diff_signal=diff(idx_signal);
%         %                 thr_cross=find(diff_signal);
%         %             end
%         %
%         %             if length(thr_cross)>2
%         %                 if length(thr_cross)>1000%added by sofia
%         %                     error('ça ne devrait pas arriver : vérifier ce qui se passe pour ce rat : %s trial %d chan %d', cfg.prefix, itrial, ichan);
%         %                     ISO_time=nan;%%ici
%         %                 else%ici
%         %                     %idx_ISO=diff(thr_cross)>=4;
%         %                     ISO_time= Activity.time{itrial}(thr_cross(end));
%         %
%         % %                     idx=  ISO_time> 0.1*max(Activity.time{itrial});
%         % %                     ISO_time=ISO_time(find(idx,1,'last'));
%         %                 end %ici
%         %             else
%         %                 idx_ISO=thr_cross;
%         %                 ISO_time=(Activity.time{itrial}(idx_ISO));
%         %             end
%         %
%         %             if  ~isnan(ISO_time)%added by sofia
%         %             idx=  ISO_time> 0.1*max(Activity.time{itrial});
%         %             ISO_time=ISO_time(find(idx,1,'first'));
%         %
%         fig_iso=figure; hold on;
%         plot(Activity.time{itrial},Activity.trial{itrial}(ichan, :), 'k');
%         plot(Activity.time{itrial},amp_signal, 'r');
%         yline(thr, 'b');
%         yline(baseline_amp, 'g');
%         scatter(ISO_time, thr, 'rx');
%         xline(ISO_time, 'r');
%         title(sprintf('ISO detection rat %s trial %d',cfg.prefix,itrial),'interpreter', 'none')
        
        
        %%%% modification par sofia 30/11/2021
        %             trier les plot selon la condition que la detection est bonne ou non "accepted" rejected
        
        %             detect_iso=fullfile(detectsavedir,'ISO',sprintf('%s',cfg.prefix));
        
%         if abs(stats_all.WoD.peak_value(ichan,itrial))> WOD_threshold
%             detect_iso=fullfile(detectsavedir,'WoD','ISO','accepted',sprintf('%s',cfg.prefix));
%         else
%             detect_iso=fullfile(detectsavedir,'WoD','ISO','rejected',sprintf('%s',cfg.prefix));
%         end
%         
%         %             if ~isfolder(detect_iso)
%         %                 mkdir(detect_iso);
%         %             end
%         
%         fname_iso=fullfile(detect_iso,sprintf('%s_WoD%i_of_%i',ichan_name,itrial,size(LFP_lpfilt.trial,2)));
%         
%         dtx_savefigure(fig_iso,fname_iso,'png','pdf','close');
        
        %%%% fin de la modification
        
    %end%added by sofia
    
%     %added by sofia
%     if isempty(ISO_time)
%         ISO_time= nan;
%     end
%     %end added by sofia
%     
%     %             plot(Activity.time{itrial},amp_signal);
%     %             yline(thr);
%     
%     stats_all.ISO(ichan,itrial)=ISO_time;
%     
%     clear ISO_time thr_cross idx_signal
    
    %% condtion pour ignorer un channel si celui_ci est trop bruité ou qu'il ne contient pas de WOD
    
    if  abs(stats_all.WoD.peak_value(ichan,itrial)) < WOD_threshold
        
        stats_all.WoR.peak_time(ichan,itrial)          = nan;
        stats_all.WoR.peak_value(ichan,itrial)         = nan;
        stats_all.WoD.peak_time(ichan,itrial)          = nan;
        stats_all.WoD.peak_value(ichan,itrial)         = nan;
        
        stats_all.WoD.min_slope_time(ichan,itrial)     = nan;
        stats_all.WoD.min_slope_value(ichan,itrial)    = nan;
        stats_all.WoR.min_slope_time(ichan,itrial)     = nan;
        stats_all.WoR.min_slope_value(ichan,itrial)    = nan;
        
        stats_all.WoD.start_time(ichan,itrial)         = nan;
        stats_all.WoD.start_slope_value(ichan,itrial)  = nan;
        stats_all.WoR.start_time(ichan,itrial)         = nan;
        stats_all.WoR.start_slope_value(ichan,itrial)  = nan;
        
        stats_all.WoD.half_width(ichan,itrial)         = nan;
        stats_all.WoR.half_width(ichan,itrial)         = nan;
        
        stats_all.Depth(ichan,itrial)                  = nan;
        stats_all.ISO(ichan,itrial)                    = nan;
        
    end
end

[origin_time, idx_origin] = min(stats_all.WoD.peak_time(:,itrial));
[origin_time_wor, idx_origin_wor] = min(stats_all.WoR.peak_time(:,itrial));


if sum(origin_time == stats_all.WoD.peak_time(:,itrial)) > 1
    idx_origin = find(origin_time == stats_all.WoD.peak_time(:,itrial), 1, 'last');
end
if sum(origin_time_wor == stats_all.WoR.peak_time(:,itrial)) > 1
    idx_origin_wor = find(origin_time_wor == stats_all.WoR.peak_time(:,itrial), 1, 'last');
end
origin_depth = stats_all.Depth(idx_origin,itrial);
origin_depth_wor= stats_all.Depth(idx_origin_wor,itrial);
stats_all.wod_origin_depth(itrial)=origin_depth;
stats_all.wod_origin_time(itrial)=origin_time;
stats_all.wor_origin_depth(itrial)=origin_depth_wor;
stats_all.wor_origin_time(itrial)=origin_time_wor;



%separate wod protocols by origin depth
limits = [700, 1300, 1700];
depth_lim{irat}(itrial) = 0;
for i = limits
    if origin_depth > i
        depth_lim{irat}(itrial) = depth_lim{irat}(itrial) + 1;
    end
end
if depth_lim{irat}(itrial)==1
    stats_all.oridepthclass(itrial)=1000;
elseif depth_lim{irat}(itrial)==2
    stats_all.oridepthclass(itrial)=1400;
elseif depth_lim{irat}(itrial)==0 || depth_lim{irat}(itrial)==3
    stats_all.oridepthclass(itrial)=0;
    
end


%% Calculate propagation speed
tdiff = stats_all.WoD.peak_time(:, itrial) - stats_all.WoD.peak_time(idx_origin, itrial);
ddiff =stats_all.Depth(:, itrial) - stats_all.Depth(idx_origin, itrial);

speed = ddiff ./ tdiff;
speed_up = mean(speed(1:idx_origin-1));
speed_dn = mean(speed(idx_origin+1:end));

stats_all.WoD.speed_up(itrial)=speed_up;
stats_all.WoD.speed_dn(itrial)=speed_dn;

%% Calculate anoxia duration

%until wod
voff_marker = MuseStruct{1}.markers.Vent_Off.synctime(itrial);
voff_time= voff_marker - starttrial(itrial) + offsettrial(itrial);

wod_time= wod_marker - starttrial(itrial) + offsettrial(itrial);

stats_all.before_wod(itrial)=wod_time-voff_time;

%total anoxia
von_marker = MuseStruct{1}.markers.Vent_On.synctime(itrial);
von_time= von_marker - starttrial(itrial) + offsettrial(itrial);

stats_all.anoxia(itrial)=von_time-voff_time;


end %itrial

%end %irat


%% Save structures

fname_out= fullfile(cfg.datasavedir,'Detection','stat_eachrat');
if ~isfolder(fname_out)
    mkdir(fname_out)
end
fname_out_each=fullfile(fname_out,sprintf('stat_eachrat%d.mat',irat));
save(fname_out_each,'stats_all')

statsfile=dir(fname_out);
statsfile=statsfile(3:end);
for ifiles=1:size(statsfile,1)
    rat_num=str2num(cell2mat(regexp(statsfile(ifiles).name,'\d+','match'))); %give the real rat number of the file within the file list
    wavedetection_allrat{rat_num}=load(fullfile(statsfile(ifiles).folder,statsfile(ifiles).name));
    wod_wavedetection_allrat{rat_num}=wavedetection_allrat{1,rat_num}.stats_all;
end

fname_dir= fullfile(cfg.datasavedir,'Detection','wod_wavedetection_allrat.mat');
save(fname_dir,'wod_wavedetection_allrat')


end




