function wod_project_sofia(irat)

%avant de lancer le projet pour plusieurs rat verifier les range de
%detection notament pour la WOR qu'il faudrait resserer : indiqué par un %FIXME


%% set parameters

try %en local
    scriptpath = matlab.desktop.editor.getActiveFilename;
catch %cluster
    scriptpath = mfilename('fullpath');
end

epicodepath = [fileparts(fileparts(fileparts(scriptpath))), filesep];

addpath (genpath([epicodepath,'shared']))
addpath (genpath([epicodepath,'external']))
addpath (genpath([epicodepath,'templates']))
addpath (genpath([epicodepath,'projects', filesep, 'wod']))
addpath (genpath([epicodepath,'projects', filesep, 'dtx']))
addpath (genpath([epicodepath,'development']))

if ispc
    addpath \\lexport\iss01.charpier\analyses\wod\fieldtrip-20200607
elseif isunix
    addpath /network/lustre/iss01/charpier/analyses/wod/fieldtrip-20200607
end

ft_defaults

config = wod_setparam_polyprobe;

ipart = 1;

%% read data per rat

% create muse marker file
if irat>0
    
    config{irat}.to_ignore = ft_getopt(config{irat}, 'to_ignore', false);
    if config{irat}.to_ignore %irat==4 , 3
        return
    end
    
    
    wod_create_Muse_Marker_file(config{irat}); %A décommenter pour les nouvelles manips
    MuseStruct        = readMuseMarkers(config{irat}, true);
    LFP               = readLFP(config{irat}, MuseStruct, true);
    
    %verifier qu'il y a bien autant de trials que de marqueurs Vent_Off
    startmarker = config{irat}.muse.startmarker.WoD_short;
    MuseStruct_concat = concatenateMuseMarkers(config{irat}, MuseStruct, true);
    if size(LFP{1}.WoD_short.trial,2) ~= size(MuseStruct_concat{1}.markers.(startmarker).synctime,2)
        error('Not the same number of trials that of marker start for %s. \nPossible solutions : \n- One marker is missing on Muse \n- Check that begin/end of each trial is not before start of file or after end of file', config{irat}.prefix(1:end-1));
    end
    
    %select channels of interest
    cfgtemp           = [];
    cfgtemp.channel   = {'all'};%{'E28LFP', 'E29LFP'};
    [data]            = ft_selectdata(cfgtemp, LFP{ipart}.WoD_short);
    
    
    % %%%%%%%
    % f=figure(1);
    % hold on
    %     for ichan= 1:size(data.label,1)
    %     plot(data.time{1}, data.trial{1}(ichan,:));
    %     end
    %  hold off
    %  close(f)
    %
    %  %LFP par structure
    %  structure=str2num(cell2mat((categories(categorical(cell2mat(config{irat}.LFP.catego))))))
    % for istruct =1:size(structure,1)
    %     idx      = []
    %     chantemp = []
    %     cat=structure(istruct)
    %     idx=find(cell2mat(config{irat}.LFP.catego)==cat);
    %         for nidx=1:length(idx)
    %         ichantemp=char(config{irat}.LFP.channel(idx(nidx)));
    %         chantemp= [chantemp, {ichantemp}];
    %         end
    %
    % cfgtemp         = [];
    % cfgtemp.channel =  chantemp;
    % [data]            = ft_selectdata(cfgtemp, LFP{ipart}.WoD_short);
    %         f=figure(istruct);
    %         hold on
    %             for ichan= 1:size(data.label,1)
    %             plot(data.time{1}, data.trial{1}(ichan,:));
    %             end
    %         hold off
    %         close(f)
    %
    %
    % end
    %
    %  %TFR par structure
    % %%%%%%%%%
    
    
    
    %% analysis per rat
    
    % Compute TFR for each rat
    wod_tfr_compute(config{irat}, MuseStruct, data);
    
    %stat
    wod_wavedetection_sofia(irat, config{irat}, data);
    
    
    %% analisis pour loop over all rats
    
elseif irat==0
    
    [config,rongeur,per_rat]=wod_per_structure(config);
    
    analysis_list = ["timefreq_baseline", "timefreq_wod_blcorrected"];
    [power_detection]=wod_powerspectrum_detection(config, rongeur, per_rat, analysis_list, true); %plot les power spectrum des chan initiation
    %pb dans le boucle itrial donc j'ai remplacé par 1
    %j'ai tenté de modifier readLFP car ne trouve pas les bon trial....
    
    proba_density(config,rongeur,per_rat)
    
    wod_correlation(config,rongeur,per_rat)
    
    %le dernier argument 'false' permet de charger les données qui ont été
    %calculées précédemment :
    [power_detection]=wod_powerspectrum_detection(config,[],[], [], false);

    wod_powerspectrum_stats(config, power_detection)
    %         analysis_names = {'timefreq_wod', 'timefreq_wod_timenorm', 'timefreq_baseline','timefreq_wod_blcorrected', 'timefreq_wod_timenorm_blcorrected', 'timefreq_baseline_blcorrected','log_timefreq_wod', 'log_timefreq_wod_timenorm', 'log_timefreq_baseline','log_timefreq_wod_blcorrected', 'log_timefreq_wod_timenorm_blcorrected','log_timefreq_baseline_blcorrected'};
    %
    %         S = {'CC','HPC', 'NC', 'PTA', 'S1', 'TH'}
    %         savedir = '\\lexport\iss01.charpier\analyses\wod\Sofia\images'
    %         %       analysis_names = {'timefreq_wod', 'timefreq_wod_timenorm', 'log_timefreq_baseline','log_timefreq_baseline_blcorrected', 'timefreq_wod_timenorm_blcorrected', 'timefreq_baseline_blcorrected','log_timefreq_wod', 'log_timefreq_wod_timenorm', 'log_timefreq_baseline','log_timefreq_wod_blcorrected', 'log_timefreq_wod_timenorm_blcorrected','log_timefreq_baseline_blcorrected'};
    %
    %         for idata = 1%:length(analysis_names);
    %         for nrat=1:length(config)
    %             if nrat==4
    %                 continue
    %             end
    %
    %             for icat=1:size(config{nrat}.classe.structure,2)
    %                 for itrial=1:size(rongeur(nrat).catego(icat).wod_time_ini,2)  %size(config{nrat}.directorylist{1},2) marche pas
    %                     interstruct_chan={};
    %
    %                     if config{nrat}.classe.structure(icat).trial ==1
    %
    %                         rongeur(nrat).catego(icat).wod_time_ini;
    %                         rongeur(nrat).catego(icat).chan_depth;
    %
    %                         cfgtemp = [];
    %                         MuseStruct        = readMuseMarkers(config{nrat}, true);
    %                         LFP               = readLFP(config{nrat}, MuseStruct, false);
    %
    %                         %verifier si le channel existe bien
    %                         cfgtemp.channel   = rongeur(nrat).catego(icat).chan_ini{1,1};
    %                         if isnan(cell2mat(cfgtemp.channel))
    %                             continue
    %                         end
    %
    %                         %les analyses possibles
    %
    %                         data_temp = load(fullfile(config{nrat}.datasavedir,[config{nrat}.prefix,analysis_names{idata},'.mat']))
    %                         for ifield = string(fieldnames(data_temp.timefreq_wod{itrial}))'
    %                             if isempty(data_temp.timefreq_wod{itrial}.(ifield))
    %                                 data_temp.timefreq_wod{itrial} = rmfield(data_temp.timefreq_wod{itrial}, ifield);
    %                             end
    %                         end
    %                         temp =  struct2cell(data_temp.timefreq_wod{itrial});
    %                         temp_ini=data_temp.timefreq_wod{itrial}.(cell2mat(cfgtemp.channel))
    %                         %powerspectrum
    %                         cfgtemp.parameter  = 'powspctrm';
    %                         %tfr_all = ft_appendfreq(cfgtemp, temp{:})
    %                         tfr_ini = ft_appendfreq(cfgtemp, temp_ini)
    %
    %                         cfgtemp = [];
    %                         cfgtemp.frequency = [0 50];%Hz
    %                         cfgtemp.latency = [0 30];  %s
    %                         cfgtemp.channel = rongeur(nrat).catego(icat).chan_ini{1,1};
    %                         data_ini = permute(tfr_ini.powspctrm, [4 3 2 1 ]);
    %                         f=figure()
    %
    %                         %plot
    %                         plot(tfr_ini.time, data_ini,'LineWidth',2);
    %                         legend()
    %                         xlabel ('time (s)','Fontsize',15)
    %                         ylabel ('power ()','Fontsize',15)
    %                         title({['Initial WOD electrode powerspectrum'], ['in structure ' S{icat} ' at depth ' num2str( rongeur(nrat).catego(icat).chan_depth) '\mum'], [config{nrat}.prefix ' trial ' num2str(config{nrat}.trial(itrial))]},'Fontsize',12)
    %                         if ~exist(fullfile(savedir,'powspctrm',cell2mat(analysis_names(idata)),config{nrat}.prefix))
    %                             mkdir(fullfile(savedir,'powspctrm',cell2mat(analysis_names(idata)),config{nrat}.prefix))
    %                         end
    %                         saveas(f,fullfile(savedir,'powspctrm',cell2mat(analysis_names(idata)),config{nrat}.prefix, ['Initial_powerspectrum' S{icat} '.jpg']))
    %                         close(f)
    %
    % %                                          interstruct_chan{icat} = cfgtemp.channel;
    % %                                          [data]            = ft_selectdata(cfgtemp, LFP{ipart}.WoD_short);
    % %                                          %plot
    % %                                          f=figure(icat);
    % %                                          plot(data.time{1}, data.trial{1});
    % %                                          %close(f)
    %                         end
    %
    %                     end
    %                 end
    %             end
    % %                     cfgtemp=[];
    % %                     cfgtemp.channel=interstruct_chan;
    % %                     [c, v, n] = ft_connectivity_corr(cfgtemp.channel,'hasjack', 0)
    %         end
    %
    
    %selectionne le channel d'interet et fait toi plez
    %                 if ~ismember(mat2str(cell2mat(rongeur(nrat).catego(istruct).chan_ini{1,1})),'NaN')
    %                  MuseStruct        = readMuseMarkers(config{nrat}, true);
    %                  LFP               = readLFP(config{nrat}, MuseStruct, false)
    %                  cfgtemp           = [];
    %                  cfgtemp.channel   = rongeur(nrat).catego(istruct).chan_ini{1,1};
    %                  cfgtemp.latency   = [0 10];
    %                  data              = ft_selectdata(cfgtemp, LFP{1}.WoD_short);
    %                  cfgtemp           = [];
    %                  analysis_names    = {'timefreq_wod', 'timefreq_wod_timenorm', 'timefreq_baseline','timefreq_wod_blcorrected', 'timefreq_wod_timenorm_blcorrected', 'timefreq_baseline_blcorrected','log_timefreq_wod', 'log_timefreq_wod_timenorm', 'log_timefreq_baseline','log_timefreq_wod_blcorrected', 'log_timefreq_wod_timenorm_blcorrected','log_timefreq_baseline_blcorrected'};
    %                  for idata = 1%:length(analysis_names)
    %                      data_temp          = load(fullfile(config{nrat}.datasavedir,[config{nrat}.prefix,analysis_names{idata},'.mat']));
    %                      temp               = struct2cell(data_temp.timefreq_wod{itrial});
    %                      cfgtemp.channel    = rongeur(nrat).catego(istruct).chan_ini{1,1};
    %                      cfgtemp.parameter  = 'powspctrm';
    %                      tfr_all            = ft_appendfreq(cfgtemp, temp{:});
    %                  end
    %
    %                  % average over chan
    %                  data_perchan= tfr_all.powspctrm;
    %                  data_perchan = permute( data_perchan, [2 3 1]);
    %                  plot(tfr_all.time, data_perchan);
    %
    %
    %                 end
    
    
    
    
    
    
    
    
end
end













%% OLD

%
% ipart=1
% for irat=1:length(rat);
%     for itrial=1:size(config{irat}.directorylist{1},2)
%
%         path    = dir(rat(irat).folder);
%         path    = path(3:end);
%         trials  = dir(fullfile(path(irat).folder,path(irat).name));
%         trials=trials(3:end)
%         Todata=dir(fullfile(trials(irat).folder,trials(irat).name))
%         Todata=Todata(3:end)
%         MuseStruc = load('Events_nev.mat')
%
%         config{irat}.directorylist{1}{1}={}
%         config{irat}.directorylist{1}{1} = trials(irat).name
%             for ichan=3:length(Todata)
%
%             LFP = readLFP(config{irat}, MuseStruct,false);
%             end
%
%     end
% end
%
%
% %% analysis by rat
%
%
%  % find concatenated LFP (see wod_concatenateLFP.m)
%     %         [~,dir_name]                       = fullfile(cfg{irat}.concatdata_path,cfg{irat}.prefix);
%     config{irat}.rawdir               = fullfile(config{irat}.concatdata_path);
%     config{irat}.directorylist{ipart} = {config{irat}.prefix};
%     %read Muse markers
%     MuseStruct = readMuseMarkers(config{irat}, true);
%     %read LFP, append electrodes, and cut into trials according to Muse Markers
%     LFP = readLFP(config{irat}, MuseStruct,false);
%
%     % add all metadata in LFP.trialinfo
%     Table_experience = readtable('\\lexport\iss01.charpier\raw\rat-wod\2probes\Data_table\Tableau_Experience_per_trial.xlsx');
%     Table_exp = table2cell(Table_experience);
%
%    path = '\\lexport\iss01.charpier\raw\rat-wod\2probes\Data_table\maxi_tab'
%    liste = dir(path)
%    liste = liste(3:end)
%    file  =  liste(irat).name
%    trial_maxitable = load(fullfile(path,file))
%    trial_maxitable =  trial_maxitable.maxidata
%    LFP{2,1}.trialinfo=struct
%
%    for itrial = 1 : size(LFP{1}.WoD_short.trial, 2)
%         for ichan=1:size(cellstr(LFP{1}.WoD_short.label),1)
%         LFP{2,1}.trialinfo(itrial).recov(ichan)    = cell2mat( Table_exp(itrial,4));
%         LFP{2,1}.trialinfo(itrial).X(ichan)        = trial_maxitable(ichan,10);
%         LFP{2,1}.trialinfo(itrial).Y(ichan)        = trial_maxitable(ichan,11);
%         LFP{2,1}.trialinfo(itrial).Z(ichan)        = trial_maxitable(ichan,12);
%         LFP{2,1}.trialinfo(itrial).STRC(ichan)     = trial_maxitable(ichan,14);
%         LFP{2,1}.trialinfo(itrial).substruc(ichan) = trial_maxitable(ichan,15);
%         LFP{2,1}.trialinfo(itrial).SilentD(ichan)  = trial_maxitable(ichan,21);
%         LFP{2,1}.trialinfo(itrial).HRunder(ichan)  = trial_maxitable(ichan,19);
%         LFP{2,1}.trialinfo(itrial).HRrecov(ichan)  = trial_maxitable(ichan,20);
%         LFP{2,1}.trialinfo(itrial).electrode(ichan)= trial_maxitable(ichan,13);
%
%         end
%     end
%
%
%     TFR = TFRtrials(config{irat}, MuseStruct, true);
%
% %     plot TFR look at fieldtrip website: ft_singleplotTFR
%
%     %choisir wod short ou wod long pour les analyses
%     %pour séparer pour les analyses temps fréquences qui doivent être
%     %faites sur wod_long et pas wod_short.
%
%     %
%     %         LFP = LFP(1).WoD_short; ?????????????
%     LFP = LFP{1}.WoD_long;
%     %
%
%     %vérifier qu'il y a bien autant de trials que de marqueurs Vent_Off
%     startmarker = config{irat}.muse.startmarker.(config{irat}.LFP.name{1});
%     if size(LFP.trial,2) ~= size(MuseStruct{1}{1}.markers.(startmarker).synctime,2)
%         error('Not the same number of trials that of marker start for %s. \nCheck that begin/end of each trial is not before start of file or after end of file', config{irat}.prefix(1:end-1));
%     end
%     % Compute TFR for each rat
%     %wod_tfr_compute(cfg{irat}, MuseStruct,LFP{1});??????????????
%
%     %Plot TFR data for each rat
%     %wod_tfr_plotrat(cfg{irat});
%
%
%
%     end



