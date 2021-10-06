function wod_project_sofia(irat)

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

anovastatpath=fullfile(config16{4}.statsavedir,'Waves_detection','anova');
tablesavedir=fullfile(fileparts(config16{4}.datasavedir),'tables','matlab');


if ~isfolder(anovastatpath)
    mkdir(anovastatpath);
end

if ~isfolder(tablesavedir)
    mkdir(tablesavedir)
end


ipart = 1;
% irat=nrat;

%% create muse marker file
% A décommenter pour les nouvelles manips
% 
% for irat = nrat%1:size(config, 2)
%    
%     wod_create_Muse_Marker_file(config{irat});
% 
% end


%% read data

MuseStruct  = readMuseMarkers(config{irat}, true);
LFP         = readLFP(config{irat}, MuseStruct, false);

% cfgtemp         = [];
% cfgtemp.channel = 'E28LFP';%{'E28LFP', 'E29LFP'};
% data            = ft_selectdata(cfgtemp, LFP{ipart}.WOD_morpho);
% 
% figure;
% plot(data.time{1}, data.trial{1});

%% analysis

% Compute TFR for each rat
wod_tfr_compute(config{irat}, MuseStruct, LFP{ipart}.WoD_short);

%Plot TFR data for each rat
wod_tfr_plotrat(config{irat});



end


%detect iso time
freq_data= wod_iso_domfreq(config,[],true);   


% end


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
%     
% end %irat
% 
% % if rat_list==0
% %     %% Analysis for all rats
% %
% %     %Detect waves, extract timings and values
% %     wod_wavedetection(cfg);
% 



