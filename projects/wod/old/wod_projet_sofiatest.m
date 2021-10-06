function wod_projet_sofia(rat_list)

% Ce script projet sert à calculer les données de chaque rat
% 'irat' est en input car il prendra la valeur du array slurm sur le cluster
% Rassembler les données de tous les rats pour faire des moyennes ou stats
% sur tous les rats : dans un autre script (si possible organisé comme
% celui-ci)

% a faire pour réorganiser les fonctions WOD :
% - mettre en input : cfg, MuseStruct, LFP, et toutes les structures qui sont calculées
% dans une autre fonction
% - retirer les addpath, les cfg = wod_setparams, les boucles irat
% - remplacer tous les cfg{irat} par cfg
% - si possible, sauvegarder les données en fin de script, et créer la
% possibilité de les charger avec l'argument 'force'

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

cfg = wod_setparam_polyprob;

% cfg = wod_setparams_32chan;
% cfgorig = cfg;

ipart= 1;


%% analysis by rat

for irat = 1
    
    %orgarnser les données en fonction des structures
    for ichan=1:size(cfg{irat}.LFP.rename,2) 
%         %moove file
%         ratdir      = fullfile(cfg{irat}.concatdata_path,cfg{irat}.prefix);
%         oldfiledir  = fullfile(cfg{irat}.concatdata_path,cfg{irat}.prefix, [cfg{irat}.prefix cfg{irat}.LFP.rename{ichan} '.ncs'] );
%         newfiledir  = fullfile(cfg{irat}.concatdata_path, cfg{irat}.prefix, cfg{irat}.LFP.structure{ichan});%,[cfg{irat}.prefix cfg{irat}.LFP.rename{ichan} '.ncs'])
%         mkdir(ratdir,char(cfg{irat}.LFP.structure{ichan}));
%         mkdir(fullfile(ratdir,char(cfg{irat}.LFP.structure{ichan})));
%         file        =[cfg{irat}.prefix cfg{irat}.LFP.rename{ichan} '.ncs'];
%         A = oldfiledir;
%         B = char(newfiledir);
%         f_raw=dir( oldfiledir);
%         movefile(fullfile(f_raw.folder,f_raw.name),B)


       
    end
    
    % find concatenated LFP (see wod_concatenateLFP.m)
    %         [~,dir_name]                       = fullfile(cfg{irat}.concatdata_path,cfg{irat}.prefix);
    cfg{irat}.rawdir               = fullfile(cfg{irat}.concatdata_path);
    cfg{irat}.directorylist{ipart} = {cfg{irat}.prefix};
    
    %read Muse markers
    MuseStruct = readMuseMarkers(cfg{irat}, true);
    
    %read LFP, append electrodes, and cut into trials according to Muse Markers
    LFP = readLFP(cfg{irat}, MuseStruct,true);
    
    % add all metadata in LFP.trialinfo
    for itrial = 1 : size(LFP.trial, 2)
        LFP.trialinfo.silenttime(itrial) = silenttime(itrial);
        etc...
    end
    

    TFR = TFRtrials(cfg{irat}, MuseStruct, true);

%     plot TFR look at fieldtrip website: ft_singleplotTFR
    
    %choisir wod short ou wod long pour les analyses
    %pour séparer pour les analyses temps fréquences qui doivent être
    %faites sur wod_long et pas wod_short.
    
    %
    %         LFP = LFP(1).WoD_short; ?????????????
    LFP = LFP{1}.WoD_long;
    %
    
    %vérifier qu'il y a bien autant de trials que de marqueurs Vent_Off
    startmarker = cfg{irat}.muse.startmarker.(cfg{irat}.LFP.name{1});
    if size(LFP.trial,2) ~= size(MuseStruct{1}{1}.markers.(startmarker).synctime,2)
        error('Not the same number of trials that of marker start for %s. \nCheck that begin/end of each trial is not before start of file or after end of file', cfg{irat}.prefix(1:end-1));
    end
    % Compute TFR for each rat
    %wod_tfr_compute(cfg{irat}, MuseStruct,LFP{1});??????????????
    
    %Plot TFR data for each rat
    %wod_tfr_plotrat(cfg{irat});
    
    
    
    
end %irat

% if rat_list==0
%     %% Analysis for all rats
%
%     %Detect waves, extract timings and values
%     wod_wavedetection(cfg);




end