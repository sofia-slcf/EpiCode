function [config] = wod_setparams_32chan

disp('setting parameters');

if ismac
    error('Platform not supported')
elseif isunix
    rootpath_analysis   = '/network/lustre/iss01/charpier/analyses/wod/Sofia';
    rootpath_data       = '/network/lustre/iss01/charpier/raw/rat-wod/Cortex/32_ch/Extra_Neuralynx';
    rootpath_concatdata = '/network/lustre/iss01/charpier/raw/rat-wod/Cortex/32_ch/concatenated_LFP';
    os                  = 'unix';
elseif ispc
    rootpath_analysis	= '\\lexport\iss01.charpier\analyses\wod\Sofia';
    rootpath_data       = '\\lexport\iss01.charpier\raw\rat-wod\Cortex\32_ch\Extra_Neuralynx';
    rootpath_concatdata = '\\lexport\iss01.charpier\raw\rat-wod\Cortex\32_ch\concatenated_LFP';
    os                  = 'windows';
else
    error('Platform not supported')
end 

datasavedir  =  fullfile(rootpath_analysis,'data');
imagesavedir =  fullfile(rootpath_analysis,'images');
imagesavedir_data= {fullfile(imagesavedir,'TFR'), fullfile(imagesavedir,'band_depth'),fullfile(imagesavedir,'band_cx')};
script_path  = mfilename('fullpath');
script_path  = fileparts(script_path);

%% config common for all patients
configcommon.concatdata_path = rootpath_concatdata;
configcommon.muse.templatemarker   = fullfile(datasavedir,'TemplateEventsWOD.mrk');%find the template file to create muse marker file

configcommon.name                  = {'WoD'};
configcommon.LFP.allchannel        = {'E32','E31','E30', 'E29', 'E28','E27', 'E26', 'E25', 'E24','E23', 'E22', 'E21', 'E20', 'E19', 'E18','E17','E16', 'E15', 'E14', 'E13', 'E12', 'E11','E10', 'E09', 'E08', 'E07', 'E06', 'E05','E04', 'E03', 'E02', 'E01'};
configcommon.LFP.name              = configcommon.name;

configcommon.muse.backupdir            = fullfile(datasavedir,'Backup_MuseMarker');
configcommon.muse.startmarker.WoD      = 'Vent_Off';   % start and end Muse marker. For defining trials
configcommon.muse.endmarker.WoD        = 'Vent_On';   % start and end Muse marker. For defining trials
configcommon.epoch.toi.WoD             = [-600, 3500];%modif
configcommon.epoch.pad.WoD             = 5;
configcommon.LFP.resamplefs            = 320;%Hz, si possible un diviseur entier de la Fs (sampling frequency) d'origine
configcommon.LFP.write                 = true; %save computed data to disk
configcommon.LFP.layer_depth           = {200,600,800,2000};



configcommon.LFP.lpfilter_wod_detection         = 4;%Hz
configcommon.LFP.wod_toisearch                  = [-1 50]; %s, were to find wod negative peak relative to the muse marker
configcommon.LFP.wor_toisearch                  = [-1 25]; %s, were to find wor positive peak relative to the muse marker
configcommon.LFP.hpfilter_wod_exclusion         = 1; %Hz

configcommon.timefreq.foi          = [1:0.5:100]; %is ritght value
configcommon.timefreq.foi_band     = {[1 9],[10 100]};%Hz
configcommon.timefreq.t_ftimwin    = 1; % in seconds
configcommon.timefreq.timestep     = 0.5;% in second, time between 2 sliding time windows. can be 'all' 2.5 right value
configcommon.timefreq.movmeanwin   = [1,1,1,100,100,100];%in sample points, one value per analysis_name
configcommon.timefreq.tapsmofrq    = 2; %2
configcommon.timefreq.toi          = [-600 600];  % time window suroundging time of interst per trial 
configcommon.timefreq.HF           = [10 100];
configcommon.timefreq.LF           = [1 9];

configcommon.circus.paramfile    = fullfile(script_path,'SpykingCircus_Sofia.params');
configcommon.circus.writedeadfile  = 'no';
configcommon.circus.reref        = 'no';
configcommon.circus.refchan      = [];
configcommon.circus.outputdir    = fullfile(rootpath_analysis, 'data', 'SpykingCircus');
configcommon.circus.hpfilter     = 'no'; % hp before writing data for SC, does not change the hp of SC
configcommon.circus.hpfreq       = 0; % even when not using
configcommon.circus.postfix      = []; % after using circus-gui-matlab's SAVE number



%% rat 1
config{1}                     = configcommon;
config{1}.datasavedir         = datasavedir;       %path where to save MuseStruct data
config{1}.imagesavedir        = imagesavedir;
config{1}.imagesavedir_data   = imagesavedir_data;
config{1}.prefix              = 'Rat-2021_03_30';                                              %patient name. Must end by "-". namepatient-
config{1}.rawdir              = fullfile(rootpath_data,'2021_03_30_WOD');                       %path to patient data

config{1}.directorylist{1}    = {'2021-03-30_14-10', '2021-03-30_15-22'}; %liste de tous les fichiers, tous les protocoles
%config{1}.LFP.channel         = {[], [], 'E16LFP', 'E15LFP', 'E14LFP', 'E13LFP', 'E12LFP', 'E11LFP', 'E10LFP', 'E09LFP', 'E08LFP', 'E07LFP', 'E06LFP', 'E05LFP', 'E04LFP', 'E03LFP','Puff'};
config{1}.LFP.channel         = { 'E29', 'E28','E27', 'E26', 'E25', 'E24','E23', 'E22', 'E21', 'E20', 'E19', 'E18','E17','E16', 'E15', 'E14', 'E13', 'E12', 'E11','E10', 'E09', 'E08'};
config{1}.LFP.rename          = {'E32','E31','E30', 'E29', 'E28','E27', 'E26', 'E25', 'E24','E23', 'E22', 'E21', 'E20', 'E19', 'E18','E17','E16', 'E15', 'E14', 'E13', 'E12', 'E11'};
config{1}.LFP.chan_depth      = {97.2, 197.2, 297.2, 397.2,	497.2, 597.2, 697.2, 797.2,	897.2,	997.2,	1097.2,	1197.2,	1297.2,	1397.2,	1497.2,	1597.2,	1697.2,	1797.2,	1897.2,	1997.2,	2097.2,	2197.2};
config{1}.LFP.origin_WoD          = {'E14', 'E14'};
config{1}.LFP.origin_WoR          = {'E11', 'E11'};
config{1}.LFP.recov               = {1,1};


config{1}.circus.channel      = {'E02', 'E04', 'E11'};       %name of the first electrode
config{1}.circus.rename       = {'E2', 'E4', 'E11'};       %name of the first electrode

%% rat 2
config{2}                     = configcommon;
config{2}.datasavedir         = datasavedir;       %path where to save MuseStruct data
config{2}.imagesavedir        = imagesavedir;
config{2}.imagesavedir_data   = imagesavedir_data;
config{2}.prefix              = 'Rat-2021_02_10';                                              %patient name. Must end by "-". namepatient-
config{2}.rawdir              = fullfile(rootpath_data,'2021_02_10_WOD');                       %path to patient data
config{2}.directorylist{1}    = {'2021-02-10_15-22','2021-02-10_17-22','2021-02-10_18-34'}; %liste de tous les fichiers, tous les protocoles
%config{2}.LFP.channel         = {[], [], 'E16LFP', 'E15LFP', 'E14LFP', 'E13LFP', 'E12LFP', 'E11LFP', 'E10LFP', 'E09LFP', 'E08LFP', 'E07LFP', 'E06LFP', 'E05LFP', 'E04LFP', 'E03LFP','Puff'};
config{2}.LFP.channel         = { 'E31','E30','E29', 'E28','E27', 'E26', 'E25', 'E24','E23', 'E22', 'E21', 'E20', 'E19', 'E18','E17','E16', 'E15', 'E14', 'E13', 'E12'};
config{2}.LFP.rename          = {'E31','E30', 'E29', 'E28','E27', 'E26', 'E25', 'E24','E23', 'E22', 'E21', 'E20', 'E19', 'E18','E17','E16', 'E15', 'E14', 'E13', 'E12'};
config{2}.LFP.chan_depth      = { 247, 347, 447, 547, 647, 747,	847, 947, 1047,	1147, 1247,	1347, 1447,	1547, 1647,	1747, 1847,	1947, 2047,	2147};
config{2}.LFP.origin_WoD          = {'E14', 'E14'};
config{2}.LFP.origin_WoR          = {'E11', 'E11'};
config{2}.LFP.recov               = {1,1};


config{2}.circus.channel      = {'E02', 'E04', 'E11'};       %name of the first electrode
config{2}.circus.rename       = {'E2', 'E4', 'E11'};


%% rat 3
config{3}                      = configcommon;
config{3}.datasavedir         = datasavedir;       %path where to save MuseStruct data
config{3}.imagesavedir        = imagesavedir;
config{3}.imagesavedir_data   = imagesavedir_data;
config{3}.prefix              = 'Rat-2021_03_02';                                              %patient name. Must end by "-". namepatient-
config{3}.rawdir              = fullfile(rootpath_data,'2021_03_02_WOD');                       %path to patient data
config{3}.directorylist{1}    = {'2021-03-02_WOD1','2021-03-02_B2','2021-03-02_WOD2'}; %liste de tous les fichiers, tous les protocoles
%config{6}.LFP.channel         = {[], [], 'E16LFP', 'E15LFP', 'E14LFP', 'E13LFP', 'E12LFP', 'E11LFP', 'E10LFP', 'E09LFP', 'E08LFP', 'E07LFP', 'E06LFP', 'E05LFP', 'E04LFP', 'E03LFP','Puff'};
config{3}.LFP.channel         = { 'E32','E31','E30','E29', 'E28','E27', 'E26', 'E25', 'E24','E23', 'E22', 'E21', 'E20', 'E19', 'E18','E17','E16', 'E15', 'E14', 'E13'};
config{3}.LFP.rename          = {'E30', 'E29', 'E28','E27', 'E26', 'E25', 'E24','E23', 'E22', 'E21', 'E20', 'E19', 'E18','E17','E16', 'E15', 'E14', 'E13', 'E12', 'E11'};
config{3}.LFP.chan_depth      = {279, 379, 479, 579, 679, 779, 879, 979, 1079, 1179, 1279, 1379, 1479, 1579, 1679, 1779, 1879, 1979, 2079, 2179};
config{3}.LFP.origin_WoD          = {'E14', 'E14'};
config{3}.LFP.origin_WoR          = {'E11', 'E11'};
config{3}.LFP.recov               = {1,0};
config{3}.epoch.toi.WoD           = [-600 2000];
config{3}.timefreq.toi           = [-600 500];
config{3}.circus.channel      = {'E11', 'E13', 'E20' ,'E21' ,'E22' ,'E23', 'E28', 'E30'};       %name of the first electrode
config{3}.circus.rename       = {};


%% rat 4
config{4}                        = configcommon;
config{4}.datasavedir         = datasavedir;       %path where to save MuseStruct data
config{4}.imagesavedir        = imagesavedir;
config{4}.imagesavedir_data   = imagesavedir_data;
config{4}.prefix              = 'Rat-2021_03_09';                                              %patient name. Must end by "-". namepatient-
config{4}.rawdir              = fullfile(rootpath_data,'2021_03_09_WOD');                       %path to patient data
config{4}.directorylist{1}    = {'2021-03-09_16-12','2021-03-09_19-18','2021-03-09_19-24'}; %liste de tous les fichiers, tous les protocoles
%config{6}.LFP.channel         = {[], [], 'E16LFP', 'E15LFP', 'E14LFP', 'E13LFP', 'E12LFP', 'E11LFP', 'E10LFP', 'E09LFP', 'E08LFP', 'E07LFP', 'E06LFP', 'E05LFP', 'E04LFP', 'E03LFP','Puff'};
config{4}.LFP.channel         = { 'E30','E29', 'E28','E27', 'E26', 'E25', 'E24','E23', 'E22', 'E21', 'E20', 'E19', 'E18','E17','E16', 'E15', 'E14', 'E13', 'E12'};
config{4}.LFP.rename          = {'E29', 'E28','E27', 'E26', 'E25', 'E24','E23', 'E22', 'E21', 'E20', 'E19', 'E18','E17','E16', 'E15', 'E14', 'E13', 'E12', 'E11'};
config{4}.LFP.chan_depth      = {400, 500, 600, 700, 800, 900, 1000 , 1100, 1200, 1300, 1400, 1500, 1600, 1700, 1800, 1900, 2000, 2100, 2200};
config{4}.LFP.origin_WoD          = {'E14', 'E14'};
config{4}.LFP.origin_WoR          = {'E11', 'E11'};
config{4}.LFP.recov               = {1,0};
config{4}.epoch.toi.WoD           =  [-600 400]; % time window suroundging time of interst trial 2 smaller than configcommon
config{4}.timefreq.toi           = [-600 400];

config{4}.circus.channel      = {'E03' ,'E04', 'E06', 'E09', 'E13' ,'E16', 'E18', 'E19', 'E20', 'E21' ,'E22', 'E23'};       %name of the first electrode
config{4}.circus.rename       = {};


%% rat 5
config{5}                        = configcommon;
config{5}.datasavedir         = datasavedir;       %path where to save MuseStruct data
config{5}.imagesavedir        = imagesavedir;
config{5}.imagesavedir_data   = imagesavedir_data;
config{5}.prefix              = 'Rat-2021_08_13';                                              %patient name. Must end by "-". namepatient-
config{5}.rawdir              = fullfile(rootpath_data,'2021_08_13_WOD');                       %path to patient data
config{5}.directorylist{1}    = {'2021-08-13_17-23','2021-08-13_18-43'}; %liste de tous les fichiers, tous les protocoles
%config{5}.LFP.channel         = {[], [], 'E16LFP', 'E15LFP', 'E14LFP', 'E13LFP', 'E12LFP', 'E11LFP', 'E10LFP', 'E09LFP', 'E08LFP', 'E07LFP', 'E06LFP', 'E05LFP', 'E04LFP', 'E03LFP','Puff'};
config{5}.LFP.channel         = { 'E32LFP','E31LFP', 'E30LFP','E29LFP', 'E28LFP', 'E27LFP', 'E26LFP','E25LFP', 'E24LFP', 'E23LFP', 'E22LFP', 'E21LFP', 'E20LFP', 'E19LFP', 'E18LFP', 'E17LFP', 'E16LFP'};
config{5}.LFP.rename          = {'E32','E31','E30','E29', 'E28','E27', 'E26', 'E25', 'E24','E23', 'E22', 'E21', 'E20','E19','E18','E17','E16'};
config{5}.LFP.chan_depth      = {110, 210, 310, 410, 510, 610, 710 , 810, 910, 1010, 1110, 1210, 1310, 1410, 1510, 1610, 1710, 1810, 1910,2010,2110,2210,2310};
config{5}.LFP.origin_WoD          = {};
config{5}.LFP.origin_WoR          = {};
config{5}.LFP.recov               = {1};


config{5}.circus.channel      = {'E16' ,'E17', 'E19'};       %name of the first electrode
config{5}.circus.rename       = {};

%% rat 6

%% rat 7

%% rat 8
