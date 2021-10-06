% set params pour probes differente
%necessite 2 tableau: 
%expe  et cut

function [config] = wod_setparamconfig_polyprob

disp('setting parameters');

if ismac
    error('Platform not supported')
elseif isunix
    rootpath_analysis        = '/network/lustre/iss01/charpier/analyses/wod/Sofia';
    rootpath_data            = '/network/lustre/iss01/charpier/raw/rat-wod/2probes/Extra_Neuralynx';
    rootpath_concatdata      = '/network/lustre/iss01/charpier/raw/rat-wod/2probes/DATA';
    Table_experience_folder  = '/network/lustre/iss01/charpier/raw/rat-wod/2probes/Data_table/Tableau_Experience_per_rat.xlsx';
    maxicut_folder           = '/network/lustre/iss01/charpier/raw/rat-wod/2probes/Data_table/maxi_cut';
    os                       = 'unix';
elseif ispc
    rootpath_analysis	     = '\\lexport\iss01.charpier\analyses\wod\Sofia';
    rootpath_data            = '\\lexport\iss01.charpier\raw\rat-wod\2probes\Extra_Neuralynx';
    rootpath_concatdata      = '\\lexport\iss01.charpier\raw\rat-wod\2probes\DATA';
    Table_experience_folder  = '\\lexport\iss01.charpier\raw\rat-wod\2probes\Data_table\Tableau_Experience_per_rat.xlsx';
    maxicut_folder           = '\\lexport\iss01.charpier\raw\rat-wod\2probes\Data_table\maxi_cut';
    os                       = 'windows';
else
    error('Platform not supported')
end 

script_path  = mfilename('fullpath');
script_path  = fileparts(script_path);

%% load les 2 tableaux de données
rat = dir(maxicut_folder);
tab = dir( maxicut_folder)
rat=rat(3:end);
tab=tab(3:end);

%% config common for all rats

configcommon = struct;

configcommon.datasavedir                  =  fullfile(rootpath_analysis,'data');
configcommon.imagesavedir                 =  fullfile(rootpath_analysis,'images');
configcommon.muse.templatemarker          = fullfile(configcommon.datasavedir ,'TemplateEventsWOD.mrk');%find the template file to create muse marker file
configcommon.muse.backupdir               = fullfile(configcommon.datasavedir ,'Backup_MuseMarker');

configcommon.name                  = {'WoD_short', 'WoD_long'};
configcommon.LFP.allchannel        = {'E32','E31','E30', 'E29', 'E28','E27', 'E26', 'E25', 'E24','E23', 'E22', 'E21', 'E20', 'E19', 'E18','E17','E16', 'E15', 'E14', 'E13', 'E12', 'E11','E10', 'E09', 'E08', 'E07', 'E06', 'E05','E04', 'E03', 'E02', 'E01','E16b', 'E15b', 'E14b', 'E13b', 'E12b', 'E11b','E10b', 'E09b', 'E08b', 'E07b', 'E06b', 'E05b','E04b', 'E03b', 'E02b', 'E01b','Events_0001'};
configcommon.LFP.name              = configcommon.name;
configcommon.LFP.resamplefs                     = 320;%Hz, si possible un diviseur entier de la Fs (sampling frequency) d'origine
configcommon.LFP.write                          = true; %save computed data to disk

% wod_short
configcommon.muse.startmarker.WoD_short      = 'Vent_Off';   % start and end Muse marker. For defining trials
configcommon.muse.endmarker.WoD_short        = 'Vent_On';   % start and end Muse marker. For defining trials
configcommon.epoch.toi.WoD_short             = [-600, 200];
configcommon.epoch.pad.WoD_short             = 5;

%wod_long
configcommon.muse.startmarker.WoD_long      = 'Vent_Off';   % start and end Muse marker. For defining trials
configcommon.muse.endmarker.WoD_long        = 'Vent_On';   % start and end Muse marker. For defining trials
configcommon.epoch.toi.WoD_long             = [-900, 3400];
configcommon.epoch.pad.WoD_long             = 5;


configcommon.LFP.lpfilter_wod_detection         = 7;%Hz
configcommon.LFP.wod_toisearch                  = [-1 50]; %s, were to find wod negative peak relative to the muse marker
configcommon.LFP.wor_toisearch                  = [-1 25]; %s, were to find wor positive peak relative to the muse marker
configcommon.LFP.hpfilter_wod_exclusion         = 1; %Hz

configcommon.timefreq.foi          = 1:2:100;%[1:2:100] is ritght value
configcommon.timefreq.foi_band     = {[1 5],[10 20],[25 50],[70 90]};%Hz
configcommon.timefreq.t_ftimwin.long    = 10;% in second, length of the time window 10 (right value)
configcommon.timefreq.t_ftimwin.short = 4; % in seconds
configcommon.timefreq.timestep.long     = 2.5;% in second, time between 2 sliding time windows. can be 'all' 2.5 right value
configcommon.timefreq.timestep.short     = 1;% in second, time between 2 sliding time windows. can be 'all'
configcommon.timefreq.movmeanwin   = [1,1,1,100,100,100];%in sample points, one value per analysis_name
configcommon.timefreq.tapsmofrq.long    = 0; %2
configcommon.timefreq.tapsmofrq.short    = 0; %2
configcommon.timefreq.toi.long           = [-900 3400];
configcommon.timefreq.toi.short           = [-900 600];

configcommon.timefreq.toi.short= [-600 200];
configcommon.timefreq.toi.long= [-900, 3400];

configcommon.circus.paramfile    = fullfile(script_path,'SpykingCircus_Sofia.params');
configcommon.circus.writedeadfile  = 'no';
configcommon.circus.reref        = 'no';
configcommon.circus.refchan      = [];
configcommon.circus.outputdir    = fullfile(rootpath_analysis, 'data', 'SpykingCircus');
configcommon.circus.hpfilter     = 'no'; % hp before writing data for SC, does not change the hp of SC
configcommon.circus.hpfreq       = 0; % even when not using
configcommon.circus.postfix      = []; % after using circus-gui-matlab's SAVE number

for irat=1:size(rat,1)

config{irat} = configcommon;
% load les données par rats
Table_experience = readtable(Table_experience_folder)
Table_experience = Table_experience(irat,:)
maxicut(irat)= load([maxicut_folder,'/',tab(irat).name]);

config{irat}.prefix = [cell2mat(Table_experience{1,1})]
config{irat}.prefix = ['Rat-' config{irat}.prefix(1:end-4)]


%MUA pour spiking circus
MUArow_A=cell2mat(Table_experience{1,23});
MUArow_B=cell2mat(Table_experience{1,24});
MUA_A = strsplit(MUArow_A,',');
MUA_B = strsplit(MUArow_B,',');
chan_A=[];
row_A=[];
electrode_A={};
MUA_Bb = {};
chan_B=[];
row_B=[];
electrode_B={};
    for i=1:length(MUA_A)
    chan_A(i)= str2num(cell2mat (erase(MUA_A(i),"E")));
    row_A(i)=find(maxicut(irat).maxicut(:,7)== chan_A(i)+16);
        if ~isnan(maxicut(irat).maxicut(row_A,7))
        electrode_A{i} = ['E' num2str((maxicut(irat).maxicut(row_A(i),14)))];
         else
            electrode_A{i}=[] 
        end
    end
    
    for i=1:length(MUA_B)
    %MUA_Bb{i} = [cell2mat(MUA_B(i)) 'b']   
    chan_B(i)= str2num(cell2mat (erase(MUA_B(i),"E")));
    row_B(i)=find(maxicut(irat).maxicut(:,7)== chan_B(i));
        if ~isnan(maxicut(irat).maxicut(row_B,7))
        electrode_B{i} = ['E' num2str((maxicut(irat).maxicut(row_B(i),14))) 'b'];
        else
            electrode_B{i}=[] 
        end
    end
   
    
    config{irat}.circus.channel = [MUA_A MUA_Bb]
    config{irat}.circus.rename  = [electrode_A electrode_B]

%% organisation des electrodes

    config{irat}.LFP.substructure  = []
    config{irat}.LFP.profondeur    = []
    config{irat}.LFP.deviation     = []

    for ichan=1:size(maxicut(irat).maxicut,1)
        
        %categories pour les structure
        STRC                                  = maxicut(irat).maxicut(ichan,15).'
        config{irat}.LFP.structure{ichan}     = cellstr(categorical(STRC , [1 2 3 4 5 6], {'CC' 'HPC' 'NC' 'PTA' 'S1' 'TH' }))
        config{irat}.LFP.substructure(ichan)  = maxicut(irat).maxicut(ichan,16)
        config{irat}.LFP.profondeur(ichan)    = maxicut(irat).maxicut(ichan,13)
        config{irat}.LFP.deviation(ichan)     = maxicut(irat).maxicut(ichan,11)
        

        if ichan<=32
            if  ~isnan(maxicut(irat).maxicut(ichan,6))
%                 config{irat}.LFP.probeA.channel{ichan}  = ['E' num2str(maxicut(irat).maxicut(ichan,7)-16) 'LFP']
%                 config{irat}.LFP.probeA.rename{ichan}   = ['E' num2str(maxicut(irat).maxicut(ichan,14)) 'LFP']
                
                if length(num2str(maxicut(irat).maxicut(ichan,7)-16))== 1
                   config{irat}.LFP.probeA.channel{ichan}  = ['E0' num2str(maxicut(irat).maxicut(ichan,7)-16) 'LFP']
                elseif length(num2str(maxicut(irat).maxicut(ichan,7)-16))== 2
                   config{irat}.LFP.probeA.channel{ichan}  = ['E' num2str(maxicut(irat).maxicut(ichan,7)-16) 'LFP'] 
                end
                
                if length(num2str(maxicut(irat).maxicut(ichan,14)))== 1
                    config{irat}.LFP.probeA.rename{ichan}   = ['E0' num2str(maxicut(irat).maxicut(ichan,14)) 'LFP']
                elseif length(num2str(maxicut(irat).maxicut(ichan,7)-16))== 2
                    config{irat}.LFP.probeA.rename{ichan}   = ['E' num2str(maxicut(irat).maxicut(ichan,14)) 'LFP'] 
                end
            else
                 config{irat}.LFP.probeA.channel{ichan}  = []
                 config{irat}.LFP.probeA.rename{ichan}   = []
            end
        end
        
        if ichan>32
            if  ~isnan(maxicut(irat).maxicut(ichan,6))
                 if length(num2str(maxicut(irat).maxicut(ichan,7)))== 1
                     config{irat}.LFP.probeB.channel{ichan}  = ['E0' num2str(maxicut(irat).maxicut(ichan,7)) 'bLFP' ]
                 elseif length(num2str(maxicut(irat).maxicut(ichan,7)))== 2
                     config{irat}.LFP.probeB.channel{ichan}  = ['E' num2str(maxicut(irat).maxicut(ichan,7)) 'bLFP' ]
                 end
                
                if length(num2str(maxicut(irat).maxicut(ichan,14)))==1
                     config{irat}.LFP.probeB.rename{ichan}   = ['E0' num2str(maxicut(irat).maxicut(ichan,14)) 'bLFP' ]
                elseif length(num2str(maxicut(irat).maxicut(ichan,14)))==2
                     config{irat}.LFP.probeB.rename{ichan}   = ['E' num2str(maxicut(irat).maxicut(ichan,14)) 'bLFP' ]
                end
                
            else
                config{irat}.LFP.probeB.channel{ichan}  = []
                config{irat}.LFP.probeB.rename{ichan}   = [] 
            end
            
        end
        
%         config{irat}.datasavedir{ichan}                  =  fullfile(rootpath_analysis,'data', config{irat}.prefix, config{irat}.LFP.structure{ichan});
%         config{irat}.imagesavedir{ichan}                 =  fullfile(rootpath_analysis,'images',config{irat}.prefix, config{irat}.LFP.structure{ichan} );
        %config{irat}.imagesavedir_data{ichan}            = {fullfile(config{irat}.imagesavedir{ichan},'TFR'), fullfile(config{irat}.imagesavedir{ichan},'band_depth'),fullfile(config{irat}.imagesavedir{ichan},'band_cx')};
%         config{irat}.concatdata_path{ichan}              = fullfile(rootpath_concatdata, config{irat}.prefix, config{irat}.LFP.structure{ichan}) ;
        config{irat}.concatdata_path                    = fullfile(rootpath_concatdata) ;
%         configcommon.muse.backupdir                      = fullfile(config{irat}.datasavedir{ichan},'TemplateEventsWOD.mrk');%find the template file to create muse marker file
        configcommon.muse.backupdir                      = fullfile(config{irat}.datasavedir,'TemplateEventsWOD.mrk');%find the template file to create muse marker file
      
%        files  =  dir(rootpath_data);
%        files  = files(3:end)
%        files(irat).name(1:end-4) = config{irat}.prefix(6:end)
%        config{irat}.rawdir = fullfile(rootpath_data, files(irat).name)
%         dirlist_temp = dir(fullfile(config{irat}.rawdir, '202*'));
%         config{irat}.filenumber{1}  =  {};
%         for isubfile = 1:size(dirlist_temp, 1)%ne selectionne que les folders contenant les data
%             if ~dirlist_temp(isubfile).isdir %suleument si dans la structure isdir == 1
%                 continue
%             end
%             config{irat}.filenumber{1}{end+1} = dirlist_temp(isubfile).name;
%         end
%         config{irat}.directorylist{1}={config{irat}.filenumber{1}{config{irat}.trial},};
        
   
            
       files  =  dir(rootpath_data);
       files  = files(3:end)
      % files(irat).name(1:end-4) = config{irat}.prefix(6:end)
       config{irat}.rawdir = fullfile(rootpath_data, files(irat).name)
        dirlist_temp = dir(fullfile(config{irat}.rawdir, '202*'));
        config{irat}.directorylist{1}  =  {};
        %config{irat}.filenumber{1}  =  {};
%         for isubfile = 1:size(dirlist_temp, 1)%ne selectionne que les folders contenant les data
%             if ~dirlist_temp(isubfile).isdir %suleument si dans la structure isdir == 1
%                 continue
%             end
%             config{irat}.filenumber{1}{end+1} = dirlist_temp(isubfile).name;
%         end
%             config{irat}.directorylist{1}={config{irat}.filenumber{1}{config{irat}.trial},};
%         
         for isubfile = 1:size(dirlist_temp, 1)%ne selectionne que les folders contenant les data
            if ~dirlist_temp(isubfile).isdir %seuleument si dans la structure isdir == 1
                continue
            end
         config{irat}.directorylist{1}{end+1} = dirlist_temp(isubfile).name;
         end
         
       
        
            
        
        
    end

   config{irat}.LFP.channel =[config{irat}.LFP.probeA.channel(~cellfun('isempty',config{irat}.LFP.probeA.channel)) config{irat}.LFP.probeB.channel(~cellfun('isempty',config{irat}.LFP.probeB.channel))]
   config{irat}.LFP.rename = [config{irat}.LFP.probeA.rename(~cellfun('isempty',config{irat}.LFP.probeA.rename)) config{irat}.LFP.probeB.rename(~cellfun('isempty',config{irat}.LFP.probeB.channel))]
  



end





end





