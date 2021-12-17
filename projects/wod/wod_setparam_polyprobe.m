% set params pour probes differente
%necessite 2 tableau:
%expe  et cut

function [config] = wod_setparam_polyprobe

disp('setting parameters');

if ismac
    error('Platform not supported')
elseif isunix
    rootpath_analysis        = '/network/lustre/iss01/charpier/analyses/wod/Sofia';
    rootpath_data            = '/network/lustre/iss01/charpier/raw/rat-wod/2probes/Extra_Neuralynx';
    rootpath_concatdata      = '/network/lustre/iss01/charpier/raw/rat-wod/2probes/LFP';
    Table_experience_folder  = '/network/lustre/iss01/charpier/raw/rat-wod/2probes/Data_infos/Tableau_Experience_per_rat.xlsx';
    maxicut_folder           = '/network/lustre/iss01/charpier/raw/rat-wod/2probes/Data_infos/maxi_cut';
    os                       = 'unix';
elseif ispc
    rootpath_analysis	     = '\\lexport\iss01.charpier\analyses\wod\Sofia';
    rootpath_data            = '\\lexport\iss01.charpier\raw\rat-wod\2probes\Extra_Neuralynx';
    rootpath_concatdata      = '\\lexport\iss01.charpier\raw\rat-wod\2probes\LFP';
    Table_experience_folder  = '\\lexport\iss01.charpier\raw\rat-wod\2probes\Data_infos\Tableau_Experience_per_rat.xlsx';
    maxicut_folder           = '\\lexport\iss01.charpier\raw\rat-wod\2probes\Data_infos\maxi_cut';
    os                       = 'windows';
else
    error('Platform not supported')
end

script_path  = mfilename('fullpath');
script_path  = fileparts(script_path);

%% load les 2 tableaux de données
rat = dir(maxicut_folder);
tab = dir( maxicut_folder);
rat=rat(3:end);
tab=tab(3:end);

%% config common for all rats

configcommon = struct;

configcommon.datasavedir        =  fullfile(rootpath_analysis,'data');
imagesavedir                    =  fullfile(rootpath_analysis,'images');
configcommon.imagesavedir       = imagesavedir;
configcommon.imagesavedir_data  = {fullfile(imagesavedir,'TFR'), fullfile(imagesavedir,'band_depth'),fullfile(imagesavedir,'TFR','Log'),fullfile(imagesavedir,'LFHF_ratio')};
configcommon.muse.templatemarker= fullfile(configcommon.datasavedir ,'TemplateEventsWOD.mrk');%find the template file to create muse marker file

%configcommon.name                  = {'WoD_short', 'WoD_long'};%, 'WOD_morpho'
%configcommon.LFP.allchannel        = {'E32','E31','E30', 'E29', 'E28','E27', 'E26', 'E25', 'E24','E23', 'E22', 'E21', 'E20', 'E19', 'E18','E17','E16', 'E15', 'E14', 'E13', 'E12', 'E11','E10', 'E09', 'E08', 'E07', 'E06', 'E05','E04', 'E03', 'E02', 'E01','E16b', 'E15b', 'E14b', 'E13b', 'E12b', 'E11b','E10b', 'E09b', 'E08b', 'E07b', 'E06b', 'E05b','E04b', 'E03b', 'E02b', 'E01b'};%,'Events_0001'
configcommon.LFP.allchannel        = {'E32LFP','E31LFP','E30LFP', 'E29LFP', 'E28LFP','E27LFP', 'E26LFP', 'E25LFP', 'E24LFP','E23LFP', 'E22LFP', 'E21LFP', 'E20LFP', 'E19LFP', 'E18LFP','E17LFP','E16LFP', 'E15LFP', 'E14LFP', 'E13LFP', 'E12LFP', 'E11LFP','E10LFP', 'E09LFP', 'E08LFP', 'E07LFP', 'E06LFP', 'E05LFP','E04LFP', 'E03LFP', 'E02LFP', 'E01LFP','E16bLFP', 'E15bLFP', 'E14bLFP', 'E13bLFP', 'E12bLFP', 'E11bLFP','E10bLFP', 'E09bLFP', 'E08bLFP', 'E07bLFP', 'E06bLFP', 'E05bLFP','E04bLFP', 'E03bLFP', 'E02bLFP', 'E01bLFP'};%,'Events_0001'
configcommon.LFP.classe            = {1,2};

%configcommon.LFP.name              = configcommon.name;
configcommon.LFP.resamplefs        = 320;%Hz, si possible un diviseur entier de la Fs (sampling frequency) d'origine
configcommon.LFP.write             = true; %save computed data to disk

% wod_short
configcommon.muse.startmarker.WoD_short      = 'Vent_Off';   % start and end Muse marker. For defining trials
configcommon.muse.endmarker.WoD_short        = 'Vent_On';   % start and end Muse marker. For defining trials
configcommon.epoch.toi.WoD_short             = [-600, 100];
configcommon.epoch.pad.WoD_short             = 0;

%wod_long
configcommon.muse.startmarker.WoD_long      = 'Vent_Off';   % start and end Muse marker. For defining trials
configcommon.muse.endmarker.WoD_long        = 'Vent_On';   % start and end Muse marker. For defining trials
configcommon.epoch.toi.WoD_long             = [-900, 3400];
% configcommon.LFP.bsfilter 					= 'no';%'yes'
% configcommon.LFP.bsfreq 					= [49 51];
configcommon.epoch.pad.WoD_long             = 0;

% %wod_morpho => EXEMPLE, à retirer
% configcommon.muse.startmarker.WOD_morpho      = 'WOD';   % start and end Muse marker. For defining trials
% configcommon.muse.endmarker.WOD_morpho        = 'WOD';   % start and end Muse marker. For defining trials
% configcommon.epoch.toi.WOD_morpho             = [-60, 60];
% configcommon.epoch.pad.WOD_morpho             = 0;

configcommon.LFP.lpfilter_wod_detection         = 7;%Hz
configcommon.LFP.wod_toisearch                  = [-1 100];%[-1 50]; %s, were to find wod negative peak relative to the muse marker
configcommon.LFP.wor_toisearch                  = [-1 25]; %s, were to find wor positive peak relative to the muse marker
configcommon.LFP.hpfilter_wod_exclusion         = 1; %Hz
%configcommon.LFP.recov = {0 0}; %pour chaque trial, 0 si pas de recovery et 1 si recovery => à régler pour chaque rat

configcommon.timefreq.foi              = 1:2:100;%[1:2:100] is ritght value
configcommon.timefreq.foi_band         = {[1 5],[10 20],[25 50],[70 90]};%Hz
configcommon.timefreq.t_ftimwin  = 1;%4; % in seconds
configcommon.timefreq.timestep   = 1;% in second, time between 2 sliding time windows. can be 'all'
% configcommon.timefreq.movmeanwin   = [1,1,1,100,100,100];%in sample points, one value per analysis_name
configcommon.timefreq.tapsmofrq   = 1; %2
configcommon.timefreq.toi         = [-600 200];
configcommon.timefreq.LF  = [1 4];
configcommon.timefreq.TF  = [5 13];
configcommon.timefreq.HF  = [14 45];
configcommon.timefreq.VHF = [55 100];


% configcommon.circus.paramfile    = fullfile(script_path,'SpykingCircus_Sofia.params');
% configcommon.circus.writedeadfile  = 'no';
% configcommon.circus.reref        = 'no';
% configcommon.circus.refchan      = [];
% configcommon.circus.outputdir    = fullfile(rootpath_analysis, 'data', 'SpykingCircus');
% configcommon.circus.hpfilter     = 'no'; % hp before writing data for SC, does not change the hp of SC
% configcommon.circus.hpfreq       = 0; % even when not using
% configcommon.circus.postfix      = []; % after using circus-gui-matlab's SAVE number



for irat =1 : size(rat,1)
    
    config{irat} = configcommon;
    
    %select pertinent data
    if irat==3||irat==4
        config{irat}.to_ignore = true;
    end

    % load les données par rats
    Table_experience = readtable(Table_experience_folder);
    Table_experience = Table_experience(irat,:);
    maxicut(irat) = load([maxicut_folder,'/',tab(irat).name]);
    
    
    config{irat}.trial= Table_experience{1,3};
    config{irat}.prefix = [cell2mat(Table_experience{1,1})];
    config{irat}.prefix = ['Rat-' config{irat}.prefix(1:end-4)];
    
    
    
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
   
    for element=1:length(MUA_A)
        chan_A(element)= str2num(cell2mat (erase(MUA_A(element),"E")));
        chanfind=find(maxicut(irat).maxicut(:,7) == chan_A(element)+16);
        if isempty(chanfind)
            row_A(element)= NaN;
        else
            row_A(element)=find(maxicut(irat).maxicut(:,7) == chan_A(element)+16);
        if ~isnan(maxicut(irat).maxicut(row_A(element),7))
            electrode_A{element} = ['E' num2str((maxicut(irat).maxicut(row_A(element),14)))];
        else
            electrode_A{element}=[];
        end
        end
    end
    
    for element=1:length(MUA_B)
        %MUA_Bb{i} = [cell2mat(MUA_B(i)) 'b']
        chan_B(element)= str2num(cell2mat (erase(MUA_B(element),"E")));
        chanfind=find(maxicut(irat).maxicut(:,7)== chan_B(element));
        if isempty(chanfind)
            row_B(element)= NaN;
        else
        row_B(element)= find(maxicut(irat).maxicut(:,7)== chan_B(element));
        
        if ~isnan(maxicut(irat).maxicut(row_B(element),7))
            electrode_B{element} = ['E' num2str((maxicut(irat).maxicut(row_B(element),14))) 'b'];
        else
            electrode_B{element}=[];
        end
        end
    end
    
    %     config{irat}.circus.channel = [MUA_A MUA_Bb];
    %     config{irat}.circus.rename  = [electrode_A electrode_B];
    
    %% organisation des electrodes
    
    config{irat}.LFP.substructure  = [];
    config{irat}.LFP.chan_depth    = [];
    config{irat}.LFP.deviation     = [];
    
    
    
    if irat>10
        firstchan=2;
    else
        firstchan=1;
    end
    
    for ichan=firstchan:size(maxicut(irat).maxicut,1)
        
        if irat>10
            ichancorrect= ichan;%-1;
        else 
            ichancorrect=ichan;
        end
        %categories pour les structure
        STRC                                         = maxicut(irat).maxicut(ichan,15).';
        config{irat}.LFP.catego{ichancorrect}        = STRC;
        config{irat}.LFP.structure{ichancorrect}     = cellstr(categorical(STRC , [1 2 3 4 5 6], {'CC' 'HPC' 'NC' 'PTA' 'S1' 'TH' }));
        config{irat}.LFP.substructure(ichancorrect)  = maxicut(irat).maxicut(ichan,16);
        config{irat}.LFP.chan_depth(ichancorrect)    = maxicut(irat).maxicut(ichan,13);
        config{irat}.LFP.deviation(ichancorrect)     = maxicut(irat).maxicut(ichan,11);
        

        
       
        if ichan<=32
            if  ~isnan(maxicut(irat).maxicut(ichan,6))
                numchan = maxicut(irat).maxicut(ichan,7)-16;%changer 16 par longueur-32
                if numchan<10
                    config{irat}.LFP.probeA.channel{ichan}  = ['E0' num2str(numchan) 'LFP'];
                else
                    config{irat}.LFP.probeA.channel{ichan}  = ['E' num2str(numchan) 'LFP'];
                end
                numrename = maxicut(irat).maxicut(ichan,14);
                if numrename<10
                    config{irat}.LFP.probeA.rename{ichan}   = ['E0'  num2str(numrename) 'LFP'];
                else
                    config{irat}.LFP.probeA.rename{ichan}   = ['E'  num2str(numrename) 'LFP'] ;
                end
            else
                config{irat}.LFP.probeA.channel{ichan}  = 'NaN';%[];
                config{irat}.LFP.probeA.rename{ichan}   = 'NaN';%[];
            end
        end
        if ichan>32
            if  ~isnan(maxicut(irat).maxicut(ichan,6))
                numchan = maxicut(irat).maxicut(ichan,7);
                if numchan<10
                    config{irat}.LFP.probeB.channel{ichan}  = ['E0' num2str(numchan) 'bLFP' ];
                else
                    config{irat}.LFP.probeB.channel{ichan}  = ['E' num2str(numchan) 'bLFP' ];
                end
                numrename = maxicut(irat).maxicut(ichan,14);
                if numrename<10
                    config{irat}.LFP.probeB.rename{ichan}   = ['E0' num2str(numrename) 'bLFP' ];
                else
                    config{irat}.LFP.probeB.rename{ichan}   = ['E' num2str(numrename) 'bLFP' ];
                end
                
            else
                config{irat}.LFP.probeB.channel{ichan}  = 'NaN';%[];
                config{irat}.LFP.probeB.rename{ichan}   = 'NaN';
            end
            
        end
        
        
      
            
        %         config{irat}.datasavedir{ichan}                  =  fullfile(rootpath_analysis,'data', config{irat}.prefix, config{irat}.LFP.structure{ichan});
        %         config{irat}.imagesavedir{ichan}                 =  fullfile(rootpath_analysis,'images',config{irat}.prefix, config{irat}.LFP.structure{ichan} );
        %config{irat}.imagesavedir_data{ichan}            = {fullfile(config{irat}.imagesavedir{ichan},'TFR'), fullfile(config{irat}.imagesavedir{ichan},'band_depth'),fullfile(config{irat}.imagesavedir{ichan},'band_cx')};
        %         config{irat}.concatdata_path{ichan}              = fullfile(rootpath_concatdata, config{irat}.prefix, config{irat}.LFP.structure{ichan}) ;
        config{irat}.concatdata_path                     = fullfile(rootpath_concatdata) ;
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
        files  = files(3:end);
        % files(irat).name(1:end-4) = config{irat}.prefix(6:end)
        config{irat}.rawdir            = fullfile(rootpath_data, files(irat).name);
        dirlist_temp                   = dir(fullfile(config{irat}.rawdir, '202*'));
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
            if ~dirlist_temp(isubfile).isdir %suleument si dans la structure isdir == 1
                continue
            end
            config{irat}.directorylist{1}{end+1} = dirlist_temp(isubfile).name;
        end
        
        
        
        
        
    end
    
    config{irat}.LFP.channel = [config{irat}.LFP.probeA.channel(~cellfun('isempty',config{irat}.LFP.probeA.channel)) config{irat}.LFP.probeB.channel(~cellfun('isempty',config{irat}.LFP.probeB.channel))];
    config{irat}.LFP.rename  = [config{irat}.LFP.probeA.rename(~cellfun('isempty',config{irat}.LFP.probeA.rename)) config{irat}.LFP.probeB.rename(~cellfun('isempty',config{irat}.LFP.probeB.channel))];
 
    if irat>10
    config{irat}.LFP.channel = ['NaN'  config{irat}.LFP.channel];    
    config{irat}.LFP.rename  = ['NaN'  config{irat}.LFP.rename];
    end
    
    %%  
  if any(contains(config{irat}.LFP.rename, 'ENaNLFP'))||any(contains(config{irat}.LFP.rename, 'ENaNbLFP'))
     index= find(contains(config{irat}.LFP.rename, 'NaN')==1);
     for idx=1:length(index)
           config{irat}.LFP.rename{index(idx)}={'NaN'}; %remplace les ENaNLFP par empty cell
           config{irat}.LFP.channel{index(idx)}={'NaN'};
     end
      
  end 
      
      
      
    
    allchan    =config{irat}.LFP.allchannel;
    selectchan = config{irat}.LFP.channel;
    
%     for element =1:length(allchan) %pour chaque channel de ce rat
%         if ismember(selectchan(element),allchan(element)) % si ce channel correspont au channel de reference alors continue 
%             continue
%         end
%         
%            if element == 1 % si dés le premier c'est pas le bon met directement un nan
%                selectchan= [ NaN selectchan];
%            else
%                selectchan= [selectchan(1:element-1) NaN selectchan(element:end)];% si le suivant n'est pas le bon
%            end
%            
%        
%     end
%     
%     config{irat}.LFP.channel=selectchan
        
  %string remplace les nan par des missing => donc trouve les missing
   
  idx = find((cellfun(@(selectchan) any(ismember(selectchan,'NaN')),selectchan))==1);


 %%%%  
   for ichan=1:size(config{irat}.LFP.channel,2)
       for iidx=1:length(idx)
           if ichan  == idx(iidx)
               config{irat}.LFP.channel(idx(iidx))      ={'NaN'};
               config{irat}.LFP.rename(idx(iidx))       ={'NaN'};
               config{irat}.LFP.catego(idx(iidx))       ={'NaN'};
               config{irat}.LFP.deviation(idx(iidx))    = NaN;
               config{irat}.LFP.chan_depth(idx(iidx))   = NaN;
               config{irat}.LFP.structure(idx(iidx))    ={'NaN'};
               config{irat}.LFP.substructure(idx(iidx)) = NaN;
           end
       end
   end
   
   channel      =   config{irat}.LFP.channel;
   rename       =   config{irat}.LFP.rename;
   catego       =   config{irat}.LFP.catego;
   structure    =   config{irat}.LFP.structure;
   substructure =   num2cell(config{irat}.LFP.substructure);
   deviation    =   num2cell(config{irat}.LFP.deviation);
   chan_depth   =   num2cell(config{irat}.LFP.chan_depth);
   
   channel(cellfun(@(channel) any(ismember(channel,'NaN')), channel))=[];
   rename(cellfun(@(rename) any(ismember(rename,'NaN')), rename))=[];
   catego(cellfun(@(catego) any(ismember(catego,'NaN')), catego))=[];
   structure(cellfun(@(structure) any(ismember(structure,'NaN')), structure))=[];
   substructure(cellfun(@(substructure) any(isnan(substructure)), substructure))=[];
   deviation(cellfun(@(deviation) any(isnan(deviation)), deviation))=[];
   chan_depth(cellfun(@(chan_depth) any(isnan(chan_depth)), chan_depth))=[];
  
   
   config{irat}.LFP.channel      = channel;
   config{irat}.LFP.rename       = rename;
   config{irat}.LFP.catego       = catego;
   config{irat}.LFP.deviation    = cell2mat(deviation);
   config{irat}.LFP.chan_depth   = cell2mat(chan_depth);%
   config{irat}.LFP.structure    = structure;
   config{irat}.LFP.substructure = cell2mat(substructure);
   
   
    for idir=1:size(config{irat}.directorylist{1},2)
        config{irat}.recov{idir}= nanmean(maxicut(irat).maxicut(:,3+idir-1)); % merde c'est 3 et 4!!!!
        config{irat}.LFP.recov{idir}= config{irat}.recov{idir};
            if config{irat}.recov{idir}==1
            config{irat}.name  = {'WoD_short', 'WoD_long'};
            elseif config{irat}.recov{idir}==0
            config{irat}.name  = {'WoD_short'} ;
            end
            config{irat}.LFP.name = config{irat}.name;
    end
     config{irat}.LFP.name = config{irat}.name;
     
      %condition en fonction du trial
     
        for itrial=1:config{irat}.trial
            if irat==15
                if itrial==2
                config{irat}.chosen_value{itrial} = [80 20]
                else
                config{irat}.chosen_value{itrial} = [10 10]
                end
            else
                config{irat}.chosen_value{itrial} = [10 10]  
            end
        end
        
     
     
    
end





end





