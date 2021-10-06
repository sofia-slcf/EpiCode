function wod_LFP(rat_list, configscript)

% Les données avec une fréquence d'échantillonnage supérieure à 30kHz sont
% sous-échantillonnés au dixième.

% code copied from writespykingcircus.m
% Must use a modified version of ft_writedata to give the good channel
% name to Muse. Each added line in ft_writedata is commented with %Paul.

% configscript : name of the setparams script, ie 'wod_setparams', or
% 'wod_setparams_32chans';

%retrouver le chemin du dossier EpiCode par rapport à ce script

try %en local
    scriptpath = matlab.desktop.editor.getActiveFilename;
catch %cluster
    scriptpath = mfilename('fullpath');
end

epicodepath = [fileparts(fileparts(fileparts(fileparts(scriptpath)))), filesep];

addpath (genpath([epicodepath,'shared']))
addpath (genpath([epicodepath,'external']))
addpath (genpath([epicodepath,'templates']))
addpath (genpath([epicodepath,'projects', filesep, 'wod']))

if ispc
    addpath \\lexport\iss01.charpier\analyses\wod\fieldtrip-20200607
elseif isunix
    addpath /network/lustre/iss01/charpier/analyses/wod/fieldtrip-20200607
end

%/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\
%Attention au risque d'écraser des fichiers marqueurs déjà existants
overwriteMuseMarkerFile = false; %false : do not write a new muse marker if one is founded.
%/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\

ft_defaults

%config = eval(configscript);%%% change 23/09/21
config = wod_setparamconfig_polyprob;

for irat = rat_list
    
    for ipart = 1 : size(config{irat}.directorylist,2)
        
    %rename rat folder yyyy_mm_dd_WOD in Rat-yyyy_mm_dd
        datadir   = config{irat}.concatdata_path
        datafolder= dir(datadir)
        datafolder= datafolder(3:end)
        foldername=datafolder(irat).name
            if foldername(1)~=('R')
            newfoldername=config{irat}.prefix
            movefile(fullfile(datadir,foldername),fullfile(datadir,newfoldername))
                %rename trial folder yyyy-mm-dd_hh-mm in Trila_1,2...
                for idir =  1:size(config{irat}.directorylist{ipart},2)
                trialfolder= config{irat}.directorylist{ipart}{idir}
                newtrialfolder=['Trial_' num2str(idir)]
                movefile(fullfile(datadir,newfoldername,trialfolder),fullfile(datadir,newfoldername,newtrialfolder))   
                end
            end
            
    % loop through different parts
            
            output_datapath = fullfile(config{irat}.concatdata_path,config{irat}.prefix)
            if ~isfolder(output_datapath)
                mkdir(output_datapath);
            end
            
            %create Muse Marker File
            add_nev = true;
            if overwriteMuseMarkerFile || ~exist(fullfile(output_datapath,'Events.mrk'), 'file')
                copyfile(config{irat}.muse.templatemarker, fullfile(output_datapath,'Events.mrk')); %écrase markers préexistants
                add_nev = true;
            end
            
      
        
        %% process events
 
            
            temp = dir(fullfile(config{irat}.rawdir,config{irat}.directorylist{ipart}{idir},'*.nev'));
            
            if size(temp,1) > 1
                error('there should be only one nev file')
            end
            if size(temp,1) == 0
                nev_data{idir} = [];
                continue
            end
            
            fname = fullfile(config{irat}.rawdir,config{irat}.directorylist{ipart}{idir},temp.name);
            temp        = dir(fullfile(config{irat}.rawdir,config{irat}.directorylist{ipart}{idir},['*',config{irat}.LFP.channel{ichan},'.ncs']));
            hdrdir{idir}  = ft_read_header(fullfile(config{irat}.rawdir,config{irat}.directorylist{ipart}{idir}, temp.name));
            
            nev_data{idir} = read_neuralynx_nev(fname,'eventformat','neuralynx_nev');
            
            for ievent = 1:size(nev_data{idir}, 1)
               
                nev_data{idir}(ievent, 1).idir = idir;
                %correct event times in case there is missing data between 2
                if idir > 1 & ~isempty(hdrdir{idir-1})
                    timesampdiff(idir) = hdrdir{idir}.FirstTimeStamp - (hdrdir{idir-1}.FirstTimeStamp + hdrdir{idir-1}.nSamples * hdrdir{idir-1}.TimeStampPerSample);
                    nev_data{idir}(ievent).TimeStamp = nev_data{idir}(ievent).TimeStamp - timesampdiff(idir);
                end
                
            end
            
            
            nev_all = nev_data{idir};
           
            for ievent = 1:size(nev_data{idir},1)
                idx = size(nev_all,1) + 1;
                for ifield = string(fieldnames(nev_data{idir})')
                    nev_all(idx, 1).(ifield) = nev_data{idir}(ievent).(ifield);
                end
            end
            
            
            
            %save events' MATLAB structure
            fname = fullfile(output_datapath,'Events_nev.mat');
            save(fname, 'nev_all');
          
            %add events to MuseStruct
            if add_nev
                cfgtemp = config{irat};
                [cfgtemp.rawdir,dir_name]     = fileparts(output_datapath);
               
                cfgtemp.directorylist{ipart}  = {dir_name};
                MuseStruct{idir} = readMuseMarkers(cfgtemp, true);
%                 MuseStruct{idir} = read_nev_Muse(cfgtemp, MuseStruct{idir}, nev_all);
%                 
                 fname = fullfile(cfgtemp.rawdir, cfgtemp.directorylist{1}{1}, 'Events.mrk');
%                 writeMuseMarkerfile(MuseStruct{1}{1}, fname);
            end

           end
            
            
            
            if isempty(idir) %correct idir if there is only 1 file
                idir = 1;
            end
            
            
        end
      
        
    end %ipart
end %irat
end %wod_concatenate