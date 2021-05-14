function ouaba_projet_pierre(slurm_task_id)



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

config = wod_setparams_ouabaine;

ipart= 1;


%% Read LFP rat by rat
if slurm_task_id >0
    for irat = slurm_task_id%size(config,2)
        
        if isempty(config{irat})
            continue
        end
        
        
        
        %find concatenated LFP (see wod_concatenateLFP.m)
        config{irat}.rawdir                = config{irat}.concatdata_path;
        config{irat}.directorylist{ipart}  = {config{irat}.prefix};
        
        %read Muse markers
        MuseStruct               = readMuseMarkers(config{irat},true);
        
        %add new markers AD_post
        %         cfgtemp=[];
        %         cfgtemp.editmarkerfile.toadd={'AD_post_START','AD_post_END'};
        %
        %         MuseStruct= editMuseMarkers(cfgtemp,MuseStruct);
        %         datapath_out = fullfile(config{irat}.rawdir, config{irat}.directorylist{1}{1}, 'Events.mrk');
        %         writeMuseMarkerfile(MuseStruct{1}{1}, datapath_out);
        
        
        
        %read LFP. T0 = Vent_Off. Each trial is one protocol
        %         if exist(name_ica, 'file')
        %           load(name_ica)
        %         else
        LFP = readLFP(config{irat}, MuseStruct, true);
        LFP = LFP{1}.(config{irat}.LFP.name{1}); %remove this 'epicode' organisation for now.
        %end
        
        %         %vérifier qu'il y a bien autant de trials que de marqueurs Vent_Off
        startmarker = config{irat}.muse.startmarker.(config{irat}.LFP.name{1});
        if size(LFP.trial,2) ~= size(MuseStruct{1}{1}.markers.(startmarker).synctime,2)
            error('Not the same number of trials that of marker start for %s. \nCheck that begin/end of each trial is not before start of file or after end of file', config{irat}.prefix(1:end-1));
        end
        
        %rename chans according to their real deepness.
        %the name is in cfg.LFP.channel, and it is renamed with the name at
        %the same index in cfg.LFP.rename
        %16 is surface, 1 is the deepest. 0 is the respi.
        n_chans = size(config{irat}.LFP.allchannel,2);
        for ichan = 1:n_chans
            if any(strcmp(config{irat}.LFP.channel,config{irat}.LFP.allchannel{ichan}))
                %search channel into config
                chan_idx = strcmp(config{irat}.LFP.channel,config{irat}.LFP.allchannel{ichan});
                new_name = config{irat}.LFP.rename{chan_idx};
                %search channel into LFP data to remane it
                %chan_idx = strcmp(LFP.label, config{irat}.LFP.allchannel{ichan});
                LFP.label{chan_idx} = new_name;
            end
        end
    end %irat
end %slurm_task_id
%% Analysis by rat


if slurm_task_id==0
%% Analysis for all rats


stats_all= ouaba_wavedetection(config,true);

%calculated_data= ouaba_propag_analysis(config,true);

% fig = figure;
% plot(1:5, 1:5);
%
% set(fig,'Renderer','Painters');
% print(fig, '-dpdf', figurepath,'-r600');

end %slurm_task_id

end %function