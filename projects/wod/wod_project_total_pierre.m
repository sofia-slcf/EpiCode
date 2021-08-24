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
        MuseStruct               = readMuseMarkers(config{irat},false);
        
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
        LFP = readLFP(config{irat}, MuseStruct, false);
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
    
    config_900 = config(3:5);
    config_2000 = config(6:8);
    % detection data for waves
    stats_all= ouaba_wavedetection(config,false);
    stats_all_900= ouaba_wavedetection(config_900, true);
    stats_all_2000= ouaba_wavedetection(config_2000, true);
    
    %extract origin time, depth and propagation speed
    calculated_data= ouaba_propag_analysis(stats_all,config,false);
    calculated_data900=ouaba_propag_analysis(stats_all_900,config_900,true);
    calculated_data2000=ouaba_propag_analysis(stats_all_2000,config_2000,true);

    
    
    %order all data according to depth
    ordered_data= ouaba_fusion_data(stats_all,config,false);
    ordered_data900= ouaba_fusion_data(stats_all_900,config_900,false);
    ordered_data2000= ouaba_fusion_data(stats_all_2000,config_2000,false);
    
    %plot delays
    
    ouaba_plot_delays(ordered_data,config);
    
    % fig = figure;
    % plot(1:5, 1:5);
    %
    % set(fig,'Renderer','Painters');
    % print(fig, '-dpdf', figurepath,'-r600');
    
end %slurm_task_id



%% Overdraws ouabaine

for irat= 1:size(config,2)

    
if isempty(config{irat})
    continue
end
    
    
temp=load(fullfile(config{irat}.datasavedir,sprintf('%sLFP_WoD',config{irat}.prefix)));
LFP=temp.LFP{1}.WoD;
clear temp

%rename channels according to depth
for ichan = 1:size(config{irat}.LFP.channel, 2)
    idx = strcmp(config{irat}.LFP.channel{ichan}, LFP.label);
    label_renamed{idx} = config{irat}.LFP.rename{ichan};
end
LFP.label = label_renamed';
clear label_renamed


chan_list=LFP.label;

for itrial=1:size(LFP.trial,2)

%find initiation channel index
idx_origin= find(stats_all{irat}.WoD.peak_time(:,itrial)==min(stats_all{irat}.WoD.peak_time(:,itrial)));

%%low-pass filter data
% cfgtemp=[];
% cfgtemp.lpfilter='yes';
% cfgtemp.lpfreq= config{irat}.LFP.lpfilter_wod_detection;
% cfgtemp.lpfilttype='fir';
% LFP_lpfilt= ft_preprocessing(cfgtemp,LFP);
% 

C_1 = autumn(max(size(chan_list,1)-idx_origin+1, idx_origin));
C_2 = flip(C_1(2:idx_origin,:));
C = vertcat(C_2,C_1);

fig_overdraw=figure;
sgtitle(sprintf('Rat_%i_AD_%i_on_%i_overdraw',irat,itrial,size(LFP.trial,2)))
for ichan=1:size(chan_list,1)

plot(LFP.time{itrial},LFP.trial{itrial}(ichan,:),'Color',C(ichan,:))
hold on

xlim([min(stats_all{irat}.WoD.peak_time(:,itrial))-5 max(stats_all{irat}.WoD.peak_time(:,itrial))+15])
end %ichan

fname_overdraw=fullfile(config{irat}.imagesavedir,'overdraw_traces','color_from_origin',sprintf('%s_AD_%i_on_%i',config{irat}.prefix,itrial,size(LFP.trial,2)));

if ~isfolder(fullfile(config{irat}.imagesavedir,'overdraw_traces','color_from_origin'))
    mkdir(fullfile(config{irat}.imagesavedir,'overdraw_traces','color_from_origin'));
end
    


dtx_savefigure(fig_overdraw,fname_overdraw,'pdf','png','close');

%% Traces of AD

fig_multitraces = figure;
sgtitle(sprintf('Rat_%i_AD_%i_on_%i_multitraces',irat,itrial,size(LFP.trial,2)))
d=6000;
a=0;
for ichan=1:size(chan_list,1)

plot(LFP.time{itrial},LFP.trial{itrial}(ichan,:)-(a*d),'Color','k','LineWidth',2);
hold on 
ylabel(sprintf('%i',config{irat}.LFP.chan_depth{ichan}));    
xlim([min(stats_all{irat}.WoD.peak_time(:,itrial))-10 max(stats_all{irat}.WoD.peak_time(:,itrial))+20])
a= a+1;

end %ichan

fname_multitraces=fullfile(config{irat}.imagesavedir,'overdraw_traces','multitraces',sprintf('%s_AD_%i_on_%i',config{irat}.prefix,itrial,size(LFP.trial,2)));

if ~isfolder(fullfile(config{irat}.imagesavedir,'overdraw_traces','multitraces'))
    mkdir(fullfile(config{irat}.imagesavedir,'overdraw_traces','multitraces'));
end

dtx_savefigure(fig_multitraces,fname_multitraces,'pdf','png','close');
end %itrial

end %irat

%% overdraw first and last ouaba and WoD traces




%% Stats
Datapath_wod = '\\lexport\iss01.charpier\analyses\wod\Antoine\data\Detection';
temp = load(fullfile(Datapath_wod, 'calculated_data.mat'));
calculated_data_wod=temp.calculated_data;
clear temp
%depth

depth_ouaba = calculated_data.WoD.origin_depth.peak_time;
depth_ouaba_900 = reshape(depth(3:5,:),1,9);
depth_ouaba_2000 = reshape(depth(6:8,:),1,9);
depth_wod = reshape(calculated_data_wod.WoD.origin_depth.peak_time,[],1);
depth_all = cat(1,depth_900,depth_2000)';

%speed
speed_up = calculated_data.WoD.speed.start_time.up;
speed_down = calculated_data.WoD.speed.start_time.down;
speed_up_900 = reshape(speed_up(3:5,:),1,9);
speed_down_900 = reshape(speed_down(3:5,:),1,9);
speed_up_2000 = reshape(speed_up(6:8,:),1,9);
speed_down_2000 = reshape(speed_down(6:8,:),1,9);

end %function