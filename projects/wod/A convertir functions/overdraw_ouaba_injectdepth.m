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

config16= wod_setparams;
config32=wod_setparams_32chan;
config_ouaba= wod_setparams_ouabaine;
cfg=[config16 config32 config_ouaba];


stats_ouaba=ouaba_wavedetection(config_ouaba,false);
stats_wod=wod_wavedetection(config16,false);

stats_all=[stats_wod stats_ouaba];


injectimagepath=fullfile(config_ouaba{3}.imagesavedir,'overdraw_traces','sperate_inject');

if ~isfolder(injectimagepath)
    mkdir(injectimagepath)
end

fig_inject=figure;hold
for irat=3: size(cfg,2)
    
    if isempty(cfg{irat})
        continue
    end
    
    %read Muse markers
    MuseStruct = readMuseMarkers(cfg{irat}, false);
    %save(fullfile(config{irat}.datasavedir,sprintf('%s-MuseStruct.mat',config{irat}.prefix)),'MuseStruct');
    
    temp=load(fullfile(cfg{irat}.datasavedir,sprintf('%sLFP_WoD',cfg{irat}.prefix)));
    LFP=temp.LFP{1}.WoD;
    clear temp
    %rename channels according to depth
    for ichan = 1:size(cfg{irat}.LFP.channel, 2)
        idx = strcmp(cfg{irat}.LFP.channel{ichan}, LFP.label);
        label_renamed{idx} = cfg{irat}.LFP.rename{ichan};
    end
    LFP.label = label_renamed';
    clear label_renamed
    
    %remove breathing and ekg channel
    cfgtemp         = [];
    cfgtemp.channel = {'all', '-E0', '-Respi', '-ECG','-Puff'};
    LFP             = ft_selectdata(cfgtemp, LFP);
    LFP_cleaned     = LFP; %save for later removing of artefacts
    
    
    C={[0.9290 0.6940 0.1250 0.5],[1 0 0 0.5]};
    
    %determine channel closest to 900 µm
    for ichan= 1:size(LFP_cleaned.label,1)
        A(ichan,1)=abs(cfg{irat}.LFP.chan_depth{ichan}-900);
    end %ichan
    
    idx_chan=find(A(:,1)<50);
    
    if cfg{irat}.LFP.inject_depth<1200
        C_plot=C{1};
    else
        C_plot=C{2};
    end
    
    for itrial= 1:size(LFP.trial,2)
        time=LFP_cleaned.time{itrial}-stats_ouaba{irat}.WoD.peak_time(idx_chan,itrial);
        plot(time,LFP_cleaned.trial{itrial}(idx_chan,:),'Color',C_plot,'LineWidth',1);
        
        xlim([-20 20]);
        ylim([-10000 2500]);
 
    end %itrial
    
    
end %irat

 fname_inject= fullfile(injectimagepath,sprintf('inject_separated',irat,itrial));
        dtx_savefigure(fig_inject,fname_inject,'png','pdf','close');