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

figpath=fullfile(config{4}.imagesavedir,'overdraw_traces','tfr_lfp');

if ~isfolder(figpath)
    mkdir(figpath);
end

config= wod_setparams;
analysis_names={'timefreq_wod_blcorrected','timefreq_baseline_blcorrected'};

for irat= 5:16
    
    if isempty(config{irat})
        continue
    end
    
    for idata=1:size(analysis_names,2)
        %read Muse markers
        MuseStruct = readMuseMarkers(config{irat}, false);
        %save(fullfile(config{irat}.datasavedir,sprintf('%s-MuseStruct.mat',config{irat}.prefix)),'MuseStruct');
        
        %read LFP, append electrodes, and cut into trials according to Muse Markers
        LFP = readLFP(config{irat}, MuseStruct, false);
        LFP = LFP{1}.(config{irat}.LFP.name{1}); %remove this 'epicode' organisation for now.
        
        %vÃ©rifier qu'il y a bien autant de trials que de marqueurs Vent_Off
        startmarker = config{irat}.muse.startmarker.(config{irat}.LFP.name{1});
        if size(LFP.trial,2) ~= size(MuseStruct{1}{1}.markers.(startmarker).synctime,2)
            error('Not the same number of trials that of marker start for %s. \nCheck that begin/end of each trial is not before start of file or after end of file', config{irat}.prefix(1:end-1));
        end
        
        %rename channels according to depth
        for ichan = 1:size(config{irat}.LFP.channel, 2)
            idx = strcmp(config{irat}.LFP.channel{ichan}, LFP.label);
            label_renamed{idx} = config{irat}.LFP.rename{ichan};
        end
        LFP.label = label_renamed';
        clear label_renamed
        
        
        
        tfr_temp = load(fullfile(config{irat}.datasavedir,[config{irat}.prefix,analysis_names{idata},'.mat']));
        tfr=tfr_temp.(analysis_names{idata});
        clear tfr_temp
        
        for itrial=1%:size(LFP.trial,2)
            %select trial 1 of LFP and tfr
            cfgtemp=[];
            cfgtemp.trials=itrial;
            LFP_trial= ft_selectdata(cfgtemp,LFP);
            
            tfr= tfr{itrial};
            
            chan_list=LFP.label;
            
            for ichan=[1 3]
                
                chan_name_tfr = sprintf('%sLFP',chan_list{ichan});
                
                data_plot = tfr.(chan_name_tfr);
                
                fig=figure;
                
                subplot(2,1,2)
                %plot TFR
                %voir les paramètres optionnels dans le descriptifs de la fonction pour
                %modifier l'aspect du TFR. Avec les paramètres par défaut :
                cfgtemp         = [];
                cfgtemp.channel = 'all';
                cfgtemp.interactive = 'no';
                cfgtemp.colormap= 'jet';
                cfgtemp.colorbar='no';
                cfgtemp.fontsize = 12;
                cfgtemp.ylim= [1 100];
                cfgtemp.masknans    = 'yes';
                ft_singleplotTFR(cfgtemp, data_plot);
                ft_pimpplot(fig, jet(5000))
                caxis([0 2]);
                xlim([-30 stats_all{irat}.ISO(ichan,itrial)]);
                
                
                subplot(2,1,1)
                plot(LFP_trial.time{1},LFP_trial.trial{1}(ichan,:),'Color','k','LineWidth',1);
                xlim([-30 stats_all{irat}.ISO(ichan,itrial)]);
                
                if ichan==1
                    ylim([-500 500]);
                else
                    ylim([-1000 1000]);
                end
                
                
                fname_fig=fullfile(figpath,sprintf('tfr_lfp_Rat_%i_%s_trial_%i_%s',irat,chan_name_tfr,itrial,analysis_names{idata}));
                dtx_savefigure(fig,fname_fig,'pdf','png','close');
                
                
            end %ichan
        end %itrial
        clear tfr LFP_trial
        
    end %idata
end %irat