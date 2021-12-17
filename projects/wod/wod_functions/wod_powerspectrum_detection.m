function [power_detection]=wod_powerspectrum_detection(config, rongeur, per_rat, analysis_list, force)

% plot les powerspectrum ( en frequence ou en range de frequences)pour chaque rat en fonction des
%differnetes régions


analysis_names = {'timefreq_wod', 'timefreq_wod_timenorm', 'timefreq_baseline','timefreq_wod_blcorrected', 'timefreq_wod_timenorm_blcorrected', 'timefreq_baseline_blcorrected','log_timefreq_wod', 'log_timefreq_wod_timenorm', 'log_timefreq_baseline','log_timefreq_wod_blcorrected', 'log_timefreq_wod_timenorm_blcorrected','log_timefreq_baseline_blcorrected'};
brain_struct = {'CC','HPC', 'NC', 'PTA', 'S1', 'TH'};
brain_struct_color = linspecer(length(brain_struct));

if ispc
    savedir = '\\lexport\iss01.charpier\analyses\wod\Sofia\images';
elseif isunix
    savedir = '/network/lustre/iss01/charpier/analyses/wod/Sofia/images';
end

do_plot=false; %false

% for idata = 4 %1:length(analysis_names);
for analysis = analysis_list %"timefreq_wod_blcorrected"
    
    %name of the output data saved on the disk :
    fname = fullfile(config{1}.datasavedir, sprintf('allrats_power_detection_%s.mat', analysis));
    
    %load precomputed data if asked :
    if exist(fname,'file') && force == false
        fprintf('Loading precomputed data results\n');
        load(fname, 'power_detection');
        return
    end
    
%     for ifield=1:6
%         power_detection(ifield).peak_LF=[];
%         power_detection(ifield).peak_TF=[];
%         power_detection(ifield).peak_HF=[];
%         power_detection(ifield).peak_VHF=[];
%         power_detection(ifield).globale=[];
%     end
    
    for irat=1:length(config)
        
        %remove rats I want to ignore
        config{irat}.to_ignore = ft_getopt(config{irat}, 'to_ignore', false);
        if config{irat}.to_ignore
            continue
        end
        
        % load precomputed data:
        cfgtemp    = [];
        MuseStruct = readMuseMarkers(config{irat}, false);
        LFP        = readLFP(config{irat}, MuseStruct, false);
        data_tfr   = load(fullfile(config{irat}.datasavedir,[config{irat}.prefix,analysis,'.mat']));
        
        %pour chaque trial
        for itrial=1:size(config{irat}.directorylist{1},2)
            
            %% clean and organize the data
            %remove empty channels
            for ifield = string(fieldnames(data_tfr.(analysis){itrial}))'
                if isempty(data_tfr.(analysis){itrial}.(ifield))
                    data_tfr.(analysis){itrial} = rmfield(data_tfr.(analysis){itrial}, ifield);
                end
            end
            
            %remove 50 Hz %FIXME à regarder
            for ifield = string(fieldnames(data_tfr.(analysis){itrial}))'
                for channame = string(fieldnames(data_tfr.(analysis){itrial}))'
                    temp = data_tfr.(analysis){itrial}.(channame);
                    freq_sel          = temp.freq < 49 | temp.freq > 51;
                    data_tfr.(analysis){itrial}.(channame).freq      = temp.freq(freq_sel);
                    data_tfr.(analysis){itrial}.(channame).powspctrm = temp.powspctrm(:,freq_sel,:);
                end
            end
            
            % gather all channels in the saame tfr structure
            %  tfr_all is a structure containg the clean powerspectrum for all channel in one trial
            temp =  struct2cell(data_tfr.(analysis){itrial});
            cfgtemp = [];
            cfgtemp.parameter  = 'powspctrm';
            tfr_all = ft_appendfreq(cfgtemp, temp{:});
            
            % selectionner les rgion du cerveau qui contiennent des data
            % Catego est l'ensemble des regions contenant les data
            str1=mat2str(cell2mat(config{1, irat}.LFP.catego));
            str2=split(str1(2:end-1));
            Catego = str2num(cell2mat(categories(categorical(str2))));
            
            %variable à rentrer
            leg = {};
            
            %% pour chaque région du cerveau
           
            %FIXME : passer à cette stratégie pour plus de clarté, enlever
            %tous les icat
            % for i_brainstruct = unique(string(config{irat}.LFP.structure))
            %selectionner les electrodes de la structure
%             channel_selected = string(config{irat}.LFP.structure) == i_brainstruct;
%             cfgtemp  = [];
%             cfgtemp.channel = config{irat}.LFP.channel(channel_selected);
%             lfp_structure = ft_selectdata(cfgtemp, LFP{1}.WoD_short);
%             tfr_structure = ft_selectdata(cfgtemp, tfr_all);
            
           for all_icat=1:length(Catego)
                icat = Catego(all_icat)
                for icat= icat
                    
                    % trier les analyses par classe de manière homogène
                    if config{irat}.classe.structure(icat).trial >1
                        %FIXME : 
                        %config{irat}.classe.HP
                        %config{irat}.classe.S1
                        %config{irat}.classe.TH
                        %config{irat}.classe.(i_brainstructure)
                        continue
                    end
                    
                    %% prendre le channel initiateur
                    %FIXME vérifier et retirer les 2 lignes ci dessous
%                     rongeur(irat).catego(icat).wod_time_ini;
%                     rongeur(irat).catego(icat).chan_depth;
                    
                    %verifier si le channel existe bien
                    channel_initiation   = rongeur(irat).catego(icat).chan_ini{1,1};
                    if isnan(cell2mat(channel_initiation))
                        continue
                    end
                    
                    %% analyses toutes les frequences de l'electrodes d'initiation d'une structure
                    
                    tfr_ini=data_tfr.(analysis){itrial}.(cell2mat(channel_initiation));
                    
                    %powerspectrum without the channel dimension
                    data_ini = permute(tfr_ini.powspctrm, [2 3 1]);
                    
                    %% plot toutes les frequences entre 0 et 50 hz du channel initiateur
                    if do_plot
                        fig1=figure;
                        plot(tfr_ini.time, data_ini,'LineWidth',1.5);
                        legend(string(round(tfr_ini.freq)));
                        xlabel ('time (s)','Fontsize',15);
                        ylabel ('power','Fontsize',15);
                        %title({['Initial WOD electrode powerspectrum'], ['in structure ' brain_struct{icat} ' at depth ' num2str( rongeur(irat).catego(icat).chan_depth) 'µm'], [config{irat}.prefix ' trial ' num2str(itrial)]},'Fontsize',12, 'interpreter', 'none');
                        title(sprintf('%s \nInitial WOD electrode power spectrum \nin structure %s at depth %dµm \n%s trial %d', analysis_names{itrial}, brain_struct{icat}, rongeur(irat).catego(icat).chan_depth(itrial), config{irat}.prefix, itrial),'Fontsize',12, 'interpreter', 'none');
                        figname = fullfile(savedir,'powspctrm',cell2mat(analysis_names(idata)),config{irat}.prefix, sprintf('Initial_powerspectrum_%s_trial%d', brain_struct{icat}, itrial));
                        dtx_savefigure(fig1, figname, 'png', 'pdf', 'close');%fonction de paul remplace la suite
                    end
                    
                    %% plot du globale power sur les channels moyennés d'une structure
                    
                    %avec la moyenne de toutes les électrodes de la structure :
                    tokeep = ~cellfun(@(c) any(isnan(c)), per_rat(irat).trial(itrial).struct(icat).chan_sel);
                    cfgtemp             = [];
                    cfgtemp.channel     = per_rat(irat).trial(itrial).struct(icat).chan_sel(tokeep);
                    cfgtemp.frequency   = 'all';
                    cfgtemp.avgoverfreq = 'yes';
                    cfgtemp.avgoverchan = 'yes';
                    cfgtemp.nanmean     = 'yes';
                    data.globale    = ft_selectdata(cfgtemp, tfr_all);
                    data.globale    = permute(data.globale.powspctrm, [2 3 1]);
                    
                                       
                    %find peak max of the global power its time
                    [power_detection(icat).peakglobale.value{irat}(itrial), power_detection(icat).peakglobale.time{irat}(itrial)] = ...
                        findpeaks(data.globale, tfr_all.time, 'NPeaks', 1, 'SortStr', 'descend');
                    
                    %FIXME
                    %save values in a structure
                    %power_detection(icat).globale = [power_detection(icat).globale  catego(icat).peakglobale.time{irat}(itrial)];
                    %power_detection(icat).globale(irat) =catego(icat).peakglobale.time{irat}(itrial);
                    
                    if do_plot
                        % plot le global power d'une structure
                        fig_globale=figure; hold on;
                        plot(tfr_all.time, data.globale,'LineWidth',1.5);
                        scatter(power_detection(icat).peakglobale.time{irat}(itrial), power_detection(icat).peakglobale.value{irat}(itrial), 'xr', 'linewidth', 2);
                        legend('globale power');
                        xlabel ('time (s)','Fontsize',15);
                        ylabel ('power','Fontsize',15);
                        axis tight
                        title(sprintf('%s \nGlobal powerspectrum averaged accross electrodes \nin structure %s  \n%s trial %d', analysis_names{itrial}, brain_struct{icat}, config{irat}.prefix, itrial),'Fontsize',12, 'interpreter', 'none');
                        figname =fullfile(savedir,'Global_powspctrm',cell2mat(analysis_names(idata)),config{irat}.prefix, sprintf('Initial_Range_powerspectrum_%s_trial%d', brain_struct{icat}, itrial));
                        dtx_savefigure(fig_globale, figname, 'png', 'pdf', 'close');
                    end
                    
                    
                    %definir TOI_detection = une fenetre de temps de detection  autour du pic du global power _ temps en sec
                    %si le peak est precoce
                    if power_detection(icat).peakglobale.time{irat}(itrial) < 10
                        toi_detection(1)=power_detection(icat).peakglobale.time{irat}(itrial) - 2;
                        toi_detection(2)= power_detection(icat).peakglobale.time{irat}(itrial) + 50;
                    else
                        toi_detection(1)= power_detection(icat).peakglobale.time{irat}(itrial) - 10;
                        toi_detection(2)= power_detection(icat).peakglobale.time{irat}(itrial) + 10;
                    end
                    
                    
                    
                    %% comparer l'evolution range de frequences au travers des différentes region
                    
                    %FIXME à clarifier ou à retirer. Renommer la variable
                    %où les données sont stockées pour la rendre explicite
                    for ifreqband = ["LF", "TF", "HF","VHF"]
                        
                        %power_detection(icat).(ifreqband)= []
                        
                        
                        
                        %avec uniquement l'électrode initiatrice de la strucutre
                        %cfgtemp             = [];
                        %cfgtemp.channel     = rongeur(irat).catego(icat).chan_ini{1,1};
                        %avec la moyenne de toutes les électrodes de la structure :
                        cfgtemp          = [];
                        tokeep           = ~cellfun(@(c) any(isnan(c)), per_rat(irat).trial(itrial).struct(icat).chan_sel);
                        cfgtemp.channel  = per_rat(irat).trial(itrial).struct(icat).chan_sel(tokeep);;
                        cfgtemp.frequency   = config{irat}.timefreq.(ifreqband);
                        cfgtemp.avgoverfreq = 'yes';
                        cfgtemp.avgoverchan = 'yes';
                        cfgtemp.nanmean     = 'yes';
                        cfgtemp.latency     = toi_detection;
                        data.(ifreqband)    = ft_selectdata(cfgtemp, tfr_all);
                        data.(ifreqband).powspctrm    = permute(data.(ifreqband).powspctrm, [2 3 1]);
                        
                        
                        %find peak max of the Range power and its time in a specific brain structure
                        [power_detection(icat).peakmean.(ifreqband).value{irat}(itrial), power_detection(icat).peakmean.(ifreqband).time{irat}(itrial)] = ...
                            findpeaks(data.(ifreqband).powspctrm, data.(ifreqband).time, 'NPeaks', 1, 'SortStr', 'descend');
                        
                        %FIXME
                        % save values in structure
                        %power_detection(icat).(ifreqband)=[power_detection(icat).(ifreqband)  catego(icat).peakmean.(ifreqband).time{irat}(itrial)];
                        
                    end %ifreqband
                    
                    %% subplot les range frequences du channel initiateur d'une région donnée dans une figure
                    if do_plot
                        
                        fig1=figure; hold on;
                        p = [];
                        iplot = 0;
                        for ifreqband = ["LF", "TF", "HF", "VHF"]
                            iplot = iplot+1;
                            subplot(2,2,iplot); hold on;
                            cfgtemp             = [];
                            cfgtemp.channel     = rongeur(irat).catego(icat).chan_ini{itrial};
                            cfgtemp.frequency   = config{irat}.timefreq.(ifreqband);
                            cfgtemp.avgoverfreq = 'yes';
                            cfgtemp.avgoverchan = 'yes'; %inutile puisque 1 seul channel
                            cfgtemp.nanmean     = 'yes'; %inutile puisque 1 seul channel
                            cfgtemp.latency     = toi_detection;
                            data_inirange.(ifreqband)    = ft_selectdata(cfgtemp, tfr_all);
                            data_inirange.(ifreqband).powspctrm    = permute(data_inirange.(ifreqband).powspctrm, [2 3 1]);
                            
                            %find peak max of the Range power and its time in the initiating channel of a specific brain structure
                            [power_detection(icat).peak_range_ini.(ifreqband).value{irat}(itrial), power_detection(icat).peak_range_ini.(ifreqband).time{irat}(itrial)] = ...
                                findpeaks(data_inirange.(ifreqband).powspctrm,data_inirange.(ifreqband).time , 'NPeaks', 1, 'SortStr', 'descend');
                            
                            cfgtemp.latency                      = [];
                            cfgtemp.latency                      = 'all';
                            data_inirange.(ifreqband)             = ft_selectdata(cfgtemp, tfr_all);
                            data_inirange.(ifreqband).powspctrm   = permute(data_inirange.(ifreqband).powspctrm, [2 3 1]);
                            
                            
                            %plot
                            p(end+1) = plot(tfr_all.time, data_inirange.(ifreqband).powspctrm ,'LineWidth',1.5);
                            scatter(power_detection(icat).peak_range_ini.(ifreqband).time{irat}(itrial), power_detection(icat).peak_range_ini.(ifreqband).value{irat}(itrial), 'xr', 'linewidth', 2);
                            title(sprintf('%s [%d-%d Hz]', ifreqband, config{irat}.timefreq.(ifreqband)(1), config{irat}.timefreq.(ifreqband)(2)));
                            xlabel ('time (s)','Fontsize',15);
                            ylabel ('power ()','Fontsize',15);
                            axis tight
                            xlim([0 180]);
                            
                        end
                        
                        sgtitle(sprintf('%s \nInitial WOD electrode Range powerspectrum \nin structure %s at depth %dµm \n%s trial %d', analysis_names{itrial}, brain_struct{icat}, rongeur(irat).catego(icat).chan_depth(itrial), config{irat}.prefix, itrial),'Fontsize',12, 'interpreter', 'none');
                        figname =fullfile(savedir,'Range_powspctrm',cell2mat(analysis_names(idata)),config{irat}.prefix, sprintf('Initial_Range_powerspectrum_%s_trial%d.jpg', brain_struct{icat}, itrial));
                        if ~isfolder(fullfile(savedir,'Range_powspctrm',cell2mat(analysis_names(idata)),config{irat}.prefix))
                            mkdir(fullfile(savedir,'Range_powspctrm',cell2mat(analysis_names(idata)),config{irat}.prefix));
                        end
                        exportgraphics(fig1,figname)
                        close all
                        
                        
                        %% plot range frequences moyenne de tous les channels d'une région donnée
                        fig1=figure; hold on;
                        p = [];
                        iplot = 0;
                        for ifreqband = ["LF", "TF", "HF", "VHF"]
                            
                            iplot = iplot+1;
                            subplot(2,2,iplot); hold on;
                            
                            cfgtemp             = [];
                            %moyenne de toutes les électrodes de la structure :
                            tokeep = ~cellfun(@(c) any(isnan(c)), per_rat(irat).trial(itrial).struct(icat).chan_sel);
                            cfgtemp.channel     = per_rat(irat).trial(itrial).struct(icat).chan_sel(tokeep);% rongeur(irat).catego(icat).chan_ini{1,1};
                            cfgtemp.frequency   = config{irat}.timefreq.(ifreqband);
                            cfgtemp.avgoverfreq = 'yes';
                            cfgtemp.avgoverchan = 'yes'; %inutile puisque 1 seul channel
                            cfgtemp.nanmean     = 'yes'; %inutile puisque 1 seul channel
                            cfgtemp.latency     = toi_detection;
                            data_mean.(ifreqband)    = ft_selectdata(cfgtemp, tfr_all);
                            data_mean.(ifreqband).powspctrm    = permute(data_mean.(ifreqband).powspctrm, [2 3 1]);
                            
                            y = data_mean.(ifreqband).powspctrm;
                            x = data_mean.(ifreqband).time;
                            
                            %                             [catego(icat).peak_mean.(ifreqband).value{irat}(itrial), catego(icat).peak_mean.(ifreqband).time{irat}(itrial)] = ...
                            %                                 findpeaks(y, x, 'NPeaks', 1, 'SortStr', 'descend');
                            
                            cfgtemp.latency = 'all';
                            data_mean.(ifreqband)    = ft_selectdata(cfgtemp, tfr_all);
                            data_mean.(ifreqband).powspctrm    = permute(data_mean.(ifreqband).powspctrm, [2 3 1]);
                            y = data_mean.(ifreqband).powspctrm;
                            x = data_mean.(ifreqband).time;
                            
                            p(end+1) = plot(x, y,'LineWidth',1.5);
                            scatter(power_detection(icat).peakmean.(ifreqband).time{irat}(itrial), power_detection(icat).peakmean.(ifreqband).value{irat}(itrial), 'xr', 'linewidth', 2);
                            
                            title(sprintf('%s [%d-%d Hz]', ifreqband, config{irat}.timefreq.(ifreqband)(1), config{irat}.timefreq.(ifreqband)(2)));
                            xlabel ('time (s)','Fontsize',15);
                            ylabel ('power ()','Fontsize',15);
                            axis tight
                            xlim([0 180]);
                            
                        end
                        
                        sgtitle(sprintf('%s \n WOD  Range powerspectrum \nin structure %s : \n%s trial %d', analysis_names{itrial}, brain_struct{icat},  config{irat}.prefix, itrial),'Fontsize',12, 'interpreter', 'none');
                        figname =fullfile(savedir,'Range_powspctrm',cell2mat(analysis_names(idata)),config{irat}.prefix, sprintf('Range_powerspectrum_%s_trial%d.jpg', brain_struct{icat}, itrial));
                        if ~isfolder(fullfile(savedir,'Range_powspctrm',cell2mat(analysis_names(idata)),config{irat}.prefix))
                            mkdir(fullfile(savedir,'Range_powspctrm',cell2mat(analysis_names(idata)),config{irat}.prefix));
                        end
                        exportgraphics(fig1,figname)
                        close all
                        
                        
                        %% Plot by range : TFR(mean over chan in structure) +  power spectrum (mean over chan in structure)+ detection
                        fig=figure('position',  get(0, 'ScreenSize'));
                        iplot = 0;
                        for ifreqband = ["LF", "TF", "HF", "VHF"]
                            iplot = iplot+1;
                            subplot(2,2,iplot); hold on;
                            %figure; hold on;
                            cfgtemp             = [];
                            %moyenne de toutes les électrodes de la structure :
                            tokeep = ~cellfun(@(c) any(isnan(c)), per_rat(irat).trial(itrial).struct(icat).chan_sel);
                            cfgtemp.channel     = per_rat(irat).trial(itrial).struct(icat).chan_sel(tokeep);% rongeur(irat).catego(icat).chan_ini{1,1};
                            cfgtemp.frequency   = config{irat}.timefreq.(ifreqband);
                            cfgtemp.avgoverchan = 'yes'; %inutile puisque 1 seul channel
                            cfgtemp.nanmean     = 'yes'; %inutile puisque 1 seul channel
                            dataplot    = ft_selectdata(cfgtemp, tfr_all);
                            %propriete TFR
                            cfgtemp              = [];
                            cfgtemp.colormap     = 'jet';
                            cfgtemp.masknans     = 'yes';
                            cfgtemp.xlim         = [0 180];
                            cfgtemp.ylim         = config{irat}.timefreq.(ifreqband);
                            ft_singleplotTFR(cfgtemp, dataplot);
                            
                            
                            %ylim([config{irat}.timefreq.(ifreqband)(1)+2, config{irat}.timefreq.(ifreqband)(2)]);
                        end
                        ft_pimpplot(fig, jet(5000))
                        % pour cacher le rajout superficiel en bas du graphique
                        for iplot = 2:4
                            subplot(2,2,iplot)
                            y = ylim;
                            ylim([y(1)+2, y(2)]);
                        end
                        
                        iplot = 0;
                        for ifreqband = ["LF", "TF", "HF", "VHF"]
                            
                            iplot = iplot+1;
                            s = subplot(2,2,iplot); hold on;
                            title(sprintf('%s [%d-%d Hz]', ifreqband, config{irat}.timefreq.(ifreqband)(1), config{irat}.timefreq.(ifreqband)(2)));
                            set(gca, 'tickdir', 'out', 'fontsize', 15);
                            %xline(catego(icat).peak_mean_ini.(ifreqband).time{irat}(itrial), 'k', 'linewidth', 2)
                            xline(power_detection(icat).peakmean.(ifreqband).time{irat}(itrial), 'k', 'linewidth', 2)
                            
                            ylabel('Freq (Hz)');
                            yyaxis right       % rajouter l'axe des power
                            
                            %moyenne de toutes les électrodes de la structure :
                            cfgtemp             = [];
                            tokeep = ~cellfun(@(c) any(isnan(c)), per_rat(irat).trial(itrial).struct(icat).chan_sel);
                            cfgtemp.channel     = per_rat(irat).trial(itrial).struct(icat).chan_sel(tokeep);% rongeur(irat).catego(icat).chan_ini{1,1};
                            cfgtemp.frequency   = config{irat}.timefreq.(ifreqband);
                            cfgtemp.avgoverfreq = 'yes';
                            cfgtemp.avgoverchan = 'yes'; %inutile puisque 1 seul channel
                            cfgtemp.nanmean     = 'yes'; %inutile puisque 1 seul channel
                            dataplot            = ft_selectdata(cfgtemp, tfr_all);
                            
                            %                             y = permute(data_inimean.(ifreqband).powspctrm, [3 2 1]);
                            %                             x = data_inimean.(ifreqband).time;
                            y = permute(dataplot .powspctrm, [3 2 1]);
                            x = dataplot.time;
                            
                            p(end+1) = plot(x, y, 'w', 'LineWidth',1);
                            
                            
                            ylabel('Mean power for the frequency band');
                            
                            xlabel('Time from Vent_Off', 'interpreter', 'none');
                        end
                        
                        sgtitle(sprintf('%s %s : trial %d', config{irat}.prefix, brain_struct{icat}, itrial), 'interpreter', 'none', 'fontsize', 22, 'fontweight', 'bold');
                        
                        figname =fullfile(savedir,'Range_TFRpowspctrm',cell2mat(analysis_names(idata)),config{irat}.prefix, sprintf('Range_TFRpowerspectrum_%s_trial%d.jpg', brain_struct{icat}, itrial));
                        if ~isfolder(fullfile(savedir,'Range_TFRpowspctrm',cell2mat(analysis_names(idata)),config{irat}.prefix))
                            mkdir(fullfile(savedir,'Range_TFRpowspctrm',cell2mat(analysis_names(idata)),config{irat}.prefix));
                        end
                        exportgraphics(fig,figname)
                        close all
                        
                        
                        
                        
                        %% plot range of frequenceies  des channel moyénné over brain  regions (1/2)
                        
                        fig2=figure(1); hold on;
                        leg{1}(icat) = plot(data.LF.time, data.LF.powspctrm,'LineWidth',1.5, 'color', brain_struct_color(icat, :));
                        scatter(power_detection(icat).peakmean.LF.time{irat}(itrial), power_detection(icat).peakmean.LF.value{irat}(itrial), 'xk', 'linewidth', 2);
                        
                        %fixme corriger les suivants :
                        fig3=figure(2); hold on;
                        leg{2}(icat) = plot(data.TF.time,data.TF.powspctrm,'LineWidth',1.5, 'color', brain_struct_color(icat, :));
                        scatter(power_detection(icat).peakmean.TF.time{irat}(itrial), power_detection(icat).peakmean.TF.value{irat}(itrial), 'xk', 'linewidth', 2);
                        
                        fig4=figure(3); hold on;
                        leg{3}(icat) = plot(data.HF.time, data.HF.powspctrm,'LineWidth',1.5, 'color', brain_struct_color(icat, :));
                        scatter(power_detection(icat).peakmean.HF.time{irat}(itrial), power_detection(icat).peakmean.HF.value{irat}(itrial), 'xk', 'linewidth', 2);
                        
                        fig5=figure(4); hold on;
                        leg{4}(icat) = plot(data.VHF.time, data.VHF.powspctrm,'LineWidth',1.5, 'color', brain_struct_color(icat, :));
                        scatter(power_detection(icat).peakmean.VHF.time{irat}(itrial), power_detection(icat).peakmean.VHF.value{irat}(itrial), 'xk', 'linewidth', 2);
                        
                        %FIXME stocker plutôt les valeurs comme ça :
                        power_detection.(i_brainstruct).peakmean.VHF.time{irat}{itrial} = 0
                        
                    end%do_plot
                    
                end %icat
            end%all_icat
            
            if do_plot
                %% plot range of frequenceies des channel moynnés over selected electrodes of a  brain  region (2/2)
                
                for ifigure=1:4
                    figure(ifigure)
                    xlabel ('time (s)','Fontsize',15);
                    ylabel ('power','Fontsize',15);
                    %legend(leg{ifigure}, brain_struct);
                    leg{ifigure}=legend( brain_struct);
                    axis tight
                    title(sprintf('%s \nPre-WOD Range powerspectrum \nover brain regions %s trial %d', analysis_names{itrial}, config{irat}.prefix, itrial),'Fontsize',12, 'interpreter', 'none');
                end
                
                if ~exist(fullfile(savedir,'Category_Range_powspctrm',analysis,config{irat}.prefix))
                    mkdir(fullfile(savedir,'Category_Range_powspctrm',analysis,config{irat}.prefix));
                end
                
                figname2 =fullfile(savedir,'Category_Range_powspctrm',analysis,config{irat}.prefix, sprintf('Category_LF_powerspectrum_trial%d.jpg', itrial));
                figname3 =fullfile(savedir,'Category_Range_powspctrm',analysis,config{irat}.prefix, sprintf('Category_TF_powerspectrum_trial%d.jpg', itrial));
                figname4 =fullfile(savedir,'Category_Range_powspctrm',analysis,config{irat}.prefix, sprintf('Category_HF_powerspectrum_trial%d.jpg', itrial));
                figname5 =fullfile(savedir,'Category_Range_powspctrm',analysis,config{irat}.prefix, sprintf('Category_VHF_powerspectrum_trial%d.jpg', itrial));
                
                dtx_savefigure(fig2, figname2, 'png', 'pdf', 'close');
                dtx_savefigure(fig3, figname3, 'png', 'pdf', 'close');
                dtx_savefigure(fig4, figname4, 'png', 'pdf', 'close');
                dtx_savefigure(fig5, figname5, 'png', 'pdf', 'close');
                
            end%do_plot
            
            %% plot pourune fréquence toute les region
            if do_plot
                for ifreqband = ["LF", "TF", "HF", "VHF"]
                    
                    fig=figure('position',  get(0, 'ScreenSize'));
                    iplot = 0;
                    
                    for all_icat=1:length(Catego)
                        icat = Catego(all_icat)
                        
                        iplot = iplot+1;
                        subplot(3,2,iplot); hold on;
                        %figure; hold on;
                        cfgtemp             = [];
                        %moyenne de toutes les électrodes de la structure :
                        tokeep = ~cellfun(@(c) any(isnan(c)), per_rat(irat).trial(itrial).struct(icat).chan_sel);
                        cfgtemp.channel     = per_rat(irat).trial(itrial).struct(icat).chan_sel(tokeep);% rongeur(irat).catego(icat).chan_ini{1,1};
                        cfgtemp.frequency   = config{irat}.timefreq.(ifreqband);
                        cfgtemp.avgoverchan = 'yes'; %inutile puisque 1 seul channel
                        cfgtemp.nanmean     = 'yes'; %inutile puisque 1 seul channel
                        dataplot    = ft_selectdata(cfgtemp, tfr_all);
                        %propriete TFR
                        cfgtemp              = [];
                        cfgtemp.colormap     = 'jet';
                        cfgtemp.masknans     = 'yes';
                        cfgtemp.xlim         = [0 180];
                        cfgtemp.ylim         = config{irat}.timefreq.(ifreqband);
                        ft_singleplotTFR(cfgtemp, dataplot);
                    end
                    ft_pimpplot(fig, jet(5000))
                    % pour cacher le rajout superficiel en bas du graphique
                    
                    for iplot = 2:all_icat
                        subplot(3,2,iplot)
                        y = ylim;
                        ylim([y(1)+2, y(2)]);
                    end
                    
                    iplot = 0;
                    for all_icat=1:length(Catego)
                        icat = Catego(all_icat)
                        
                        iplot = iplot+1;
                        s = subplot(3,2,iplot); hold on;
                        %title(sprintf('%s [%d-%d Hz]', ifreqband, config{irat}.timefreq.(ifreqband)(1), config{irat}.timefreq.(ifreqband)(2)));
                        title(sprintf('%s', brain_struct{icat}));
                        set(gca, 'tickdir', 'out', 'fontsize', 15);
                        xline(power_detection(icat).peakmean.(ifreqband).time{irat}(itrial), 'k', 'linewidth', 2)
                        ylabel(sprintf('Freq (Hz), [%d-%d Hz]', config{irat}.timefreq.(ifreqband)(1), config{irat}.timefreq.(ifreqband)(2)));
                        yyaxis right       % rajouter l'axe des power
                        
                        %moyenne de toutes les électrodes de la structure :
                        cfgtemp             = [];
                        tokeep = ~cellfun(@(c) any(isnan(c)), per_rat(irat).trial(itrial).struct(icat).chan_sel);
                        cfgtemp.channel     = per_rat(irat).trial(itrial).struct(icat).chan_sel(tokeep);
                        cfgtemp.frequency   = config{irat}.timefreq.(ifreqband);
                        cfgtemp.avgoverfreq = 'yes';
                        cfgtemp.avgoverchan = 'yes';
                        cfgtemp.nanmean     = 'yes';
                        dataplot            = ft_selectdata(cfgtemp, tfr_all);
                        
                        y = permute(dataplot .powspctrm, [3 2 1]);
                        x = dataplot.time;
                        p(end+1) = plot(x, y, 'w', 'LineWidth',1);
                        ylabel(sprintf('Mean power'))% for the frequency band %s %s' ,(ifreqband), brain_struct{icat}));
                        xlabel('Time from Vent_Off', 'interpreter', 'none');
                    end
                    
                    sgtitle(sprintf('%s %s  [%d-%d Hz]: trial %d', config{irat}.prefix, (ifreqband),  config{irat}.timefreq.(ifreqband)(1), config{irat}.timefreq.(ifreqband)(2),  itrial), 'interpreter', 'none', 'fontsize', 22, 'fontweight', 'bold');
                    
                    figname =fullfile(savedir,'Range_TFRpowspctrm_allstructures',cell2mat(analysis_names(idata)),config{irat}.prefix, sprintf('Range_TFRpowerspectrum_%s_trial%d.jpg', (ifreqband), itrial));
                    if ~isfolder(fullfile(savedir,'Range_TFRpowspctrm_allstructures',cell2mat(analysis_names(idata)),config{irat}.prefix))
                        mkdir(fullfile(savedir,'Range_TFRpowspctrm_allstructures',cell2mat(analysis_names(idata)),config{irat}.prefix));
                    end
                    exportgraphics(fig,figname)
                    close all
                end
                
                
            end
            
            
            
        end%itrial
    end%irat
    
    save(fname, 'power_detection');
    
    return
    
    

    
    
    
    
    
    %% OLD
    
    
    
    %% organiser les donner pour un plot statistic
    for icat=1:6%size(config.classe.structure,2)
        nLF(icat)=length(power_detection(icat).peak_LF);
        nTF(icat)=length(power_detection(icat).peak_TF);
        nHF(icat)=length(power_detection(icat).peak_HF);
        nVHF(icat)=length(power_detection(icat).peak_VHF)
    end
    N=max([max(nLF),max(nTF),max(nHF),max(nVHF)]);
    data=NaN(N,6) ;
    jcat=1;
    for icat=1:6%size(config.classe.structure,2)
        data(1:length(power_detection(icat).peak_LF),jcat)=power_detection(icat).peak_LF';
        data(1:length(power_detection(icat).peak_TF),jcat+1)=power_detection(icat).peak_TF';
        data(1:length(power_detection(icat).peak_HF),jcat+2)=power_detection(icat).peak_HF';
        data(1:length(power_detection(icat).peak_VHF),jcat+3)=power_detection(icat).peak_VHF';
        jcat=jcat+4;
    end
    
    
    
    %% plot les peak de range du power spectrum across brain regions
    
    sd=nanstd(data);
    m=nanmean(data);
    
    fig5=figure(5)
    x=data;
    n=length(data); xx=([1:24])';
    r=repmat(xx,1,n)';
    g=r(:)';
    positions = [1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24];
    h=boxplot(x,g, 'positions', positions);
    set(h,'linewidth',2)
    set(gca,'xtick',[mean(positions(1:4)) mean(positions(5:8)) mean(positions(9:12)) mean(positions(13:16)) mean(positions(17:20)) mean(positions(21:24))])
    set(gca,'xticklabel',brain_struct,'Fontsize',20)
    color = repmat({'c'; 'y' ;'r';'b'} ,6,1)'
    h = findobj(gca,'Tag','Box');
    for j=1:length(h)
        patch(get(h(j),'XData'),get(h(j),'YData'),cell2mat(color(j)),'FaceAlpha',0.5);
    end
    yt = get(gca, 'YTick');
    axis([xlim    0  ceil(max(yt)*1.2)])
    xt=positions
    text(xt(1)-0.3,m(1)+3500,{[num2str(m(1))],['\pm', num2str(sd(1))]})
    for p=2:length(xt)
        text(xt(p)-0.3,m(p)+2000,{[num2str(m(p))],['\pm', num2str(sd(p))]})
    end
    hold on
    lgd=legend({'LF','TF','HF','VHF'},'Fontsize',16)
    xlabel ('Structures')
    ylabel ('Peak delay from Vent-off','Fontsize',15)
    title({['Comparaison accross brain regions of the '], ['Power spectrum peak delay from Vent Off of 3 frenquency ranges N=22']},'Fontsize',15)
    figname =fullfile(savedir,'Statistics',['peak_power_spectrum']);
    dtx_savefigure(fig5, figname, 'png', 'pdf', 'jpeg', 'close');
    hold off
    
    %% statistics ANOVAN des variable frequce /structure
    
    groupeHF = repmat({'HF'} ,size(data,1),1)'
    groupeLF = repmat({'LF'} ,size(data,1),1)'
    groupeTF = repmat({'TF'} ,size(data,1),1)'
    groupeTF = repmat({'VHF'} ,size(data,1),1)'
    
    groupe_cat= repmat(brain_struct ,size(data,1)*4,1)
    groupe= repmat([groupeHF; groupeTF; groupeLF; groupeVHF] ,size(brain_struct,2),1)'
    
    Groupe_cat= reshape(groupe_cat,[],1)';
    Groupe=reshape(groupe,[],1)';
    Y=reshape(data,1,[]);
    
    [p,tbl,stats,terms] = anovan(Y,{Groupe,Groupe_cat},'model','interaction',...
        'varnames',{'Groupe','Groupe_cat'},'display','on')
    results = multcompare(stats,'Dimension',[1 2])
    
    %% Anova 1 facteur comparaison d'une frequence entre les structures
    groupe_cat= repmat(brain_struct ,size(data,1),1);
    Groupe=reshape(groupe_cat,[],1)';
    dataHF_stat =[];
    dataTF_stat =[];
    dataLF_stat =[];
    dataVHF_stat =[];
    i=1;j=2;k=3;
    while i<=size(data,2)
        dataHF_stat  = [dataHF_stat data(:,i)'];
        dataTF_stat  = [dataTF_stat data(:,j)'];
        dataLF_stat  = [dataLF_stat data(:,k)'];
        dataVHF_stat = [dataVHF_stat data(:,m)'];
        i=i+3;
        j=j+3;
        k=k+3;
        m=m+4;
    end
    [p_HF,tbl_HF,statsHF,terms] = anovan(dataHF_stat,{Groupe},'model','interaction',...
        'varnames',{'Groupe'},'display','on');
    resultsHF = multcompare(statsHF,'Dimension',[1]);
    [p_LF,tbl_LF,statsLF,terms] = anovan(dataLF_stat,{Groupe},'model','interaction',...
        'varnames',{'Groupe'},'display','on');
    resultsLF = multcompare(statsLF,'Dimension',[1]);
    [p_TF,tbl_TF,statsTF,terms] = anovan(dataTF_stat,{Groupe},'model','interaction',...
        'varnames',{'Groupe'},'display','on');
    resultsTF = multcompare(statsTF,'Dimension',[1]);
    [p_VHF,tbl_VHF,statsVHF,terms] = anovan(dataVHF_stat,{Groupe},'model','interaction',...
        'varnames',{'Groupe'},'display','on');
    resultsVHF = multcompare(statsVHF,'Dimension',[1]);
    
    
end %idata
