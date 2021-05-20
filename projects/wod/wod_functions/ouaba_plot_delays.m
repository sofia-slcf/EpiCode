function wod_plot_delay(ordered_data,cfg)

        


%% load data

for iwave= ["WoD"]
    for itime= ["peak_time" "min_slope_time" "start_time"]
        %determine image paths according to data
        
        delayimagespath=fullfile(cfg{3}.imagesavedir,'delays');
        
        if cfg{3}.LFP.inject_depth<2000 && size(cfg,2)<8
            delayimagespath=fullfile(cfg{3}.imagesavedir,'delays','inject_900');
        elseif size(cfg,2) <8
            delayimagespath=fullfile(cfg{3}.imagesavedir,'delays','inject_2000');
        end

        fname_medraw=fullfile(delayimagespath,sprintf('%s_%s_median_raw',iwave,itime));
        fname_raw=fullfile(delayimagespath,sprintf('%s_%s_raw',iwave,itime));
        fname_meanraw=fullfile(delayimagespath,sprintf('%s_%s_mean_raw',iwave,itime));
        
        fname_norm=fullfile(delayimagespath,sprintf('%s_%s_norm',iwave,itime));
        fname_mednorm=fullfile(delayimagespath,sprintf('%s_%s_median_norm',iwave,itime));
        fname_meannorm=fullfile(delayimagespath,sprintf('%s_%s_mean_norm',iwave,itime));
        
        
        
        
        %interpolate missing values
        for irat = 1:size(ordered_data.(iwave).(itime), 1)
            startinterp = find(~isnan(ordered_data.(iwave).(itime)(irat,:)), 1, 'first');
            endinterp = find(~isnan(ordered_data.(iwave).(itime)(irat,:)), 1, 'last');
            ordered_data_interp.(iwave).(itime)(irat,startinterp:endinterp) = fillmissing(ordered_data.(iwave).(itime)(irat,startinterp:endinterp), 'linear'); %ou pchip, ou spline
            ordered_data_interp.Depth(irat,startinterp:endinterp)= fillmissing(ordered_data.Depth(irat,startinterp:endinterp), 'linear');
        end
        
        %replace zeros by NaNs
        idx_zeros= find(ordered_data_interp.(iwave).(itime)(:,:)==0);
        idx_zeros_depth=find(ordered_data_interp.Depth(:,:)==0);
        ordered_data_interp.(iwave).(itime)(idx_zeros)=nan;
        ordered_data_interp.Depth(idx_zeros_depth)=nan;
        
        %calculate median and standard deviation
        for ichan= 1:size(ordered_data.(iwave).(itime),2)
            data_med(1,ichan)= nanmedian(ordered_data_interp.(iwave).(itime)(:,ichan));
            data_std(1,ichan)= std(ordered_data_interp.(iwave).(itime)(:,ichan),'omitnan');
            data_med_depth(1,ichan)=nanmedian(ordered_data_interp.Depth(:,ichan));
            data_mean(1,ichan)=nanmean(ordered_data_interp.(iwave).(itime)(:,ichan));
            data_sem=data_std / (sqrt(length(ordered_data_interp.(iwave).(itime))));
            data_mean_depth(1,ichan)=nanmean(ordered_data_interp.Depth(:,ichan));

        end %ichan
        
        %% plot raw timings
        fig_raw=figure;
        sgtitle(sprintf('%s %s time',iwave,itime));

        
        C_prot={[0.5 0.5 0.5],[0.85 0.325 0.098]};
        C_med= {'k',[0.635 0.078 0.184]};
        
        for irat= 1:size(ordered_data.(iwave).(itime), 1)
            
            if iwave=="WoD"
                C=C_prot{1};
                C_bis=C_med{1};
            else
                C=C_prot{2};
                C_bis=C_med{2};
            end
            %plot each protocol
            plot(ordered_data_interp.(iwave).(itime)(irat,:),ordered_data_interp.Depth(irat,:),'Color',C,'LineWidth',1);
            hold on;
            
        end %irat
        %plot median on top
        plot(data_med,data_med_depth,'Color',C_bis,'LineWidth',1.5);
        ylim([0 2200]);
        xlim([0 max(max(ordered_data_interp.(iwave).(itime)(:,:)))+10]);

        if ~isfolder(delayimagespath)
            mkdir(delayimagespath);
        end
        %save figure
        dtx_savefigure(fig_raw,fname_raw,'png','pdf','close');
        
        %% Plot median +/- SD filled
        
        fig_medraw=figure;hold
        sgtitle(sprintf('%s %s time raw median',iwave,itime));

        
        
        C_prot={[0.5 0.5 0.5],[0.85 0.325 0.098]};
        C_med= {'k',[0.635 0.078 0.184]};
        C_std= {[0.5 0.5 0.5]; [1  0.549 0.412]};
        
        C=C_prot{1};
        C_bis=C_med{1};
        C_ter=C_std{1};
        
        
        plot(data_med,data_mean_depth,'Color',C_bis,'LineWidth',1.5);
        patch([data_med- data_std, data_med(end:-1:1)+ data_std(end:-1:1)], [data_med_depth, data_med_depth(end:-1:1)], C_ter, 'facealpha', 0.3, 'edgecolor', 'none');

        ylim([0 2200]);
        xlim([-20 max(data_med)+10]);

        if ~isfolder(delayimagespath)
            mkdir(delayimagespath);
        end
        %save figure
        dtx_savefigure(fig_medraw,fname_medraw,'png','pdf','close');
        
        %% Plot mean +/- SEM filled
        
        fig_meanraw=figure;hold
        sgtitle(sprintf('%s %s time raw mean',iwave,itime));

        
        
        C_prot={[0.5 0.5 0.5],[0.85 0.325 0.098]};
        C_med= {'k',[0.635 0.078 0.184]};
        C_std= {[0.5 0.5 0.5]; [1  0.549 0.412]};

            C=C_prot{1};
            C_bis=C_med{1};
            C_ter=C_std{1};

        
        plot(data_mean,data_mean_depth,'Color',C_bis,'LineWidth',1.5);
        patch([data_mean- data_sem, data_mean(end:-1:1)+ data_sem(end:-1:1)], [data_mean_depth, data_mean_depth(end:-1:1)], C_ter, 'facealpha', 0.3, 'edgecolor', 'none');

        ylim([0 2200]);
        xlim([0 max(data_mean)+10]);
        
        if ~isfolder(delayimagespath)
            mkdir(delayimagespath);
        end
        %save figure
        dtx_savefigure(fig_meanraw,fname_meanraw,'png','pdf','close');
        
        %% plot normalized timings
        
        %normalize data to origin
        for irat=1:size(ordered_data.(iwave).(itime), 1)
            data_plot_norm(irat,:)=ordered_data_interp.(iwave).(itime)(irat,:) -min(ordered_data_interp.(iwave).(itime)(irat,:));
        end %irat
        
        %calculate median and standard deviation
        for ichan= 1:size(ordered_data.(iwave).(itime),2)
            data_med_norm(1,ichan)= nanmedian(data_plot_norm(:,ichan));
            data_std_norm(1,ichan)= std(ordered_data_interp.(iwave).(itime)(:,ichan),'omitnan');
            data_mean_norm(1,ichan)=nanmean(data_plot_norm(:,ichan));
            data_sem_norm=data_std_norm / (sqrt(length(data_plot_norm)));
        end %ichan
        
        fig_norm=figure;
                sgtitle(sprintf('%s %s time normalized',iwave,itime));

        for irat=1:size(ordered_data.(iwave).(itime), 1)
            

                C=C_prot{1};
                C_bis=C_med{1};

            %plot each protocol
            plot(data_plot_norm(irat,:),ordered_data_interp.Depth(irat,:),'Color',C,'LineWidth',1);
            
            hold on;
            
        end %irat
        %plot median on top
        plot(data_med_norm,data_med_depth,'Color',C_bis,'LineWidth',1.5);
        ylim([0 2200]);
        xlim([-10 max(max(ordered_data.(iwave).(itime)(:,:)))+10]);
        
        
        if ~isfolder(delayimagespath)
            mkdir(delayimagespath);
        end
        %save figure
        dtx_savefigure(fig_norm,fname_norm,'png','pdf','close');
        
        
        %% Plot median +/- Std normalized
        
        fig_mednorm=figure;hold
        sgtitle(sprintf('%s %s time normalized median',iwave,itime));
        
        C_prot={[0.5 0.5 0.5],[0.85 0.325 0.098]};
        C_med= {'k',[0.635 0.078 0.184]};
        C_std= {'k'; 'r'};

            C=C_prot{1};
            C_bis=C_med{1};
            C_ter=C_std{1};


        plot(data_med_norm,data_med_depth,'Color',C_bis,'LineWidth',1.5);
%         plot((data_med_norm+(data_std_norm/2)),data_med_depth,'Color',C_ter,'LineWidth',0.5);
%         plot((data_med_norm-(data_std_norm/2)),data_med_depth,'Color',C_ter,'LineWidth',0.5);
%         
        patch([data_med_norm- data_std_norm, data_med_norm(end:-1:1)+ data_std_norm(end:-1:1)], [data_med_depth, data_med_depth(end:-1:1)], C_ter, 'facealpha', 0.3, 'edgecolor', 'none');

        ylim([0 2200]);
        xlim([-20 max(data_med_norm)+10]);

        if ~isfolder(delayimagespath)
            mkdir(delayimagespath);
        end
        %save figure
        dtx_savefigure(fig_mednorm,fname_mednorm,'png','pdf','close');
        
        
        %% Plot mean +/- Sem normalized
        
        fig_meannorm=figure;hold
        sgtitle(sprintf('%s %s time normalized mean',iwave,itime));
        
        C_prot={[0.5 0.5 0.5],[0.85 0.325 0.098]};
        C_med= {'k',[0.635 0.078 0.184]};
        C_std= {[0.5 0.5 0.5]; [1  0.549 0.412]};

            C=C_prot{1};
            C_bis=C_med{1};
            C_ter=C_std{1};

        plot(data_mean_norm,data_mean_depth,'Color',C_bis,'LineWidth',1.5);
        patch([data_mean_norm- data_sem_norm, data_mean_norm(end:-1:1)+ data_sem_norm(end:-1:1)], [data_mean_depth, data_mean_depth(end:-1:1)], C_ter, 'facealpha', 0.3, 'edgecolor', 'none');

        ylim([0 2200]);
        xlim([-10 max(data_mean_norm)+10]);

        if ~isfolder(delayimagespath)
            mkdir(delayimagespath);
        end
        %save figure
        dtx_savefigure(fig_meannorm,fname_meannorm,'png','pdf','close');

    end %itime
end %iwave




end %function
 