


%% load data


temp= load(fullfile(cfg{4}.datasavedir,'Detection',sprintf('wod_wavedetection_allprobes.mat')));
ordered_data= temp.ordered_data;




for iwave= ["WoD" "WoR"]
    for itime= ["peak_time" "min_slope_time" "start_time"]
        
        
        %interpolate missing values
        for irat = 1:size(ordered_data.(iwave).(itime), 1)
            startinterp = find(~isnan(ordered_data.(iwave).(itime)(irat,:)), 1, 'first');
            endinterp = find(~isnan(ordered_data.(iwave).(itime)(irat,:)), 1, 'last');
            ordered_data_interp.(iwave).(itime)(irat,startinterp:endinterp) = fillmissing(ordered_data.(iwave).(itime)(irat,startinterp:endinterp), 'pchip'); %ou pchip, ou spline
            %ordered_data_interp.Depth(irat,startinterp:endinterp)= fillmissing(ordered_data.Depth(irat,startinterp:endinterp), 'pchip');
        end
        
        %replace zeros by NaNs
        idx_zeros= find(ordered_data_interp.(iwave).(itime)(:,:)==0);
        ordered_data_interp.(iwave).(itime)(idx_zeros)=nan;
        
        %calculate median and standard deviation
        for ichan= 1:size(ordered_data.(iwave).(itime),2)
            data_med(1,ichan)= nanmedian(ordered_data_interp.(iwave).(itime)(:,ichan));
            data_std(1,ichan)= std(ordered_data_interp.(iwave).(itime)(:,ichan),'omitnan');
            data_med_depth(1,ichan)=nanmedian(ordered_data_interp.Depth(:,ichan));
        end %ichan
        
        %% plot raw timings
        fig_raw=figure;
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
        
        
        fname_raw=fullfile(cfg{4}.imagesavedir,'delays',sprintf('%s_%s_raw',iwave,itime));
        
        if ~isfolder(fullfile(cfg{4}.imagesavedir,'delays'))
            mkdir(fullfile(cfg{4}.imagesavedir,'delays'));
        end
        %save figure
        dtx_savefigure(fig_raw,fname_raw,'png','pdf','close');
        
        
        %% plot normalized timings
        
        %normalize data to origin
        for irat=1:size(ordered_data.(iwave).(itime), 1)
            data_plot_norm(irat,:)=ordered_data_interp.(iwave).(itime)(irat,:) -min(ordered_data_interp.(iwave).(itime)(irat,:)); 
        end %irat
        
        %calculate median and standard deviation
        for ichan= 1:size(ordered_data.(iwave).(itime),2)
            data_med_norm(1,ichan)= nanmedian(data_plot_norm(:,ichan));
            data_std_norm(1,ichan)= std(ordered_data_interp.(iwave).(itime)(:,ichan),'omitnan');
        end %ichan
        
        fig_norm=figure;
        
        for irat= 1:size(ordered_data.(iwave).(itime), 1)
            
            if iwave=="WoD"
                C=C_prot{1};
                C_bis=C_med{1};
            else
                C=C_prot{2};
                C_bis=C_med{2};
            end
            %plot each protocol
            %plot(data_plot_norm(irat,:),ordered_data_interp.Depth(irat,:),'Color',C,'LineWidth',1);
                        plot(data_plot_norm(irat,:),A,'Color',C,'LineWidth',1);

            hold on;
            
        end %irat
        %plot median on top
        %plot(data_med_norm,data_med_depth,'Color',C_bis,'LineWidth',1.5);
                plot(data_med_norm,A,'Color',C_bis,'LineWidth',1.5);

        
        fname_raw=fullfile(cfg{4}.imagesavedir,'delays',sprintf('%s_%s_raw',iwave,itime));
        
        if ~isfolder(fullfile(cfg{4}.imagesavedir,'delays'))
            mkdir(fullfile(cfg{4}.imagesavedir,'delays'));
        end
        %save figure
        dtx_savefigure(fig_raw,fname_raw,'png','pdf','close');
        
        
    end %itime
end %iwave


