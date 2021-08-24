function wod_plot_freqdata(ordered_freqdata,cfg)

analysis_names={'timefreq_wod','timefreq_wod_timenorm','timefreq_wod_blcorrected','timefreq_wod_timenorm_blcorrected'};
freqdatapath=fullfile(cfg{4}.imagesavedir,'freq_band_peaks');




%% prepare data to plot

for idata=1:size(analysis_names,2)
    for iana=["peak_time" "peak_value"]
        for ifield = string(fieldnames(ordered_freqdata.(analysis_names{idata}).(iana)))'
            %interpolate missing values
            for irat = 1:size(ordered_freqdata.(analysis_names{idata}).(iana).(ifield), 1)
                startinterp = find(~isnan(ordered_freqdata.(analysis_names{idata}).(iana).(ifield)(irat,:)), 1, 'first');
                endinterp = find(~isnan(ordered_freqdata.(analysis_names{idata}).(iana).(ifield)(irat,:)), 1, 'last');
                ordered_freqdata_interp.(analysis_names{idata}).(iana).(ifield)(irat,startinterp:endinterp) = fillmissing(ordered_freqdata.(analysis_names{idata}).(iana).(ifield)(irat,startinterp:endinterp), 'linear'); %ou pchip, ou spline
                ordered_freqdata_interp.Depth(irat,startinterp:endinterp)= fillmissing(ordered_freqdata.Depth(irat,startinterp:endinterp), 'linear');
                ordered_freqdata_interp.Iso(irat,startinterp:endinterp)=fillmissing(ordered_freqdata.Iso_time(irat,startinterp:endinterp), 'linear');
            end
            
            %replace zeros by NaNs
            idx_zeros= find(ordered_freqdata_interp.(analysis_names{idata}).(iana).(ifield)(:,:)==0);
            idx_zeros_depth=find(ordered_freqdata_interp.Depth(:,:)==0);
            ordered_freqdata_interp.(analysis_names{idata}).(iana).(ifield)(idx_zeros)=nan;
            ordered_freqdata_interp.Depth(idx_zeros_depth)=nan;
            ordered_freqdata_interp.Iso(idx_zeros)=nan;
            
            
        end %ifield
        
        
        %% Make histogram with Kernel distribution
        
        bands= fieldnames(ordered_freqdata.(analysis_names{idata}).(iana))';
        Color_band=winter(size(bands,2));
        Color_ISO='k';
        fig_histo=figure;hold
        for ifield = string(fieldnames(ordered_freqdata.(analysis_names{idata}).(iana)))'
            
            color_idx=strcmp(ifield,bands);
            
            %put all timings in one column
            alldata= reshape(ordered_freqdata_interp.(analysis_names{idata}).(iana).(ifield),[],1);
            range_data=[min(alldata):1:max(alldata)];
            nbins= floor(size(range_data,2)/2);
            
            h=histfit(alldata,40,'kernel');
            set(h(2),'color',Color_band(color_idx,:));
            
            if iana=="peak_value" & idata==3
                xlim([-5 20]);
            end
            if iana=="peak_time"
            
                
                
                
            data_iso=reshape(ordered_freqdata_interp.Iso,[],1);
            range_data_iso=[min(data_iso):1:max(data_iso)];
            nbins= floor(size(range_data_iso,2)/2);
            
            h_iso=histfit(data_iso,40,'kernel');
            set(h(2),'color',Color_ISO);
            end
        end %ifield
        fname_histo= fullfile(freqdatapath,'histograms',sprintf('%s_%s',iana,analysis_names{idata}));
        
        if ~isfolder(fullfile(freqdatapath,'histograms'))
            mkdir(fullfile(freqdatapath,'histograms'));
        end
        
        dtx_savefigure(fig_histo,fname_histo,'pdf','png','close');
        
        
        
        %% Plot frequency peak data versus depth
        
        
        %% Mean +/- SEM
        fig_mean= figure; hold;
        for ifield = string(fieldnames(ordered_freqdata.(analysis_names{idata}).(iana)))'
            
            color_idx=strcmp(ifield,bands);
           
            
            %calculate median,mean,standard deviation and SEM
            for ichan= 1:size(ordered_freqdata.(analysis_names{idata}).(iana).(ifield),2)
                data_std(1,ichan)= std(ordered_freqdata_interp.(analysis_names{idata}).(iana).(ifield)(:,ichan),'omitnan');
                data_mean(1,ichan)=nanmean(ordered_freqdata_interp.(analysis_names{idata}).(iana).(ifield)(:,ichan));
                data_sem=data_std / (sqrt(length(ordered_freqdata_interp.(analysis_names{idata}).(iana).(ifield))));
                data_mean_depth(1,ichan)=nanmean(ordered_freqdata_interp.Depth(:,ichan));
                data_stdISO(1,ichan)=std(ordered_freqdata_interp.Iso(:,ichan),'omitnan');
                data_meanISO(1,ichan)=nanmean(ordered_freqdata_interp.Iso(:,ichan));
                data_semISO(1,ichan)=data_stdISO(1,ichan)/(sqrt(length(ordered_freqdata_interp.Iso(:,ichan))));
                
            end %ichan
            
            plot(data_mean,data_mean_depth,'Color',Color_band(color_idx,:),'LineWidth',1.5);
            patch([data_mean- data_sem, data_mean(end:-1:1)+ data_sem(end:-1:1)], [data_mean_depth, data_mean_depth(end:-1:1)], Color_band(color_idx,:), 'facealpha', 0.3, 'edgecolor', 'none');
            
            
            if iana=="peak_time"
            plot(data_meanISO,data_mean_depth,'Color','k','LineWidth',1.5);
            patch([data_meanISO- data_semISO, data_meanISO(end:-1:1)+ data_semISO(end:-1:1)], [data_mean_depth, data_mean_depth(end:-1:1)], Color_ISO, 'facealpha', 0.3, 'edgecolor', 'none');
            end
            
            if idata==2 || 4 && iana=="peak_time"
                xlim([0 1]);
            end
            
            
        end %ifield
        fname_mean=fullfile(freqdatapath,'freq_depth',sprintf('mean_sem_%s_%s',iana,analysis_names{idata}));
        
        if ~isfolder(fullfile(freqdatapath,'freq_depth'))
            mkdir(fullfile(freqdatapath,'freq_depth'));
        end
        
        dtx_savefigure(fig_mean,fname_mean,'png','pdf','close');
        
        %% Med +/- SD
        
        fig_med= figure; hold;
        for ifield = string(fieldnames(ordered_freqdata.(analysis_names{idata}).(iana)))'
            
            color_idx=strcmp(ifield,bands);
           for ichan= 1:size(ordered_freqdata.(analysis_names{idata}).(iana).(ifield),2)
                data_med(1,ichan)= nanmedian(ordered_freqdata_interp.(analysis_names{idata}).(iana).(ifield)(:,ichan));
                data_std(1,ichan)= std(ordered_freqdata_interp.(analysis_names{idata}).(iana).(ifield)(:,ichan),'omitnan');
                data_med_depth(1,ichan)=nanmedian(ordered_freqdata_interp.Depth(:,ichan));
                data_medISO(1,ichan)= nanmedian(ordered_freqdata_interp.Iso(:,ichan));
                data_stdISO(1,ichan)= std(ordered_freqdata_interp.Iso(:,ichan),'omitnan');
           end %ichan
            
            plot(data_med,data_med_depth,'Color',Color_band(color_idx,:),'LineWidth',1.5);
            patch([data_med- data_std, data_med(end:-1:1)+ data_std(end:-1:1)], [data_med_depth, data_med_depth(end:-1:1)], Color_band(color_idx,:), 'facealpha', 0.3, 'edgecolor', 'none');
            
            if iana=="peak_time"
                plot(data_medISO,data_med_depth,'Color',Color_band(color_idx,:),'LineWidth',1.5);
                patch([data_medISO- data_stdISO, data_medISO(end:-1:1)+ data_stdISO(end:-1:1)], [data_med_depth, data_med_depth(end:-1:1)], Color_band(color_idx,:), 'facealpha', 0.3, 'edgecolor', 'none');
            end
            
        end %ifield
        
        fname_med=fullfile(freqdatapath,'freq_depth',sprintf('med_sd_%s_%s',iana,analysis_names{idata}));
        
        if ~isfolder(fullfile(freqdatapath,'freq_depth'))
            mkdir(fullfile(freqdatapath,'freq_depth'));
        end
        
        dtx_savefigure(fig_med,fname_med,'png','pdf','close');
            
            
    end %iana
end %idata
end %function