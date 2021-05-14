function wod_plot_freqdata(ordered_freqdata,cfg)

analysis_names={'timefreq_wod','timefreq_wod_blcorrected'};
freqdatapath=fullfile(cfg{4}.imagesavedir,'freq_band_peaks');


%% prepare data to plot

for idata=1:size(analysis_names,2)
    for iana=["peak_time" "peak_value"]
        for ifield = string(fieldnames(ordered_freqdata.(analysis_names{idata}).(iana)))'
%interpolate missing values
        for irat = 1:size(ordered_data.(iwave).(itime), 1)
            startinterp = find(~isnan(ordered_freqdata.(analysis_names{idata}).(iana).(ifield)(irat,:)), 1, 'first');
            endinterp = find(~isnan(ordered_freqdata.(analysis_names{idata}).(iana).(ifield)(irat,:)), 1, 'last');
            ordered_freqdata_interp.(analysis_names{idata}).(iana).(ifield)(irat,startinterp:endinterp) = fillmissing(ordered_freqdata.(analysis_names{idata}).(iana).(ifield)(irat,startinterp:endinterp), 'linear'); %ou pchip, ou spline
            ordered_freqdata_interp.Depth(irat,startinterp:endinterp)= fillmissing(ordered_freqdata.Depth(irat,startinterp:endinterp), 'linear');
        end
        
        %replace zeros by NaNs
        idx_zeros= find(ordered_freqdata_interp.(iwave).(itime)(:,:)==0);
        idx_zeros_depth=find(ordered_freqdata_interp.Depth(:,:)==0);
        ordered_freqdata_interp.(iwave).(itime)(idx_zeros)=nan;
        ordered_freqdata_interp.Depth(idx_zeros_depth)=nan;
        
        %calculate median,mean,standard deviation and SEM
        for ichan= 1:size(ordered_freqdata.(analysis_names{idata}).(iana).(ifield),2)
            data_med(1,ichan)= nanmedian(ordered_freqdata_interp.(analysis_names{idata}).(iana).(ifield)(:,ichan));
            data_std(1,ichan)= std(ordered_freqdata_interp.(analysis_names{idata}).(iana).(ifield)(:,ichan),'omitnan');
            data_med_depth(1,ichan)=nanmedian(ordered_freqdata_interp.Depth(:,ichan));
            data_mean(1,ichan)=nanmean(ordered_freqdata_interp.(analysis_names{idata}).(iana).(ifield)(:,ichan));
            data_sem=data_std / (sqrt(length(ordered_freqdata_interp.(analysis_names{idata}).(iana).(ifield))));
            data_mean_depth(1,ichan)=nanmean(ordered_freqdata_interp.Depth(:,ichan));

        end %ichan
        
        
        
        
        
        
        
        
        end %ifield
    end %iana  
end %idata
