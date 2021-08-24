function wod_freq_data(ordered_freqdata,cfg)

freqdataimages=fullfile(cfg{4}.imagesavedir,'frequency_data');

if ~isfolder(freqdataimages)
    mkdir(freqdataimages);
end

analysis_names={'timefreq_wod_timenorm','timefreq_wod_timenorm_blcorrected'};

for idata= 1: size(analysis_names,2)
    
    for ivalue= ["peak_time"]
        
        fig_depth=figure;hold
        for iband= ["HF" "LF"]
            
            for ichan= 1:size(ordered_freqdata.(analysis_names{idata}).(ivalue).(iband),2)
                mean_freq(ichan,1)=nanmean(ordered_freqdata.(analysis_names{idata}).(ivalue).(iband)(:,ichan));
                std_freq(ichan,1)=std(ordered_freqdata.(analysis_names{idata}).(ivalue).(iband)(:,ichan),'omitnan');
                sem_freq(ichan,1)=std_freq(ichan,1)/size(ordered_freqdata.(analysis_names{idata}).(ivalue).(iband)(:,ichan),1);
                mean_depth(ichan,1)=nanmean(ordered_freqdata.Depth(:,ichan));
            end %ichan
            
            
            %interpolate data
            t= 1:size(mean_freq,1);
            t_interp= linspace(t(1),t(end),60);
            
            mean_freqinterp=pchip(t,mean_freq,t_interp);
            sem_freqinterp=pchip(t,sem_freq,t_interp);
            std_freqinterp=pchip(t,std_freq,t_interp);
            mean_depthinterp= pchip(t,mean_depth,t_interp);
            C_1={'k','r'};
            C_2={'k','r'};
            
            mean_depthinterp_flip=fliplr(mean_depthinterp);
            
            if iband== "HF"
                
                C_mean= C_1{1};
                C_sem= C_2{1};
                
            else
                
                C_mean= C_1{2};
                C_sem= C_2{2};
                
            end
            
            
            plot(mean_freqinterp,mean_depthinterp_flip,C_mean);
            patch([mean_freqinterp- sem_freqinterp, mean_freqinterp(end:-1:1)+ sem_freqinterp(end:-1:1)], [mean_depthinterp_flip, mean_depthinterp_flip(end:-1:1)], C_sem, 'facealpha', 0.3, 'edgecolor', 'none');
            xlim([0 1])
            
        end %iband
        fname_depth=fullfile(freqdataimages,sprintf('freq_data_depth_%s',analysis_names{idata}));
        dtx_savefigure(fig_depth,fname_depth,'png','pdf','close');
        
        fig_hist=figure;hold
        for iband=["HF" "LF"]
            if iband== "HF"
                
                C_mean= C_1{1};
                C_sem= C_2{1};
                
            else
                
                C_mean= C_1{2};
                C_sem= C_2{2};
                
            end
            A=reshape(ordered_freqdata.(analysis_names{idata}).(ivalue).(iband),[],1);
            histo= histfit(A,30,'weibull');
            set(histo(2),'color',C_mean);
            xlim([0 1]);
        end %iband
        
        fname_dist=fullfile(freqdataimages,sprintf('freq_data_distrib_%s',analysis_names{idata}));
        dtx_savefigure(fig_hist,fname_dist,'png','pdf','close');
        
        
    end %ivalue
end %idata

end %function
