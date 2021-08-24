function freq_data= wod_iso_domfreq(cfg,freq_data,force)

%load data
fname_out = fullfile(cfg{4}.datasavedir,'freq_data', sprintf('freq_data_allrat.mat'));
freq_datadir=fullfile(cfg{4}.datasavedir,'freq_data', sprintf('dom_freq_allrat.mat'));

if exist(fname_out, 'file') && force == false && exist(freq_datadir, 'file') && force == false
    load(fname_out, 'freq_data');
    load(freq_datadir,'dom_freq_alldata');
    return
end



analysis_names={'timefreq_wod','timefreq_wod_timenorm','timefreq_wod_blcorrected','timefreq_wod_timenorm_blcorrected'};
freq_imagedir=fullfile(cfg{4}.imagesavedir,'freq_band_peaks','dominant_frequencies');
if ~isfolder(freq_imagedir)
    mkdir(freq_imagedir);
end

fname_freqdata = fullfile(cfg{4}.datasavedir,'freq_data', sprintf('freq_data_allrat.mat'));



for idata= 1: size(analysis_names,2)
    for irat= 1:size(cfg,2)
        
        if isempty(cfg{irat})
            continue
        end
        
        cfgcommon = cfg{4}; %same common config for all rats
        
        fprintf('Load timefreq data for rat %d/%d\n', irat,size(cfg,2));
        data_temp = load(fullfile(cfg{irat}.datasavedir,[cfg{irat}.prefix,analysis_names{idata},'.mat']));
        
        data_rat = [];
        count_trials                    = 0;
        
        for itrial = 1:size(data_temp.(analysis_names{idata}),2)
            
            count_trials = count_trials +1;
            wod_rat(count_trials) = count_trials; %store rat id for each wod (some rats have several wod)
            chan_list                                   = fieldnames(data_temp.(analysis_names{idata}){itrial});
            
            data_temp.(analysis_names{idata}){itrial}.Puff=[];
            
            data_rat= data_temp.(analysis_names{idata}){itrial};
            
            
            for ichan= 1: numel(chan_list)
                
                chan_name = chan_list{ichan};
                %replace empty channels by nans
                hasdata = find(~structfun(@isempty,data_rat));
                if isempty(hasdata)
                    continue
                end
                
                if isempty(data_rat.(chan_name))
                    data_rat.(chan_name)                       = data_rat.(chan_list{hasdata(1)});
                    data_rat.(chan_name).powspctrm             = ones(size(data_rat.(chan_name).powspctrm));
                end
                
                %% Detect Isoelectric time
               
                irat
                itrial
                chan_name
                analysis_names{idata}
                
                GlobalPower.(chan_name)=data_rat.(chan_name);
                GlobalPower.(chan_name).powspctrm=squeeze(data_rat.(chan_name).powspctrm);
                GlobalPower.(chan_name).powspctrm=mean(GlobalPower.(chan_name).powspctrm,1);
                
                %determine threshold to cross
                endbaseline_sample=floor(0.1*length(GlobalPower.(chan_name).powspctrm));
                baseline=nanmean(GlobalPower.(chan_name).powspctrm(1:endbaseline_sample));
                thr_iso=0.1*baseline;
                
                idx_signal= GlobalPower.(chan_name).powspctrm>thr_iso;
                diff_signal=diff(idx_signal);
                thr_cross=find(diff_signal);
                nanflag=sum(isnan(GlobalPower.(chan_name).powspctrm(thr_cross)));
                
                
                
                if isempty(thr_cross) || nanflag>=1
                    thr_iso=0.5*baseline;
                    idx_signal= GlobalPower.(chan_name).powspctrm>thr_iso;
                    diff_signal=diff(idx_signal);
                    thr_cross=find(diff_signal);
                end
                
                if length(thr_cross)>1
                    idx_ISO=diff(thr_cross)>=2;
                    ISO_time= GlobalPower.(chan_name).time(thr_cross(idx_ISO));
                    idx=  ISO_time> 0.1*max(GlobalPower.(chan_name).time) & ISO_time< 0.9*max(GlobalPower.(chan_name).time);
                    ISO_time=ISO_time(find(idx,1,'last'));
                else
                    idx_ISO=thr_cross;
                    ISO_time=(GlobalPower.(chan_name).time(idx_ISO));
                end
                
                if irat==24 
                    ISO_time= GlobalPower.(chan_name).time(thr_cross(end));
                
                end
                    
                
                if  GlobalPower.(chan_name).powspctrm(1:end)==1
                    freq_data{irat}.Iso_time(ichan,itrial)= nan;
                else
                    
                    fig=figure;
                    sgtitle(sprintf('%s_Rat_%i_Trial_%i',chan_name,irat,itrial));
                    plot(GlobalPower.(chan_name).time,GlobalPower.(chan_name).powspctrm)
                    xline(ISO_time);
                    yline(thr_iso);
                    fname_iso= fullfile(cfg{4}.imagesavedir_data{2},'detection',sprintf('%sWOD%i_ISO_time_%s',cfg{irat}.prefix,itrial,chan_name));
                    
                    dtx_savefigure(fig,fname_iso,'pdf','png','close');
                    
                    freq_data{irat}.(analysis_names{idata}).ISO(ichan,itrial)=  ISO_time;
                end
                
                clear idx_signal diff_signal thr_cross  thr_iso idx_ISO ISO_time idx baseline endbaseline_sample
                
                %% Detect dominant frequency and plot it vs time by channel
                
                data_allfreq= data_rat;
                data_allfreq.(chan_name).powspctrm=squeeze(data_allfreq.(chan_name).powspctrm);
                
                dom_pow=max(data_allfreq.(chan_name).powspctrm,[],1);
                idx_domfreq=data_allfreq.(chan_name).powspctrm==dom_pow;
                
                for icol=1:size(data_allfreq.(chan_name).powspctrm,2)
                    if isempty(find(idx_domfreq(:,icol),1))
                        dom_freq(icol)=nan;
                    else
                        dom_freq(icol)=find(idx_domfreq(:,icol),1);
                    end
                end %icol
                
                %replace 50 and 100 Hz by 0 bc in Iso 50Hz and harmonics
                %are dominant frequency
                
                idx=dom_freq== 50 | dom_freq==100;
                idx_around= dom_freq==49 | dom_freq==51;
                dom_freq(idx)=0;
                dom_freq(idx_around)=0;
                
                %smooth result
                dom_freq=movmean(dom_freq,10);
                
                
                fig_domfreq=figure;
                plot(data_allfreq.(chan_name).time,dom_freq)
                ylim([0 100]);
                
                fname_domfreq=fullfile(freq_imagedir,'by_channel',sprintf('%sWoD%i_%s_%s',cfg{irat}.prefix,itrial,chan_name,analysis_names{idata}));
                dtx_savefigure(fig_domfreq,fname_domfreq,'png','close');
                
                dom_freq_alldata{irat}.(analysis_names{idata}){itrial}.(chan_name).trial=dom_freq;
                clear dom_freq dom_pow
                
                
            end %ichan
            clear dom_freq dom_pow
            
        end %itrial
        clear dom_freq dom_pow
        
    end %irat
end %idata


%% Save variables
save(fname_out,'freq_data');
save(freq_datadir,'dom_freq_alldata');

end %function