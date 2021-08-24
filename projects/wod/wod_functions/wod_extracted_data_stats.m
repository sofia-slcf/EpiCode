function Stats= wod_extracted_data_stats(calculated_data,cfg,force)


detect_stats_dir= fullfile(cfg{4}.statsavedir,'Waves_detections');
fname_out=fullfile(detect_stats_dir,'calculateddata_stats.mat');

if exist(fname_out, 'file') && force == false
    load(fname_out, 'Stats');
    return
end

%% Make statiscal analysis between trials

% origin time and depth
for iwave=["WoD" "WoR"]
    for iana= ["origin_time" "origin_depth"]
        count=0;
        for itime= ["peak_time" "min_slope_time" "start_time"]
            count= count+1;
            %separate trials as different arrays
            first_trial=calculated_data.(iwave).(iana).(itime)(:,1);
            second_trial=calculated_data.(iwave).(iana).(itime)(:,2);
            %t-test between trials
            p=signrank(first_trial,second_trial);
            p_val_trials.(iwave).(iana)(count,1)=p;
        end %itime
        %p values corrections
        [h, crit_p, adj_ci_cvrg, adj_p]=fdr_bh(p_val_trials.(iwave).(iana));
        adj_pval_trials.(iwave).(iana)=adj_p;
    end %iana
end %iwave
clear first_trial second_trial

%same operation for average propagation speed
for iwave=["WoD" "WoR"]
    for itime= ["peak_time" "min_slope_time" "start_time"]
        count=0;
        
        for isens=["up" "down"]
            
            count=count+1;
            %separate trials as different arrays
            first_trial=calculated_data.(iwave).speed.(itime).(isens)(:,1);
            second_trial=calculated_data.(iwave).speed.(itime).(isens)(:,2);
            %t-test between trials
            p=ranksum(first_trial,second_trial);
            p_val_trials.(iwave).speed.(itime).(isens)(count,1)=p;
            
        end %isens
    end %itime
end %iwave


%% Make statistical analysis between waves

for iana= ["origin_time" "origin_depth"]
    count=0;
    for itime= ["peak_time" "min_slope_time" "start_time"]
        count= count+1;
        wod_first_trial=calculated_data.WoD.(iana).(itime)(:,1);
        wod_second_trial=calculated_data.WoD.(iana).(itime)(:,2);
        wod_alltrials=vertcat(wod_first_trial,wod_second_trial);
        
        wor_first_trial=calculated_data.WoR.(iana).(itime)(:,1);
        wor_second_trial=calculated_data.WoR.(iana).(itime)(:,2);
        wor_alltrials=vertcat(wor_first_trial,wor_second_trial);
        
        p=signrank(wod_alltrials,wor_alltrials);
        p_val_waves.(iana)(count,1)=p;
        
    end %itime
    [h, crit_p, adj_ci_cvrg, adj_p]=fdr_bh(p_val_waves.(iana));
    adj_pval_waves.(iana)=adj_p;
end %iana

%% Make analysis between directions on speed

for iwave=["WoD" "WoR"]
    for itime= ["peak_time" "min_slope_time" "start_time"]

            %t-test between directions
            p=ranksum(calculated_data.(iwave).speed.(itime).up(:,1),calculated_data.(iwave).speed.(itime).down(:,1));
            p_val_speed.(iwave).speed.(itime)=p;
            
    end %itime
end %iwave

%% Make analysis between waves in same direction

for itime= ["peak_time" "min_slope_time" "start_time"]
    for isens=["up" "down"]
        %separate trials as different arrays
        first_trial=calculated_data.(iwave).speed.(itime).(isens)(:,1);
        second_trial=calculated_data.(iwave).speed.(itime).(isens)(:,2);
        %t-test between directions
        p=ranksum(calculated_data.WoD.speed.(itime).(isens)(:,1),calculated_data.WoR.speed.(itime).(isens)(:,1));
        p_val_speed_waves.(itime).(isens)=p;
    end %isens
end %itime



detect_stats_dir= fullfile(cfg{4}.statsavedir,'Waves_detections');

if ~isfolder(detect_stats_dir)
    mkdir(detect_stats_dir);
end

%store p-values in a single structure
Stats.pval_waves=p_val_waves;
Stats.adj_pval_waves=adj_pval_waves;

Stats.pval_speed= p_val_speed;
Stats.pval_speedwaves=p_val_speed_waves;

Stats.pval_trials=p_val_trials;
Stats.adj_pval_trials=adj_pval_trials;


save(fname_out,'Stats');

end %function