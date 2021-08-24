function calculated_data= wod_propag_analysis(cfg,force)


detectiondatapath= fullfile(cfg{4}.datasavedir,'Detection');


fname_out=fullfile(detectiondatapath,'calculated_data.mat');
if exist(fname_out, 'file') && force == false
    load(fname_out, 'calculated_data');
    return
end


%% loading data
temp= load(fullfile(detectiondatapath,'wod_wavedetection_allrat'));
stats_all=temp.stats_all;
clear temp

for irat = 1:size(cfg,2)
    
    if isempty(cfg{irat})
        continue
    end

    irat_name = sprintf('Rat_%i',irat);
    
    for iwave=["WoD" "WoR"]
        for itime= ["peak_time" "min_slope_time" "start_time"]
            for itrial= 1:size(stats_all{irat}.(iwave).peak_time,2)
                %% WOD Find origin and store origin depth and timings
                
                %excluse trials without WoR
                if strcmp(iwave, "WoR") && all(isnan(stats_all{irat}.(iwave).(itime)(:,itrial)))
                    continue
                end 
                
                [origin_time, idx_origin] = min(stats_all{irat}.(iwave).(itime)(:,itrial));
                if sum(origin_time == stats_all{irat}.(iwave).(itime)(:,itrial)) > 1
                    idx_origin = find(origin_time == stats_all{irat}.(iwave).(itime)(:,itrial), 1, 'last');
                end
                
                origin_depth = stats_all{irat}.Depth(idx_origin,itrial);
                calculated_data.(iwave).origin_time.(itime)(irat,itrial)    = origin_time;
                calculated_data.(iwave).origin_depth.(itime)(irat,itrial)   = origin_depth;
                
                
                if iwave== "WoD"
                
                %separate protocols by origin depth
                limits = [700, 1300, 1700];
                depth_lim{irat}(itrial) = 0;
                for i = limits
                    if origin_depth > i
                        depth_lim{irat}(itrial) = depth_lim{irat}(itrial) + 1;
                    end
                end
                
                if depth_lim{irat}(itrial)==1
                    delays_1000{irat}.WoD.peak_time(:,itrial)=stats_all{irat}.WoD.peak_time(:,itrial);
                    delays_1000{irat}.Depth(:,itrial)=stats_all{irat}.Depth(:,itrial);
                elseif depth_lim{irat}(itrial)==2
                    delays_1400{irat}.WoD.peak_time(:,itrial)=stats_all{irat}.WoD.peak_time(:,itrial);
                    delays_1400{irat}.Depth(:,itrial)=stats_all{irat}.Depth(:,itrial);
                elseif depth_lim{irat}(itrial)==0 || depth_lim{irat}(itrial)==3
                    delays_others{irat}.WoD.peak_time(:,itrial)=stats_all{irat}.WoD.peak_time(:,itrial);
                    delays_others{irat}.Depth(:,itrial)=stats_all{irat}.Depth(:,itrial);
                end
                end
                
                %% WOD Calculate propagation speed

                
                tdiff = stats_all{irat}.(iwave).(itime)(:, itrial) - stats_all{irat}.(iwave).(itime)(idx_origin, itrial);
                ddiff =stats_all{irat}.Depth(:, itrial) - stats_all{irat}.Depth(idx_origin, itrial);
                 
                speed = ddiff ./ tdiff;
                speed_up = mean(speed(1:idx_origin-1));
                speed_dn = mean(speed(idx_origin+1:end));
                
                calculated_data.(iwave).speed.(itime).up(irat,itrial)=speed_up;
                calculated_data.(iwave).speed.(itime).down(irat,itrial)=speed_dn;

            end %itrial
        end %itime
    end %iwave
end %irat

 

%% Replace all 0 values by NaNs

%in non-speed data
for iwave=["WoD" "WoR"]
    for iana= ["origin_time" "origin_depth"]
        for itime= ["peak_time" "min_slope_time" "start_time"]
            idx=calculated_data.(iwave).(iana).(itime)(:,:)==0;
            calculated_data.(iwave).(iana).(itime)(idx)=NaN;
        end %itime
    end %iana
end %iwave

%take absolute values of speed and remove zeros
for iwave=["WoD" "WoR"]
    for itime= ["peak_time" "min_slope_time" "start_time"]
        for isens= ["down" "up"]
            idx=calculated_data.(iwave).speed.(itime).(isens)(:,:)==0;
            calculated_data.(iwave).speed.(itime).(isens)= abs(calculated_data.(iwave).speed.(itime).(isens));
            calculated_data.(iwave).speed.(itime).(isens)(idx)=NaN;
        end %isens
    end %itime
end %iwave

for iwave=["WoD" "WoR"]
            idx=calculated_data.(iwave).speed.(itime).(isens)(:,:)==0;
            delays.(iwave).speed.(itime).(isens)= abs(calculated_data.(iwave).speed.(itime).(isens));
            calculated_data.(iwave).speed.(itime).(isens)(idx)=NaN;
end %iwave



save(fullfile(detectiondatapath,'calculated_data.mat'),'calculated_data');
save(fullfile(detectiondatapath,'delays_1000.mat'),'delays_1000');
save(fullfile(detectiondatapath,'delays_1400.mat'),'delays_1400');
save(fullfile(detectiondatapath,'delays_others.mat'),'delays_others');



