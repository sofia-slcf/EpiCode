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

for irat= 1:size(cfg,2)
    
    if isempty(cfg{irat})
        continue
    end
    
    
    irat_name=sprintf('Rat_%i',irat);
        for itime= ["peak_time" "min_slope_time" "start_time"]
            for itrial= 1:size(stats_all{irat}.WoD.peak_time,2)
                %% WOD Find origin and store origin depth and timings
                
                %excluse trials without WoR
                if iwave== "WoR" & isnan(stats_all{irat}.WoD.(itime)(:,itrial))
                    continue
                end
                
                A=stats_all{irat}.WoD.(itime)(:,itrial);
                B= stats_all{irat}.Depth(:,itrial);
                %find minimum timing
                origin_time=min(A);
                %store minimum timing in structure
                calculated_data.WoD.origin_time.(itime)(irat,itrial)=origin_time;
                
                %find index of origin and get origin depth
                idx_origin=find(A==origin_time);
                origin_depth=stats_all{irat}.Depth(idx_origin,itrial);
                %store origin depth
                
                %security for 2 origin depth
                if size(origin_depth,1)>1
                    origin_depth=origin_depth(2);
                    A(idx_origin(1))==nan;
                    idx_origin=idx_origin(2);
                end
                
                calculated_data.WoD.origin_depth.(itime)(irat,itrial)=origin_depth;
                
                %% WOD Calculate propagation speed
                
                %Calculate propagation speed from origin in both ways
                for ichan= 1:size(A,1)
                    
                    if ichan== idx_origin
                        continue
                    end
                    
                    if ichan<idx_origin
                        Speed_up(ichan,1)= (B(idx_origin,1)-B(ichan,1))/(A(ichan,1)-A(idx_origin,1));
                    end
                    
                    
                    if ichan > idx_origin
                        Speed_down(ichan,1)= (B(ichan,1)-B(idx_origin,1))/(A(ichan,1)-A(idx_origin,1));
                    end
                end %ichan
                
                %replace zeros by NaN
                exist Speed_up
                bool_up= ans;
                exist Speed_down
                bool_down= ans;
                if bool_up==1
                    Nul_values_up= find(Speed_up(:,1)==0);
                    Speed_up(Nul_values_up,1)=NaN;
                end
                
                if bool_down==1
                    Nul_values_down= find(Speed_down(:,1)==0);
                    Speed_down(Nul_values_down,1)=NaN;
                end
                %Average and store values
                if bool_up==1
                    calculated_data.WoD.speed.(itime).up(irat,itrial)=nanmean(Speed_up);
                end
                if bool_down==1
                    calculated_data.WoD.speed.(itime).down(irat,itrial)=nanmean(Speed_down);
                end
                
                %% WOD Calculate Instantaneous speed
                
                %upwards
                for ichan= 1:idx_origin
                    if ichan== idx_origin
                        continue
                    end
                    Speed_instan(ichan,1)= (B(ichan+1,1)-B(ichan,1))/(A(ichan,1)-A(ichan+1,1));
                end %ichan
                
                %downwards
                for ichan= idx_origin:size(A,1)
                    if ichan==size(A,1)
                        continue
                    end
                    Speed_instan(ichan,1)= (B(ichan+1,1)-B(ichan,1))/(A(ichan+1,1)-A(ichan,1));
                end %ichan
                
                Speed_instant{irat}.WoD.(itime)(:,itrial)=Speed_instan;
                
                
                
                clear A origin_time origin_depth bool_up bool_down Speed_up Speed_down Speed_instan
                
            end %itrial
        end %itime
end %irat

%% Replace all 0 values by NaNs

%in non-speed data
for iwave=["WoD" "WoR"]
    for iana= ["origin_time" "origin_depth"]
        for itime= ["peak_time" "min_slope_time" "start_time"]
            idx=find(calculated_data.WoD.(iana).(itime)(:,:)==0);
            calculated_data.WoD.(iana).(itime)(idx)=NaN;
        end %itime
    end %iana
end %iwave

%take absolute values of speed and remove zeros
for iwave=["WoD" "WoR"]
    for itime= ["peak_time" "min_slope_time" "start_time"]
        for isens= ["down" "up"]
            idx=find(calculated_data.WoD.speed.(itime).(isens)(:,:)==0);
            calculated_data.WoD.speed.(itime).(isens)= abs(calculated_data.WoD.speed.(itime).(isens));
            calculated_data.WoD.speed.(itime).(isens)(idx)=NaN;
        end %isens
    end %itime
end %iwave


save(fullfile(detectiondatapath,'calculated_data.mat'),'calculated_data');
save(fullfile(detectiondatapath,'instant_speed.mat'),'Speed_instant')




