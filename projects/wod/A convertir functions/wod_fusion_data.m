function ordered_data=wod_fusion_data(data,cfg,force)

%load data
fname_out = fullfile(cfg{4}.datasavedir,'Detection', sprintf('wod_wavedetection_allprobes.mat'));
if exist(fname_out, 'file') && force == false
    load(fname_out, 'stats');
    return
end


%% Arrange structure for 32 chan 

%load data for 16 and 32 chans

stats_concat.ratname = [];
depth_start = 0;
depth_end = 3200; %µm
depth_step = 100;

icol = 0;
for idepth = depth_start:depth_step:depth_end %step de 250 car les électrodes sont espacées de 250µm
    icol = icol+1;
    irow = 0;
    for irat = 1:size(data, 2)
        for itrial = 1:size(data{irat}.wod_peak_time, 2)
            irow = irow+1;
            sel = abs(data{irat}.Depth(:, itrial) - idepth) < depth_step/2;
            for ifield = string(fieldnames(data{irat}))'
                if strcmp(ifield, 'ratname')
                    stats_concat.ratname = [stats_concat.ratname; string(data{irat}.ratname)];
                    continue
                end
                if sum(sel) == 1
                    stats_concat.(ifield)(irow, icol) = data{irat}.(ifield)(sel, itrial);
                    
                elseif sum(sel) == 0
                    stats_concat.(ifield)(irow, icol) = nan;
                elseif sum(sel) > 1
                    error('it should have only one electrode for one deepness');
                end
            end
        end
    end
end
%interpolate missing data
for irat = 1:size(wodtimings_all, 1)
    startinterp = find(~isnan(wodtimings_all(irat,:)), 1, 'first');
    endinterp = find(~isnan(wodtimings_all(irat,:)), 1, 'last');
    wodtimings_all(irat,startinterp:endinterp) = fillmissing(wodtimings_all(irat,startinterp:endinterp), 'linear'); %ou pchip, ou spline
end