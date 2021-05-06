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
        if isempty(cfg{irat})
            continue
        end
        for itrial = 1:size(data{irat}.WoD.peak_time, 2)
            irow = irow+1;
            sel = abs(data{irat}.Depth(:, itrial) - idepth) < depth_step/2;
            for iwave=["WoD" "WoR"]
                for ifield = string(fieldnames(data{irat}.(iwave)))'
%                     if strcmp(ifield, 'ratname')
%                         ordered_data.ratname = [ordered_data.ratname; string(data{irat}.ratname)];
%                         continue
%                     end
                    if sum(sel) == 1
                        ordered_data.(iwave).(ifield)(irow, icol) = data{irat}.(iwave).(ifield)(sel, itrial);
                        ordered_data.Depth(irow,icol)= data{irat}.Depth(sel,itrial);
                    elseif sum(sel) == 0
                        ordered_data.(iwave).(ifield)(irow, icol) = nan;
                        ordered_data.Depth(irow,icol)= nan;

                    elseif sum(sel) > 1
                        error('it should have only one electrode for one deepness');
                    end
                    
                end
            end
        end
    end
end

save(fname_out,'ordered_data')


