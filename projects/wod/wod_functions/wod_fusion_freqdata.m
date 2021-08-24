function ordered_freqdata=wod_fusion_freqdata(data,cfg,data_delays,force)

%load data
fname_out = fullfile(cfg{4}.datasavedir,'Detection', sprintf('freq_data_allprobes.mat'));
if exist(fname_out, 'file') && force == false
    load(fname_out, 'ordered_freqdata');
    return
end

analysis_names={'timefreq_wod','timefreq_wod_timenorm','timefreq_wod_blcorrected','timefreq_wod_timenorm_blcorrected'};

%% Arrange structure for 32 chan

%load data for 16 and 32 chans
depth_start = 10;
depth_end = 2210; %µm
depth_step = 100; % real value 100

icol = 0;
for idepth = depth_start:depth_step:depth_end %step de 250 car les électrodes sont espacées de 250µm
    icol = icol+1;
    irow = 0;
    for irat = 1:size(data, 2)
        if isempty(cfg{irat})
            continue
        end
        for itrial = 1:size(data{irat}.(analysis_names{1}).peak_time.HF, 2)
            irow = irow+1;
            sel = abs(data{irat}.Depth(:, itrial) - idepth) < depth_step/2;
            for idata=1:size(analysis_names,2)
                for iana= ["peak_time" "peak_value"]
                for ifield = string(fieldnames(data{irat}.(analysis_names{idata}).(iana)))'
%                     if strcmp(ifield, 'ratname')
%                         ordered_data.ratname = [ordered_data.ratname; string(data{irat}.ratname)];
%                         continue
%                     end
                    if sum(sel) == 1
                        ordered_freqdata.(analysis_names{idata}).(iana).(ifield)(irow, icol) = data{irat}.(analysis_names{idata}).(iana).(ifield)(sel, itrial);
                        ordered_freqdata.Depth(irow,icol)= data{irat}.Depth(sel,itrial);
                        ordered_freqdata.Iso_time(irow,icol)= data_delays{irat}.ISO(sel,itrial);
                    elseif sum(sel) == 0
                        ordered_freqdata.(analysis_names{idata}).(iana).(ifield)(irow, icol) = nan;
                        ordered_freqdata.Depth(irow,icol)= nan;
                        ordered_freqdata.Iso_time(irow,icol)=nan;
                    elseif sum(sel) > 1
                        error('it should have only one electrode for one deepness');
                    end
                    
                end %ifield
                end %iana
            end %idata
        end
    end
end

save(fname_out,'ordered_freqdata')


