function ordered_data=wod_fusion_data_separated(data,force,cfg,oridepth)

%load data
fname_out = fullfile(cfg{4}.datasavedir,'Detection', sprintf('ordered_data_allprobes_%s.mat',string(oridepth)));
if exist(fname_out, 'file') && force == false
    load(fname_out, 'ordered_data');
    return
end


%% Arrange structure for 32 chan

%load data for 16 and 32 chans

depth_start = 10;
depth_end = 2210; %µm
depth_step = 100;

icol = 0;
for idepth = depth_start:depth_step:depth_end %step de 250 car les électrodes sont espacées de 250µm
    icol = icol+1;
    irow = 0;
    for irat = 1:size(data, 2)
        if isempty(data{irat})
            continue
        end
        for itrial = 1:size(data{irat}.WoD, 2)
            if data{irat}.Depth(3,itrial)==0
                continue
            end
            
            irow = irow+1;
            sel = abs(data{irat}.Depth(:, itrial) - idepth) < depth_step/2;
            %                     if strcmp(ifield, 'ratname')
            %                         ordered_data.ratname = [ordered_data.ratname; string(data{irat}.ratname)];
            %                         continue
            %                     end
            if sum(sel) == 1
                ordered_data.WoD(irow, icol) = data{irat}.WoD.peak_time(sel, itrial);
                ordered_data.Depth(irow,icol)= data{irat}.Depth(sel,itrial);
            elseif sum(sel) == 0
                ordered_data.WoD(irow, icol) = nan;
                ordered_data.Depth(irow,icol)= nan;
                
            elseif sum(sel) > 1
                error('it should have only one electrode for one deepness');
            end
            
        end
    end
end

save(fname_out,'ordered_data');


