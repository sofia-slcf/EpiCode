
function [config,rongeur,per_rat]=wod_per_structure(config)


% load le gros tableau des données de détections
fname_out = fullfile(config{1}.datasavedir,'Detection', sprintf('wod_wavedetection_allrat.mat'));
load(fname_out, 'wod_wavedetection_allrat');
for irat=1:length(config) %pour chaque rat sauf le 4
%     if irat==4
%         continue
%     end
config{irat}.to_ignore = ft_getopt(config{irat}, 'to_ignore', false);
if config{irat}.to_ignore
    continue
end


    structure=str2num(cell2mat((categories(categorical(cell2mat(config{irat}.LFP.catego)))))); % ouput : les structure enregistrées dans les donées de tel rat
    for itrial=1:size(wod_wavedetection_allrat{1, irat}.anoxia,2)%config{irat}.trial%ne marche pas pour le rat 6
        for istruct = 1:6%size(structure,1)%pour chaque structure
            index_structure = find(cell2mat(config{irat}.LFP.catego) == istruct); % donne l'index des structures d'interets
            per_rat(irat).trial(itrial).struct(istruct).chan_sel        =[];
            per_rat(irat).trial(itrial).struct(istruct).depth_sel       =[];
            per_rat(irat).trial(itrial).struct(istruct).rename_sel      =[];
            per_rat(irat).trial(itrial).struct(istruct).silent_sel      =[];
            per_rat(irat).trial(itrial).struct(istruct).wod_time        =[];
            per_rat(irat).trial(itrial).struct(istruct).wod_value       =[];
            per_rat(irat).trial(itrial).struct(istruct).wod_slope_time  =[];
            per_rat(irat).trial(itrial).struct(istruct).wod_slope_value =[];
            per_rat(irat).trial(itrial).struct(istruct).wod_start_time  =[];
            per_rat(irat).trial(itrial).struct(istruct).wod_start_value =[];
            per_rat(irat).trial(itrial).struct(istruct).wod_half_width  =[];
            
            
            if ~isempty(index_structure)
                for iindex=1:length(index_structure)% selectionne les parametre  pour cette strucuture d'interet pour ce rat et pour ce trial
                    per_rat(irat).trial(itrial).struct(istruct).wod_time(iindex)          = wod_wavedetection_allrat{1,irat}.WoD.peak_time(index_structure(iindex),itrial);
                    per_rat(irat).trial(itrial).struct(istruct).wod_value(iindex)         = wod_wavedetection_allrat{1,irat}.WoD.peak_value(index_structure(iindex),itrial);
                    per_rat(irat).trial(itrial).struct(istruct).wod_slope_time(iindex)    = wod_wavedetection_allrat{1,irat}.WoD.min_slope_time(index_structure(iindex),itrial);
                    per_rat(irat).trial(itrial).struct(istruct).wod_slope_value(iindex)   = wod_wavedetection_allrat{1,irat}.WoD.min_slope_value(index_structure(iindex),itrial);
                    per_rat(irat).trial(itrial).struct(istruct).wod_start_time(iindex)    = wod_wavedetection_allrat{1,irat}.WoD.start_time(index_structure(iindex),itrial);
                    per_rat(irat).trial(itrial).struct(istruct).wod_start_value(iindex)   = wod_wavedetection_allrat{1,irat}.WoD.start_slope_value(index_structure(iindex),itrial);
                    per_rat(irat).trial(itrial).struct(istruct).wod_half_width(iindex)    = wod_wavedetection_allrat{1,irat}.WoD.half_width(index_structure(iindex),itrial);
                    per_rat(irat).trial(itrial).struct(istruct).chan_sel{iindex}          = config{irat}.LFP.channel{index_structure(iindex)};
                    per_rat(irat).trial(itrial).struct(istruct).depth_sel(iindex)         = config{irat}.LFP.chan_depth(index_structure(iindex));
                    per_rat(irat).trial(itrial).struct(istruct).rename_sel{iindex}        = config{irat}.LFP.rename{index_structure(iindex)};
                    per_rat(irat).trial(itrial).struct(istruct).silent_sel(iindex)        = wod_wavedetection_allrat{1, irat}.ISO(index_structure(iindex),itrial);
                    per_rat(irat).trial(itrial).struct(istruct).wod_speed_up(iindex)      = wod_wavedetection_allrat{1,irat}.WoD.speed_up(itrial);%(index_structure(iindex))
                    per_rat(irat).trial(itrial).struct(istruct).wod_speed_dn(iindex)      = wod_wavedetection_allrat{1,irat}.WoD.speed_dn(itrial);%(index_structure(iindex))
                end
                
                idx  = find(min( per_rat(irat).trial(itrial).struct(istruct).wod_time));%trouve la wod la plus petite et l'index correspondants
                
                if length(idx)>1
                    idx=idx(1);%arbitraire
                end
                %selectionne les parametre pour la wod la plus petite de la structure
                rongeur(irat).catego(istruct).wod_time_ini(:,itrial)          = per_rat(irat).trial(itrial).struct(istruct).wod_time(idx);
                rongeur(irat).catego(istruct).wod_value_ini(:,itrial)         = per_rat(irat).trial(itrial).struct(istruct).wod_value(idx);
                rongeur(irat).catego(istruct).wod_slope_time_ini(:,itrial)    = per_rat(irat).trial(itrial).struct(istruct).wod_slope_time(idx);
                rongeur(irat).catego(istruct).wod_slope_value_ini(:,itrial)   = per_rat(irat).trial(itrial).struct(istruct).wod_slope_value(idx);
                rongeur(irat).catego(istruct).wod_start_time_ini(:,itrial)    = per_rat(irat).trial(itrial).struct(istruct).wod_start_time(idx);
                rongeur(irat).catego(istruct).wod_start_value_ini(:,itrial)   = per_rat(irat).trial(itrial).struct(istruct).wod_start_value(idx);
                rongeur(irat).catego(istruct).wod_half_width_ini(:,itrial)    = per_rat(irat).trial(itrial).struct(istruct).wod_half_width(idx);
                rongeur(irat).catego(istruct).wod_start_value_ini(:,itrial)   = per_rat(irat).trial(itrial).struct(istruct).wod_start_value(idx);
                rongeur(irat).catego(istruct).wod_half_width_ini(:,itrial)    = per_rat(irat).trial(itrial).struct(istruct).wod_half_width(idx);
                rongeur(irat).catego(istruct).wod_speed_up(:,itrial)          = wod_wavedetection_allrat{1,irat}.WoD.speed_up(1);
                rongeur(irat).catego(istruct).wod_speed_dn(:,itrial)          = wod_wavedetection_allrat{1,irat}.WoD.speed_dn(1);
                rongeur(irat).catego(istruct).chan_ini{:,itrial}              = per_rat(irat).trial(itrial).struct(istruct).chan_sel(idx);
                rongeur(irat).catego(istruct).chan_depth(:,itrial)            = per_rat(irat).trial(itrial).struct(istruct).depth_sel(idx);
                rongeur(irat).catego(istruct).chan_rename{:,itrial}           = per_rat(irat).trial(itrial).struct(istruct).rename_sel(idx);
                rongeur(irat).catego(istruct).chan_silent(:,itrial)           = per_rat(irat).trial(itrial).struct(istruct).silent_sel(idx);
            else
                rongeur(irat).catego(istruct).wod_time_ini(:,itrial)          = nan;
                rongeur(irat).catego(istruct).wod_value_ini(:,itrial)         = nan;
                rongeur(irat).catego(istruct).wod_slope_time_ini(:,itrial)    = nan;
                rongeur(irat).catego(istruct).wod_slope_value_ini(:,itrial)   = nan;
                rongeur(irat).catego(istruct).wod_start_time_ini(:,itrial)    = nan;
                rongeur(irat).catego(istruct).wod_start_value_ini(:,itrial)   = nan;
                rongeur(irat).catego(istruct).wod_half_width_ini(:,itrial)    = nan;
                rongeur(irat).catego(istruct).wod_start_value_ini(:,itrial)   = nan;
                rongeur(irat).catego(istruct).wod_half_width_ini(:,itrial)    = nan;
                rongeur(irat).catego(istruct).wod_speed_up(:,itrial)          = nan;
                rongeur(irat).catego(istruct).wod_speed_dn(:,itrial)          = nan;
                rongeur(irat).catego(istruct).chan_ini{:,itrial}              = {nan};
                rongeur(irat).catego(istruct).chan_depth(:,itrial)            = nan;
                rongeur(irat).catego(istruct).chan_rename{:,itrial}           = {nan};
                rongeur(irat).catego(istruct).chan_silent(:,itrial)           = nan;
            end
            
            %rajouter une condition 'classe' dans les config : visera à trier les données
            if  rongeur(irat).catego(istruct).wod_time_ini(:,itrial)> nanstd( per_rat(irat).trial(itrial).struct(istruct).depth_sel) &  rongeur(irat).catego(istruct).wod_time_ini(:,itrial)<nanstd( per_rat(irat).trial(itrial).struct(istruct).depth_sel)
                config{irat}.classe.structure(istruct).trial(itrial)=2;
            else
                config{irat}.classe.structure(istruct).trial(itrial)=1;
            end
            
            
            
            
            
        end%istruct
    end%trial
end%irat

end