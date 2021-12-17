function [per_rat_normwodmin,per_rat_normdelta]=wod_time_norm(config,rongeur,per_rat)

%normalise les temps des wod

for irat=1:length(config)
    for itrial=1:size(per_rat(irat).trial,2)
        WODmin_cat=[];
        for icat=1:6
            if ~isempty(per_rat(irat).trial(itrial).struct(icat).wod_time)
                WODmin_cat(icat)=nanmin(per_rat(irat).trial(itrial).struct(icat).wod_time);
            else
                WODmin_cat(icat)= nan;
            end
            WODmin_abs=nanmin(WODmin_cat);
            WODmin_max=nanmax(WODmin_cat);
            animal(irat).trial(itrial).deltaWOD=WODmin_max-WODmin_abs;
            
            per_rat_normwodmin(irat).trial(itrial).struct(icat).wod_time = per_rat(irat).trial(itrial).struct(icat).wod_time-WODmin_abs;
            per_rat_normdelta(irat).trial(itrial).struct(icat).wod_time  = per_rat_normwodmin(irat).trial(itrial).struct(icat).wod_time/animal(irat).trial(itrial).deltaWOD;
        
            if isempty(per_rat_normwodmin(irat).trial(itrial).struct(icat).wod_time)
                per_rat_normwodmin(irat).trial(itrial).struct(icat).wod_time=nan;
            end
            
            if isempty(per_rat_normdelta(irat).trial(itrial).struct(icat).wod_time)
                per_rat_normdelta(irat).trial(itrial).struct(icat).wod_time=nan;
            end
                
        end
    end
end

