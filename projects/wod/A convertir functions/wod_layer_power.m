function wod_layer_power(cfg,freq_data,force)

fname_out = fullfile(cfg{4}.datasavedir,'freq_data', sprintf('pooled_layer_power_allrat.mat'));
if exist(fname_out, 'file') && force == false
    load(fname_out, 'pool_pow');
    return
end

analysis_names={'timefreq_wod','timefreq_wod_blcorrected'};

boxplotpath=fullfile(cfg{4}.imagesavedir,'freq_band_peaks','boxplots');

if ~isfolder(boxplotpath)
    mkdir(boxplotpath);
end


%% Pool power of peak for layer 5 and 2/3
for idata=1:size(analysis_names,2)
    for irat= 1:size(cfg,2)
        
        
        if isempty(cfg{irat})
            continue
        end
        for itrial= 1:size(freq_data{irat}.(analysis_names{idata}).peak_value.HF,2)
            for iband= ["HF" "LF"]
                L23Idx= cell2mat(cfg{irat}.LFP.chan_depth)>cfg{irat}.LFP.layer_depth{1} & cell2mat(cfg{irat}.LFP.chan_depth)<cfg{irat}.LFP.layer_depth{2};
                L5Idx= cell2mat(cfg{irat}.LFP.chan_depth)>cfg{irat}.LFP.layer_depth{3} & cell2mat(cfg{irat}.LFP.chan_depth)<cfg{irat}.LFP.layer_depth{4};
                
                pool_pow.L23.(analysis_names{idata}).(iband)(irat,itrial)=nanmean(freq_data{irat}.(analysis_names{idata}).peak_value.(iband)(L23Idx,itrial));
                pool_pow.L5.(analysis_names{idata}).(iband)(irat,itrial)=nanmean(freq_data{irat}.(analysis_names{idata}).peak_value.(iband)(L5Idx,itrial));
                
            end %iband
        end %itrial
    end %irat
end %idata

%remove zeros and replace by NaNs
for idata=1:size(analysis_names,2)
    for ilayer=["L23" "L5"]
        for iband= ["HF" "LF"]
            idx_zeros=pool_pow.(ilayer).(analysis_names{idata}).(iband)==0;
            pool_pow.(ilayer).(analysis_names{idata}).(iband)(idx_zeros)=nan;
        end %iband
    end %ilayer
end %idata

%% Make table with data

powertable=table.empty;
irow=0;
for idata=1:size(analysis_names,2)
    for ilayer=["L23" "L5"]
        for iband=["HF" "LF"]
            for irat=1:size(pool_pow.(ilayer).(analysis_names{idata}).(iband),1)
                
                if isempty(cfg{irat})
                    continue
                end
                
                for itrial=1:size(pool_pow.(ilayer).(analysis_names{idata}).(iband),2)
                    irow=irow+1;
                    
                    powertable.power(irow)=pool_pow.(ilayer).(analysis_names{idata}).(iband)(irat,itrial);
                    powertable.normalisation(irow)=idata;
                    
                    if ilayer=="L23"
                        layerIdx=2;
                    else
                        layerIdx=5;
                    end
                    
                    if iband=="HF"
                        bandIdx=1;
                    else
                        bandIdx=2;
                    end
                    
                    powertable.Layer(irow)=layerIdx;
                    powertable.pow_band(irow)=bandIdx;
                    powertable.trial(irow)=itrial;
                end %itrial
            end %irat
        end %iband
    end %ilayer
end %idata

%% Make ANOVAs for factors

anovastatpath=fullfile(cfg{4}.statsavedir,'freq_data','anova');



%HF surge
mdl_power= fitlm(powertable, 'power ~ normalisation + Layer + trial+ pow_band+ normalisation:Layer + normalisation:trial + normalisation:pow_band +Layer:trial +Layer:pow_band +trial:pow_band');
stats_power           = anova(mdl_power,'component');
fname_power=fullfile(anovastatpath,'anova_power_layerpooled');

writetable(powertable,fname_power,'WriteRowNames',true,'FileType','spreadsheet');

%% Make 2 by 2 comparisons

%between layers
pval_powertable=table.empty;
%Compare peak values according between layers for HF

irow=0;
for idata=1:size(analysis_names,2)
    for iband= 1:2
        irow=irow+1;

sel= powertable.pow_band==iband & powertable.Layer==2 & powertable.normalisation==idata;
datachan1= powertable.power(sel,:);


sel_2=powertable.pow_band==iband & powertable.Layer==5 & powertable.normalisation==idata;
datachan2= powertable.power(sel_2,:);

    p= ranksum(datachan1,datachan2);
    
    pval_powertable.(analysis_names{idata})(irow)=p;
    pval_powertable.pow_band(irow)=iband;


% %correct all p-values
% [~, ~, ~, adj_p]=fdr_bh(p);
    end %iband
end %idata
fname_pvalpower=fullfile(cfg{4}.statsavedir,'freq_data','pval_power_layers');
writetable(pval_powertable,fname_pvalpower,'FileType','Spreadsheet')

clear pval_powertable p 


%between bands
pval_powertable=table.empty;
%Compare peak values according between layers for HF

irow=0;
for idata=1:size(analysis_names,2)
    for ilayer=[2 5]
        irow=irow+1;

sel= powertable.pow_band==1 & powertable.Layer==ilayer & powertable.normalisation==idata;
datachan1= powertable.power(sel,:);


sel_2=powertable.pow_band==2 & powertable.Layer==ilayer & powertable.normalisation==idata;
datachan2= powertable.power(sel_2,:);

    p= ranksum(datachan1,datachan2);
    
    pval_powertable.(analysis_names{idata})(irow)=p;
    pval_powertable.Layer(irow)=ilayer;


% %correct all p-values
% [~, ~, ~, adj_p]=fdr_bh(p);
    end %ilayer
end %idata
fname_pvalpower=fullfile(cfg{4}.statsavedir,'freq_data','pval_power_bands');
writetable(pval_powertable,fname_pvalpower,'FileType','Spreadsheet')

clear pval_powertable p 


%% Make categorical barplot 
sel= powertable.pow_band==1 & powertable.Layer==2 & powertable.normalisation==2;
sel_2=powertable.pow_band==1 & powertable.Layer==5 & powertable.normalisation==2;
pow_hf= powertable.power(sel,:);
pow_hf(:,2)=powertable.power(sel_2,:);

clear sel sel_2

sel=powertable.pow_band==2 & powertable.Layer==2 & powertable.normalisation==2;
sel_2=powertable.pow_band==2 & powertable.Layer==5 & powertable.normalisation==2;
pow_lf= powertable.power(sel,:);
pow_lf(:,2)=powertable.power(sel_2,:);

clear sel sel_2


x=[1:4];

val=[nanmean(pow_hf(:,1)) nanmean(pow_hf(:,1)); nanmean(pow_lf(:,1)) nanmean(pow_lf(:,2))];
errhigh=[std(pow_hf(:,1),'omitnan') std(pow_hf(:,2),'omitnan'); std(pow_lf(:,1),'omitnan') std(pow_lf(:,2),'omitnan')];
errlow=[std(pow_hf(:,1),'omitnan') std(pow_hf(:,2),'omitnan'); std(pow_lf(:,1),'omitnan') std(pow_lf(:,2),'omitnan')];

b=bar(val,'grouped');
hold on

[ngroups,nbars]=size(val);

groupwidth=min(0.8,nbars/(nbars+1.5));

for i=1:nbars
x=(1:ngroups)-groupwidth/2 + (2*i-1)*groupwidth/(2*nbars);
er=errorbar(x,val(:,i),errhigh(:,i),'k','linestyle','none');
er.CapSize=15;
end

ylim([0 20])
set(gca,'xticklabel',{'HF','LF'})%This line should replace the numbers with your categorical label

fig=gcf;
fname_boxplot_pow=fullfile(boxplotpath,'pow_by_band_by_layer');
dtx_savefigure(fig,fname_boxplot_pow,'pdf','png','close');

end