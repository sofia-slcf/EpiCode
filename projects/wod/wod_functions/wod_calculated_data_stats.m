function wod_calculated_data_stats(cfg,data)



anovastatpath=fullfile(cfg{4}.statsavedir,'Waves_detection','anova');

if ~isfolder(anovastatpath)
    mkdir(anovastatpath);
end

boxplotpath=fullfile(cfg{4}.imagesavedir,'delays','boxplots');

if ~isfolder(boxplotpath)
    mkdir(boxplotpath);
end



%% Make ANOVA for separated depths and trial number for WoD

irow=0;
calculated_table=table.empty;
for irat= 1:size(data,2)
    if isempty(data{irat})
        continue
    end
    
    for itrial=1:size(data{irat}.WoD.peak_time,2)
        irow=irow+1;
        
        calculated_table.origin_time(irow)      =data{irat}.wod_origin_time(itrial);
        calculated_table.speed_up(irow)         =data{irat}.WoD.speed_up(itrial);
        calculated_table.speed_dn(irow)         =data{irat}.WoD.speed_dn(itrial);
        calculated_table.trial(irow)            =itrial;
        calculated_table.depthclass(irow)       =data{irat}.oridepthclass(itrial);
        calculated_table.oridepth(irow)         =data{irat}.wod_origin_depth(itrial);
        calculated_table.anoxia(irow)           =data{irat}.anoxia(itrial);
        calculated_table.oridepthwor(irow)      =data{irat}.wor_origin_depth(itrial);
        calculated_table.origin_time_wor(irow)  =data{irat}.wor_origin_time(itrial);
        calculated_table.prewod_anoxia(irow)    =data{irat}.before_wod(itrial);
    end %itrial
end %irat

nan_flag=find(isnan(calculated_table.speed_up));
calculated_table([nan_flag],:)=[];
nan_flag2=find(isnan(calculated_table.speed_dn));
calculated_table([nan_flag2],:)=[];

%origin time
mdl_origintime= fitlm(calculated_table, 'origin_time ~ trial + depthclass + depthclass:trial');
stats_origintime           = anova(mdl_origintime,'component');
fname_origintime=fullfile(anovastatpath,'anova_origin_time');

%origin time of WoR
mdl_origintimewor= fitlm(calculated_table, 'origin_time_wor ~ trial + depthclass + depthclass:trial');
stats_origintimewor           = anova(mdl_origintimewor,'component');
fname_origintime_wor=fullfile(anovastatpath,'anova_origin_time_wor');

%upward speed
mdl_speedup= fitlm(calculated_table, 'speed_up ~ trial + depthclass + depthclass:trial');
stats_speedup          = anova(mdl_speedup,'component');
fname_speedup=fullfile(anovastatpath,'anova_speedup');

%downward speed
mdl_speeddn= fitlm(calculated_table, 'speed_dn ~ trial + depthclass + depthclass:trial');
stats_speeddn          = anova(mdl_speeddn,'component');
fname_speeddn=fullfile(anovastatpath,'anova_speeddn');

%depth of origin
mdl_oridepth= fitlm(calculated_table, 'oridepth ~ trial');
stats_oridepth          = anova(mdl_oridepth,'component');
fname_oridepth=fullfile(anovastatpath,'anova_oridepth');

%depth of origin of WoR
mdl_oridepthwor= fitlm(calculated_table, 'oridepthwor ~ trial + depthclass + depthclass:trial');
stats_oridepthwor          = anova(mdl_oridepthwor,'component');
fname_oridepthwor=fullfile(anovastatpath,'anova_oridepth_wor');

%anoxia duration
mdl_anoxia= fitlm(calculated_table, 'anoxia ~ trial+ depthclass + depthclass:trial');
stats_anoxia          = anova(mdl_oridepth,'component');
fname_anoxia=fullfile(anovastatpath,'anova_anoxia');



%% Save stats tables
writetable(stats_origintime,fname_origintime,'WriteRowNames',true,'FileType','spreadsheet');
writetable(stats_speedup,fname_speedup,'WriteRowNames',true,'FileType','spreadsheet');
writetable(stats_speeddn,fname_speeddn,'WriteRowNames',true,'FileType','spreadsheet');
writetable(stats_oridepth,fname_oridepth,'WriteRowNames',true,'FileType','spreadsheet');
writetable(stats_anoxia,fname_anoxia,'WriteRowNames',true,'FileType','spreadsheet');
writetable(stats_origintimewor,fname_origintime_wor,'WriteRowNames',true,'FileType','spreadsheet');
writetable(stats_oridepthwor,fname_oridepthwor,'WriteRowNames',true,'FileType','spreadsheet');


%% Make t-test comparison between groups and trials

statstable=table.empty;

for idata=["origin_time" "speed_up" "speed_dn" "oridepth" "anoxia" "origin_time_wor" "oridepthwor" "prewod_anoxia"]
    
    %between groups
    sel=calculated_table.depthclass==1000;
    data_1000=calculated_table.(idata)(sel,:);
    
    sel_2=calculated_table.depthclass==1400;
    data_1400=calculated_table.(idata)(sel_2,:);
    
    p=ranksum(data_1000,data_1400);
    clear sel sel_2
    
    %between trials
    sel=calculated_table.trial==1;
    data_t1=calculated_table.(idata)(sel,:);
    
    sel_2=calculated_table.trial==2;
    data_t2=calculated_table.(idata)(sel_2,:);
    
    p2=ranksum(data_t1,data_t2);
    clear sel sel_2
    
    %store in table
    statstable.(idata)=[p p2]';
    
    
end %idata

RowNames={'groups','trials'}';
fname_pval= fullfile(anovastatpath,'pval_calculated_data');
statstable.vs=RowNames;
writetable(statstable,fname_pval,'FileType','spreadsheet');
clear statstable

%% Make boxplots of comparisons


%for origin time and depth
for idata= ["origin_time" "oridepth" "speed_up" "speed_dn" "anoxia" "origin_time_wor" "oridepthwor"]
    
    %separate by depth class
    boxplot(calculated_table.(idata),calculated_table.depthclass,'Colors','k','Symbol','o');
    
    fig=gcf;
    fname_boxplot_depth=fullfile(boxplotpath,sprintf('%s_per_depth',idata));
    dtx_savefigure(fig,fname_boxplot_depth,'pdf','png','close');
    
    %separate by trial
    boxplot(calculated_table.(idata),calculated_table.trial,'Colors','k','Symbol','o');
    
    fig=gcf;
    fname_boxplot_trial=fullfile(boxplotpath,sprintf('%s_per_trial',idata));
    dtx_savefigure(fig,fname_boxplot_trial,'pdf','png','close');
    
end %idata

%grouped boxplot for speed

% Create example data

sel_1= calculated_table.depthclass==1000;
speed_1000= abs(calculated_table.speed_up(sel_1,:));
speed_1000(:,2)=abs(calculated_table.speed_dn(sel_1,:));

sel_2= calculated_table.depthclass==1400;
speed_1400= abs(calculated_table.speed_up(sel_2,:));
speed_1400(:,2)=abs(calculated_table.speed_dn(sel_2,:));

% prepare data
speed_cell=cell(2,2);
for ii=1:size(speed_cell,1)
    speed_1000_c{ii}=speed_1000(:,ii);
    speed_1400_c{ii}=speed_1400(:,ii);
end
speed_cell=vertcat(speed_1000_c,speed_1400_c);

xlab={'Upward speed','Downward speed'};
col=[255,64,64, 200;
    255,127,36, 200];
col=col/255;

multiple_boxplot(speed_cell',xlab,{'Superficial origin', 'Deep origin'},col')

fig=gcf;
fname_boxplot_trial=fullfile(boxplotpath,'speed_per_way_per group');
dtx_savefigure(fig,fname_boxplot_trial,'pdf','png','close');

%Categorical scatter plot for origin depth

fig_scatter=figure;
data_depth= calculated_table.oridepth;
[row, col]=size(data_depth);

xdata = repmat(1:col, row, 1);
scatter(xdata,data_depth,'jitter','on','jitterAmount',0.1);

ylim([0 2200]);

fname_scatter=fullfile(boxplotpath,'oridepth_scatter');
dtx_savefigure(fig_scatter,fname_scatter,'pdf','png','close');

clear row col xdata data_depth

%% Make table with all measured data from origin channels

%make table with all comparison variables
statstable = table.empty;
irow=0;
for irat=1:size(cfg,2)
    if isempty(cfg{irat})
        continue
    end

    for itrial=1:size(data{irat}.WoD.peak_time,2)
        irow= irow+1;
        sel=data{irat}.Depth(:,itrial)==data{irat}.wod_origin_depth(itrial);
        
        statstable.peakval(irow)=data{irat}.WoD.peak_value(sel,itrial);
        statstable.minslope(irow)=data{irat}.WoD.min_slope_value(sel,itrial);
        
        statstable.halfwidth(irow)=data{irat}.WoD.half_width(sel,itrial);
        statstable.Wod_delay(irow)=data{irat}.WoD.peak_time(sel,itrial);
        statstable.oriclass(irow)=data{irat}.oridepthclass(itrial);
        statstable.trial(irow)=itrial;
        statstable.chandepth(irow)=data{irat}.wod_origin_depth(itrial);
    end %itrial
end %irat


%% Make ANOVA of parameters

%peak_value
mdl_peakval= fitlm(statstable, 'peakval ~ chandepth + trial + oriclass + chandepth:trial + chandepth:oriclass + trial:oriclass');
stats_peakval          = anova(mdl_peakval,'component');
fname_peakval=fullfile(anovastatpath,'anova_peakval_origin');

%min_slope
mdl_minslope= fitlm(statstable, 'minslope ~ chandepth + trial + oriclass + chandepth:trial + chandepth:oriclass + trial:oriclass');
stats_minslope          = anova(mdl_minslope,'component');
fname_minslope=fullfile(anovastatpath,'anova_minslope_origin');

%half-width
mdl_halfwidth= fitlm(statstable, 'halfwidth ~ chandepth + trial + oriclass + chandepth:trial + chandepth:oriclass + trial:oriclass');
stats_halfwidth          = anova(mdl_halfwidth,'component');
fname_halfwidth=fullfile(anovastatpath,'anova_halfwidth_origin');

%% Save stats tables
writetable(stats_peakval,fname_peakval,'WriteRowNames',true,'FileType','spreadsheet');
writetable(stats_minslope,fname_minslope,'WriteRowNames',true,'FileType','spreadsheet');
writetable(stats_halfwidth,fname_halfwidth,'WriteRowNames',true,'FileType','spreadsheet');

%% Make 2 by 2 comparisons

pvaltable=table.empty;

for idata=["peakval" "minslope" "halfwidth"]
    %between groups
    sel=statstable.oriclass==1000;
    data_1000=statstable.(idata)(sel,:);
    
    sel_2=statstable.oriclass==1400;
    data_1400=statstable.(idata)(sel_2,:);
    
    p=ranksum(data_1000,data_1400);
    clear sel sel_2
    
    %between trials
    sel=statstable.trial==1;
    data_t1=statstable.(idata)(sel,:);
    
    sel_2=statstable.trial==2;
    data_t2=statstable.(idata)(sel_2,:);
    
    p2=ranksum(data_t1,data_t2);
    clear sel sel_2

    pvaltable.(idata)=[p p2]';
end %idata

Row_Names={'groups','trials'}';
fname_pval= fullfile(anovastatpath,'pval_measured_data');
pvaltable.vs=Row_Names;
writetable(pvaltable,fname_pval,'FileType','spreadsheet');

%% Make boxplots of comparisons

%for origin time and depth
for idata= ["peakval" "minslope" "halfwidth"]
    
    %separate by depth class
    boxplot(statstable.(idata),statstable.oriclass,'Colors','k','Symbol','o');
    
    fig=gcf;
    fname_boxplot_depth=fullfile(boxplotpath,sprintf('%s_per_depth',idata));
    dtx_savefigure(fig,fname_boxplot_depth,'pdf','png','close');
    
    %separate by trial
    boxplot(statstable.(idata),statstable.trial,'Colors','k','Symbol','o');
    
    fig=gcf;
    fname_boxplot_trial=fullfile(boxplotpath,sprintf('%s_per_trial',idata));
    dtx_savefigure(fig,fname_boxplot_trial,'pdf','png','close');
    
end %idata


end