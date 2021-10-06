%% set parameters
try %en local
    scriptpath = matlab.desktop.editor.getActiveFilename;
catch %cluster
    scriptpath = mfilename('fullpath');
end

epicodepath = [fileparts(fileparts(fileparts(fileparts(scriptpath)))), filesep];

addpath (genpath([epicodepath,'development']))
addpath (genpath([epicodepath,'shared']))
addpath (genpath([epicodepath,'external']))
addpath (genpath([epicodepath,'templates']))
addpath (genpath([epicodepath,'projects', filesep, 'wod']))
addpath (genpath([epicodepath,'projects', filesep, 'dtx']))
addpath (genpath([epicodepath,'projects', filesep, 'wod',filesep,'wod_functions']))
addpath (genpath([epicodepath,'projects', filesep, 'wod',filesep,'Intra']))
addpath (genpath([epicodepath,'projects', filesep, 'wod',filesep,'DC']))


if ispc
    addpath \\lexport\iss01.charpier\analyses\wod\fieldtrip-20200607
elseif isunix
    addpath /network/lustre/iss01/charpier/analyses/wod/fieldtrip-20200607
end

ft_defaults

config=DC_setparams;

tablesavedir=fullfile(fileparts(fileparts(config{4}.datasavedir)),'tables','matlab');
boxplotpath=fullfile(config{1}.imagesavedir,'boxplots');
statssavedir=fullfile(fileparts(fileparts(config{4}.datasavedir)),'stats','DC');
averagetracepath=fullfile(config{1}.imagesavedir,'average_trace');

if ~isfolder(boxplotpath)
    mkdir(boxplotpath);
end


if ~isfolder(averagetracepath)
    mkdir(averagetracepath);
end

if ~isfolder(statssavedir)
    mkdir(statssavedir);
end
%% GET DATA
chanlist = ["Vm", "DC"];

for iprot=1:size(config,2)
    matlabdata_savedir=fullfile(config{iprot}.datasavedir,'Matlab_struct');
    if ~isfile(fullfile(matlabdata_savedir,sprintf('%s_Vm.mat',config{iprot}.directorylist{1}{1})))
        continue
    end
    
    %get raw superficial DC
    temp= load(fullfile(matlabdata_savedir,sprintf('%s_Vm.mat',config{iprot}.directorylist{1}{1})));
    Superficial_raw{iprot}=temp.data_raw;
    clear temp
    
    if ~isfile(fullfile(matlabdata_savedir,sprintf('%s_DC.mat',config{iprot}.directorylist{1}{1})))
        continue
    end
    
    %get raw deep DC
    temp= load(fullfile(matlabdata_savedir,sprintf('%s_DC.mat',config{iprot}.directorylist{1}{1})));
    Deep_raw{iprot}=temp.data_raw;
    clear temp
    
    if ~isfile(fullfile(matlabdata_savedir,sprintf('%s_Vm_filt.mat',config{iprot}.directorylist{1}{1})))
        continue
    end
    
    %get filtered superficial DC
    temp= load(fullfile(matlabdata_savedir,sprintf('%s_Vm_filt.mat',config{iprot}.directorylist{1}{1})));
    Superficial_filt{iprot}=temp.data_filt;
    clear temp
    
    
    if ~isfile(fullfile(matlabdata_savedir,sprintf('%s_DC_filt.mat',config{iprot}.directorylist{1}{1})))
        continue
    end
    
    
    %get filtered deep DC
    temp= load(fullfile(matlabdata_savedir,sprintf('%s_DC_filt.mat',config{iprot}.directorylist{1}{1})));
    Deep_filt{iprot}=temp.data_filt;
    clear temp
    
    
    %get events
    temp= load(fullfile(matlabdata_savedir,sprintf('%s_Vm_events.mat',config{iprot}.directorylist{1}{1})));
    Events{iprot}=temp.CEDStruct;
    clear temp
end %ineuron
%% Get AD amplitude, time and min and max slope
for iprot= 1:size(config,2)
    
    %select AD
    t_1=Events{iprot}.markers.VentOff.synctime+20;
    t_2=Events{iprot}.markers.WoD.synctime+10;
    t_sel=[t_1 t_2];
    
    cfgtemp=[];
    cfgtemp.latency=t_sel;
    AD_sup{iprot}=ft_selectdata(cfgtemp,Superficial_raw{iprot});
    AD_deep{iprot}=ft_selectdata(cfgtemp,Deep_raw{iprot});
    
    clear t_1 t_2 t_sel
    %find peaks
    [v_peak_sup, t_peak_sup]=findpeaks(-AD_sup{iprot}.trial{1},AD_sup{iprot}.time{1},'NPeaks',1,'SortStr','descend');
    [v_peak_deep, t_peak_deep]=findpeaks(-AD_deep{iprot}.trial{1},AD_deep{iprot}.time{1},'NPeaks',1,'SortStr','descend');
    
    amp_peak_sup= -v_peak_sup;
    time_peak_sup=t_peak_sup;
    amp_peak_deep= -v_peak_deep;
    time_peak_deep=t_peak_deep;
    
    %store values
    statsDC{iprot}.Sup.amp_peak=amp_peak_sup;
    statsDC{iprot}.Sup.time_peak=time_peak_sup-Events{iprot}.markers.VentOff.synctime;
    statsDC{iprot}.Deep.amp_peak=amp_peak_deep;
    statsDC{iprot}.Deep.time_peak=time_peak_deep-Events{iprot}.markers.VentOff.synctime;
    
    
    clear amp_peak_sup amp_peak_deep time_peak_sup time_peak_deep v_peak_deep v_peak_sup t_peak_sup t_peak_deep
    
    %make derivative of signal
    AD_supSlope{iprot}=AD_sup{iprot};
    AD_deepSlope{iprot}=AD_deep{iprot};
    AD_supSlope{iprot}.trial{1}=ft_preproc_derivative(AD_sup{iprot}.trial{1},1);
    AD_deepSlope{iprot}.trial{1}=ft_preproc_derivative(AD_deep{iprot}.trial{1},1);
    
    %get peak of negative slope
    [~,max_slope_sup]=findpeaks(-AD_supSlope{iprot}.trial{1},AD_supSlope{iprot}.time{1},'NPeaks',1,'SortStr','descend');
    [~,max_slope_deep]=findpeaks(-AD_deepSlope{iprot}.trial{1},AD_deepSlope{iprot}.time{1},'NPeaks',1,'SortStr','descend');
    
    %store value
    statsDC{iprot}.Sup.min_slope=max_slope_sup;
    statsDC{iprot}.Deep.min_slope=max_slope_deep;
    
    clear max_slope_sup max_slope_deep
    
    %select rising phase
    t_1=Events{iprot}.markers.VentOn.synctime;
    t_2=Events{iprot}.markers.VentOn.synctime+60;
    t_sel=[t_1 t_2];
    
    cfgtemp=[];
    cfgtemp.latency=t_sel;
    Risin_sup{iprot}=ft_selectdata(cfgtemp,Superficial_raw{iprot});
    Risin_deep{iprot}=ft_selectdata(cfgtemp,Deep_raw{iprot});
    
    clear t_1 t_2 t_sel
    
    %derivate signal
    Risin_sup{iprot}.trial{1}=ft_preproc_derivative(Risin_sup{iprot}.trial{1},1);
    Risin_deep{iprot}.trial{1}=ft_preproc_derivative(Risin_deep{iprot}.trial{1},1);
    
    %get peak of positive slope
    [~,max_slope_sup]=findpeaks(Risin_sup{iprot}.trial{1},Risin_sup{iprot}.time{1},'NPeaks',1,'SortStr','descend');
    [~,max_slope_deep]=findpeaks(Risin_deep{iprot}.trial{1},Risin_deep{iprot}.time{1},'NPeaks',1,'SortStr','descend');
    
    %store value
    statsDC{iprot}.Sup.max_slope=max_slope_sup;
    statsDC{iprot}.Deep.max_slope=max_slope_deep;
    
    clear max_slope_sup max_slope_deep
end
%% Make table with data

%make table with all comparison variables
    datatable = table.empty;
    irow=0;
    for irat=1:size(config,2)
        if isempty(config{irat})
            continue
        end
        
        for idepth=["Sup" "Deep"]
            
            if idepth=="Sup"
                depth_idx=1;
            else
                depth_idx=2;
            end
            
            irow= irow+1;
            datatable.rat(irow)=irat;
            if depth_idx==1
            datatable.depth(irow)=config{irat}.DC.sup_dep;
            else
            datatable.depth(irow)=config{irat}.DC.dep_dep;
            end
            datatable.amp_peak(irow)=abs(statsDC{irat}.(idepth).amp_peak);
            datatable.time_peak(irow)=statsDC{irat}.(idepth).time_peak;
            datatable.min_slope(irow)=statsDC{irat}.(idepth).min_slope;
            datatable.max_slope(irow)=statsDC{irat}.(idepth).max_slope;
            datatable.depthclass(irow)=depth_idx;
        end %idepth
    end %irat
    
    fname_DC=fullfile(tablesavedir,'table_DC');
    writetable(datatable,fname_DC,'FileType','Spreadsheet')
%% Make non parametric comparisons of means
pValtable=table.empty;

sel_1=datatable.depthclass==1;
sel_2=datatable.depthclass==2;

irow=0;
for idata= ["amp_peak" "time_peak" "min_slope" "max_slope"]
    irow=irow+1;
    sup=datatable.(idata)(sel_1);
    deep=datatable.(idata)(sel_2);
    
    p(irow)=signrank(sup,deep,'method','exact');
    
    pValtable.(idata)=p(irow);
end %idata

%correct p values
[~,~,~,adj_p]=fdr_bh(p);


fname_Pval=fullfile(statssavedir,'table_pVal');
writetable(pValtable,fname_Pval,'FileType','Spreadsheet');

clear fname_Pval pValtable
%% Make boxplots of comparisons

for idata= ["amp_peak" "time_peak" "min_slope" "max_slope"]
    
    toplot.(idata)=[data_L5.(idata)',data_L23.(idata)'];
    
    
    %Make boxplots of comparisons
    boxplot(toplot.(idata),'Colors','k','Symbol','o','Labels',{'Layer 5','Layer 2/3'});
    
    fig=gcf;
    fname_boxplot_depth=fullfile(boxplotpath,sprintf('%s_per_layer',idata));
    dtx_savefigure(fig,fname_boxplot_depth,'pdf','png','close');
    
end %idata
%% Make average traces

for iprot=1:size(config,2)
    
    %select only AD
    t_1= Events{iprot}.markers.VentOff.synctime;
    t_2_sup= statsDC{iprot}.Sup.time_peak+t_1+10;
    t_2_deep= statsDC{iprot}.Deep.time_peak+t_1+10;
    t_sel_sup=[t_1 t_2_sup];
    t_sel_deep=[t_1 t_2_deep];
    
    
    cfgtemp=[];
    cfgtemp.latency=t_sel_sup;
    Sup_raw_AD{iprot}=ft_selectdata(cfgtemp,Superficial_raw{iprot});
    cfgtemp.latency= t_sel_deep;
    Deep_raw_AD{iprot}=ft_selectdata(cfgtemp,Deep_raw{iprot});
    clear t_1 t_2_sup t_sel_sup t_sel_deep t_2_deep
    
    %normalize time tVoff=0 tpic= 1
    t_old_sup                                                = Sup_raw_AD{iprot}.time{1};
    t_new_sup                                                = (t_old_sup-min(Sup_raw_AD{iprot}.time{1}))/(max(Sup_raw_AD{iprot}.time{1})-min(Sup_raw_AD{iprot}.time{1}));
    Sup_raw_AD_new{iprot}                                    = Sup_raw_AD{iprot};
    Sup_raw_AD_new{iprot}.time{1}                            = t_new_sup;
    
    t_old_deep                                                = Deep_raw_AD{iprot}.time{1};
    t_new_deep                                                = (t_old_deep-min(Deep_raw_AD{iprot}.time{1}))/(max(Deep_raw_AD{iprot}.time{1})-min(Deep_raw_AD{iprot}.time{1}));
    Deep_raw_AD_new{iprot}                                    = Deep_raw_AD{iprot};
    Deep_raw_AD_new{iprot}.time{1}                            = t_new_deep;
    
    
    cfgtemp = [];
    cfgtemp.time = {0:1/180:1};
    Sup_raw_AD_new{iprot} = ft_resampledata(cfgtemp, Sup_raw_AD_new{iprot});
    Deep_raw_AD_new{iprot} = ft_resampledata(cfgtemp, Deep_raw_AD_new{iprot});
    
    tomean_sup(iprot,:)=Sup_raw_AD_new{iprot}.trial{1};
    tomean_deep(iprot,:)=Deep_raw_AD_new{iprot}.trial{1};
    time=Sup_raw_AD_new{iprot}.time{1};
    
end %iprot

%Make substranction of traces sup - deep
for icol=1:size(tomean_sup,2)
    for irow=1:size(tomean_sup,1)
    diff_trace(irow,icol)=tomean_sup(irow,icol)-tomean_deep(irow,icol);
    end
end


%Calculate mean traces
mean_diff=mean(diff_trace,1);
mean_sup=mean(tomean_sup,1);
mean_deep=mean( tomean_deep,1);
std_diff=std(diff_trace,0,1);
std_sup=std(tomean_sup,0,1);
std_deep=std(tomean_deep,0,1);
sem_sup=std_sup/sqrt(size(tomean_sup,1));
sem_deep=std_deep/sqrt(size(tomean_deep,1));
sem_diff=std_diff/sqrt(size(diff_trace,1));
%% Make comparison of traces


zero_trace=zeros(size(diff_trace));

pVal_table=table.empty;
for icol=1:size(tomean_sup,2)
   p(icol)=signrank(tomean_sup(:,icol),tomean_deep(:,icol),'method','exact');
   p_diff(icol)=ranksum(diff_trace(:,icol),zero_trace(:,icol));
end
%correct p-values
[~,~,~,adj_p]=fdr_bh(p);
[~,~,~,adj_p_diff]=fdr_bh(p_diff);

NoSigni=adj_p>0.01;
NoSigniDiff=adj_p_diff>0.01;
pVal_plot=ones(1,size(adj_p,2))*10;
pVal_plot(NoSigni)=nan;

pVal_plotDiff=ones(1,size(adj_p_diff,2))*20;
pVal_plotDiff(NoSigniDiff)=nan;

pVal_table.pval=p';
pVal_table.adj_Pval=adj_p';
pVal_table.diff_pval=p_diff';
pVal_table.diff_adj_pval=adj_p_diff';

fname_pVal=fullfile(statssavedir,'average_traces');
writetable(pVal_table,fname_pVal,'FileType','Spreadsheet');
%% Plot Average traces

fig_L23=figure;hold;
plot(time,mean_sup,'Color',[0.8500 0.3250 0.0980],'LineWidth',1);
plot(time,pVal_plot,'r','LineWidth',2);
patch( [time(1,:),time(1,end:-1:1)],[mean_sup- sem_sup, mean_sup(end:-1:1)+ sem_sup(end:-1:1)], [0.8500 0.3250 0.0980], 'facealpha', 0.3, 'edgecolor', 'none');

ylim([-30 10]);

fname_23=fullfile(averagetracepath,'DC_superficial_mean_sem');
dtx_savefigure(fig_L23,fname_23,'png','pdf','close');

fig_L5=figure;hold;
plot(time(1,:),mean_deep,'Color','k','LineWidth',1);
plot(time,pVal_plot,'r','LineWidth',2);
patch( [time(1,:),time(1,end:-1:1)],[mean_deep- sem_deep, mean_deep(end:-1:1)+ sem_deep(end:-1:1)], 'k', 'facealpha', 0.3, 'edgecolor', 'none');

ylim([-30 10]);
fname_5=fullfile(averagetracepath,'DC_deep_mean_sem');
dtx_savefigure(fig_L5,fname_5,'png','pdf','close');

fig_diff=figure;hold;
plot(time(1,:),mean_diff,'Color','r','LineWidth',1);
patch( [time(1,:),time(1,end:-1:1)],[mean_diff- sem_diff, mean_diff(end:-1:1)+ sem_diff(end:-1:1)], 'r', 'facealpha', 0.3, 'edgecolor', 'none');

plot(time(1,:),mean_deep,'Color','k','LineWidth',1);
patch( [time(1,:),time(1,end:-1:1)],[mean_deep- sem_deep, mean_deep(end:-1:1)+ sem_deep(end:-1:1)], 'k', 'facealpha', 0.3, 'edgecolor', 'none');

plot(time,mean_sup,'Color',[0.8500 0.3250 0.0980],'LineWidth',1);
patch( [time(1,:),time(1,end:-1:1)],[mean_sup- sem_sup, mean_sup(end:-1:1)+ sem_sup(end:-1:1)], [0.8500 0.3250 0.0980], 'facealpha', 0.3, 'edgecolor', 'none');

fname_diff=fullfile(averagetracepath,'DC_bothtraces_subtracted');
dtx_savefigure(fig_diff,fname_diff,'png','pdf','close');
%% Plot single trials and overdraw filtered signal
singletrialpath=fullfile(config{1}.imagesavedir,'filter_DC_overdraw');

if ~isfolder(singletrialpath)
        mkdir(singletrialpath);
    end

for iprot=1:size(config,2)
   
   fig_sup=figure;hold 
   plot(Superficial_raw{iprot}.time{1},Superficial_raw{iprot}.trial{1},'Color',[0.8500 0.3250 0.0980],'LineWidth',1);
   plot(Superficial_filt{iprot}.time{1},Superficial_filt{iprot}.trial{1},'Color','r','LineWidth',1,'LineStyle','--');
   
   fname_sup=fullfile(singletrialpath,sprintf('protocol_%i_superficial',iprot));
   dtx_savefigure(fig_sup,fname_sup,'pdf','png','close');
   
   fig_deep=figure;hold 
   plot(Deep_raw{iprot}.time{1},Deep_raw{iprot}.trial{1},'Color','k','LineWidth',1);
   plot(Deep_filt{iprot}.time{1},Deep_filt{iprot}.trial{1},'Color','r','LineWidth',1,'LineStyle','--');
   
   fname_deep=fullfile(singletrialpath,sprintf('protocol_%i_deep',iprot));
   dtx_savefigure(fig_deep,fname_deep,'pdf','png','close');

end
