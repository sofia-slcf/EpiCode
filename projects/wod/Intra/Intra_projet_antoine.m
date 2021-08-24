function Intra_projet_antoine(slurm_task_id)


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
    addpath \\lexport\iss01.charpier\echanges\scripts-paul\Spike2_vers_MATLAB
elseif isunix
    addpath /network/lustre/iss01/charpier/analyses/wod/fieldtrip-20200607
    addpath /network/lustre/iss01/charpier/echanges/scripts-paul/Spike2_vers_MATLAB
    
end
ft_defaults

%remove fieldtrip's output
ft_warning off
ft_notice off
ft_info off
ft_debug off
global ft_default
ft_default.checkconfig = 'silent';
ft_default.checkpath = 'once';
ft_default.showcallinfo = 'no';
ft_default.trackcallinfo = 'no';
ft_default.tracktimeinfo = 'no';

config=Intra_setparams;
configDC=DC_setparams;

calculated_datapath=fullfile(config{1}.datasavedir,'calculated_data');

if ~isfolder(calculated_datapath)
    mkdir(calculated_datapath)
end


if slurm_task_id==0
    
    %% GET DATA
    
    for ineuron=1:size(config,2)
        matlabdata_savedir=fullfile(config{ineuron}.rawdir,'matlab_structures');
        %get Vm
        temp= load(fullfile(matlabdata_savedir,sprintf('%s_Vm.mat',config{ineuron}.prefix)));
        Vm{ineuron}=temp.data;
        clear temp
        
        %get events
        temp= load(fullfile(matlabdata_savedir,sprintf('%s_events.mat',config{ineuron}.prefix)));
        Events{ineuron}=temp.CEDStruct;
        clear temp
        
        %get EEG Intra
        temp= load(fullfile(matlabdata_savedir,sprintf('%s_EEG-S1-L',config{ineuron}.prefix)));
        EEG{ineuron}=temp.data;
        clear temp
        
        
    end %ineuron
    %% Baseline correct and remove APs from Vm
    
    for ineuron=1:size(config,2)
        Vm_APless{ineuron}=Vm{ineuron};
        
        % find AP
        [~, t_ap]=findpeaks(Vm{ineuron}.trial{1},Vm{ineuron}.time{1},'MinPeakProminence',20);
        
        % replace neighbouring samples of AP timings by NaNs Fsample= 1000Hz
        for i= 1:size(t_ap,2)
            idx_ap=find( Vm_APless{ineuron}.time{1}==t_ap(i));
            
            remove_ap= [idx_ap-3 :1: idx_ap+3];
            
            if idx_ap-3 < 0 || idx_ap+3> length(Vm_APless{ineuron}.time{1})
                remove_ap= [idx_ap-1 :1: idx_ap+1];
            end
            
            if size(remove_ap,2)==3
                continue
            end
            
            %replace +/- 3 samples aound AP by NaNs
            Vm_APless{ineuron}.trial{1}(1,remove_ap)=NaN;
        end %ap
        
        %interpolate
        Vm_APless{ineuron}.trial{1}=fillmissing(Vm_APless{ineuron}.trial{1},'pchip');
        
        %Determine baseline
        t1= Events{ineuron}.markers.VentOff.synctime-60;
        t2= Events{ineuron}.markers.VentOff.synctime;
        t_sel=[t1 t2];
        
        cfgtemp=[];
        cfgtemp.latency=t_sel;
        baseline{ineuron}= ft_selectdata(cfgtemp,Vm_APless{ineuron});
        clear t1 t2 t_sel
        
        
        %average baseline Vm
        baseline_mean=nanmean(baseline{ineuron}.trial{1});
        %correct Vm with mean baseline Vm
        Vm_APless_norm{ineuron}=Vm_APless{ineuron};
        Vm_APless_norm{ineuron}.trial{1}=Vm_APless{ineuron}.trial{1}-baseline_mean;
        
        %store mean Vm in data
        stats_intra{ineuron}.Vm_rest=baseline_mean;
    end %ineuron
    %% Overdraw V_off aligned smoothed Vm
    
    color={[0.5 0.5 0.5],[0.8500 0.3250 0.0980]};
    fig_ad=figure;hold
    for ineuron=1:size(config,2)
        
        t_1=Events{ineuron}.markers.VentOff.synctime-60;
        t_2=Events{ineuron}.markers.WoD.synctime+20;
        t_sel=[t_1 t_2];
        
        %select plot window
        cfgtemp=[];
        cfgtemp.latency=t_sel;
        data_plot= ft_selectdata(cfgtemp,Vm_APless_norm{ineuron});
        
        %filter plot data
        cfgtemp=[];
        cfgtemp.lpfilter='yes';
        cfgtemp.lpfreq=5;
        cfgtemp.lpfilttype='fir';
        data_plot=ft_preprocessing(cfgtemp,data_plot);
        
        %smoothing
        data_plot.trial{1}=movmean(data_plot.trial{1},100);
        
        
        data_plot.time{1}=data_plot.time{1}-Events{ineuron}.markers.VentOff.synctime;
        
        if config{ineuron}.Intra.dep >= 1000
            c=color{1};
        else
            c=color{2};
        end
        
        plot(data_plot.time{1},data_plot.trial{1},'Color',c,'LineWidth',1.5);
        
    end %ineuron
    fname_ad= fullfile(config{1}.imagesavedir,'smooth_depol','voff_aligned');
    
    if ~isfolder(fullfile(config{1}.imagesavedir,'smooth_depol'))
        mkdir(fullfile(config{1}.imagesavedir,'smooth_depol'));
    end
    
    dtx_savefigure(fig_ad,fname_ad,'png','pdf','close');
    %% Calculate deltaVm AD and AD timing
    
    for ineuron=1:size(config,2)
        
        %take Vm AD
        t_1=Events{ineuron}.markers.WoD.synctime+3;
        t_2=Events{ineuron}.markers.WoD.synctime+5;
        t_sel=[t_1 t_2];
        
        cfgtemp=[];
        cfgtemp.latency=t_sel;
        AD_Intra=ft_selectdata(cfgtemp,Vm_APless{ineuron});
        
        clear t_1 t_2 t_sel
        
        Vm_AD= mean(AD_Intra.trial{1});
        dVm= Vm_AD-stats_intra{ineuron}.Vm_rest;
        stats_intra{ineuron}.depol=dVm;
        
        %take time of max AD
        t_1=Events{ineuron}.markers.VentOff.synctime+30;
        t_2=Events{ineuron}.markers.WoD.synctime+5;
        t_sel=[t_1 t_2];
        
        cfgtemp=[];
        cfgtemp.latency=t_sel;
        AD_search= ft_selectdata(cfgtemp,Vm_APless{ineuron});
        
        idx_AD=AD_search.trial{1}> Vm_AD;
        idx_AD=find(idx_AD,1,'first');
        
        time_AD=AD_search.time{1}(idx_AD);
        
        fig=figure;
        plot(AD_search.time{1},AD_search.trial{1})
        xline(time_AD)
        
        detectionpath=fullfile(config{ineuron}.imagesavedir,'detection');
        
        if ~isfolder(detectionpath)
            mkdir(detectionpath);
        end
        fnamefig=fullfile(detectionpath,sprintf('AD_Vm_measure_neuron_%i',ineuron));
        dtx_savefigure(fig,fnamefig,'png','pdf','close');
        
        stats_intra{ineuron}.AD_time=time_AD-Events{ineuron}.markers.VentOff.synctime;
        
        
        %% Calculate max slope of AD
        
        %smooth AD signal
        AD_search_smooth=AD_search;
        AD_search_smooth.trial{1}=movmean(AD_search.trial{1}(1,:),1000);
        
        AD_slope=ft_preproc_derivative(AD_search_smooth.trial{1},1);
        plot(AD_search.time{1},AD_search_smooth.trial{1})
        plot(AD_search.time{1},AD_slope)
        
        t=AD_search_smooth.time{1};
        t_1=t< Events{ineuron}.markers.WoD.synctime+10;
        t_2= t> Events{ineuron}.markers.WoD.synctime-10;
        t_sel= t_1 & t_2;
        
        [v_peak_slope, t_peak_slope] = findpeaks( AD_slope(t_sel),t(t_sel),'NPeaks',1,'SortStr','descend','WidthReference','Halfheight');
        clear t t_1 t_2 t_sel
        
        
        stats_intra{ineuron}.AD_slope=v_peak_slope;
    end
    %% Determine Time of ISO
    for ineuron=1:size(config,2)
        [up,lo]=envelope(baseline{ineuron}.trial{1},10000,'peak');
        baseline_amp=mean(up-lo);
        clear up lo
        
        [up,lo]=envelope(Vm_APless{ineuron}.trial{1},10000,'peak');
        amp_signal=up-lo;
        clear up lo
        thr=0.1*baseline_amp;
        
        
        %         plot(Vm_APless{ineuron}.time{1},amp_signal);
        %         hold
        %         yline(thr)
        
        idx_signal= amp_signal>thr(1);
        diff_signal=diff(idx_signal);
        thr_cross=find(diff_signal);
        idx_ISO=diff(thr_cross)>10000;
        ISO_time= Vm_APless{ineuron}.time{1}(thr_cross(idx_ISO));
        idx=  ISO_time> 0.2*max(Vm_APless{ineuron}.time{1});
        ISO_time=ISO_time(find(idx,1,'first'));
        
        clear idx_signal diff_signal thr_cross idx_ISO idx thr_cross
        
        stats_intra{ineuron}.iso_delay=ISO_time-Events{ineuron}.markers.VentOff.synctime;
        stats_intra{ineuron}.iso_time=ISO_time;
        
        
        
        fig_detect=figure;
        plot(Vm_APless{ineuron}.time{1},Vm_APless{ineuron}.trial{1})
        xline(ISO_time);
        xlim([Events{ineuron}.markers.VentOff.synctime+10,Events{ineuron}.markers.WoD.synctime])
        
        fname_detect=fullfile(detectionpath,sprintf('Time_ISO_neuron_%i',ineuron));
        dtx_savefigure(fig_detect,fname_detect,'pdf','png','close');
    end %ineuron
    
    fname_dataall=fullfile(calculated_datapath,'data_allneurons.mat');
    save(fname_dataall,'stats_intra');
    %% Make average trace of normalized signal
    
    for ineuron=1:size(config,2)
        
        %select only AD
        t_1= Events{ineuron}.markers.VentOff.synctime;
        t_2= Events{ineuron}.markers.pic.synctime;
        t_sel=[t_1 t_2];
        
        cfgtemp=[];
        cfgtemp.latency=t_sel;
        AD_tomean{ineuron}=ft_selectdata(cfgtemp,Vm_APless_norm{ineuron});
        clear t_1 t_2 t_sel
        
        %normalize time tVoff=0 tpic= 1
        t_old                                                = AD_tomean{ineuron}.time{1};
        t_new                                                = (t_old-min(AD_tomean{ineuron}.time{1}))/(max(AD_tomean{ineuron}.time{1})-min(AD_tomean{ineuron}.time{1}));
        AD_tomean_new{ineuron}                               = AD_tomean{ineuron};
        AD_tomean_new{ineuron}.time{1}                       = t_new;
        
        cfgtemp = [];
        cfgtemp.time = {0:1/180:1};
        AD_tomean_new{ineuron} = ft_resampledata(cfgtemp, AD_tomean_new{ineuron});
        
        length_data(ineuron)=size(AD_tomean_new{ineuron}.time{1},2);
        
        if config{ineuron}.Intra.rename{1}=='Intra_dep'
            tomean_L5.Vm(ineuron,:)=AD_tomean_new{ineuron}.trial{1};
            tomean_L5.time(ineuron,:)=AD_tomean_new{ineuron}.time{1};
        elseif config{ineuron}.Intra.rename{1}=='Intra_sup'
            tomean_L23.Vm(ineuron,:)=AD_tomean_new{ineuron}.trial{1};
            tomean_L23.time(ineuron,:)=AD_tomean_new{ineuron}.time{1};
        end
        
    end %ineuron
    
    %Delete empty rows
    for i=1:size(tomean_L23.Vm,1)
        sumRowVm=sum(tomean_L23.Vm(i,:));
        if sumRowVm==0
            tomean_L23.Vm(i,:)=[];
            tomean_L23.time(i,:)=[];
        end
    end
    
    for i=1:size(tomean_L5.Vm,1)
        sumRowVm=sum(tomean_L5.Vm(i,:));
        if sumRowVm==0
            tomean_L5.Vm(i,:)=[];
            tomean_L5.time(i,:)=[];
        end
    end
    
    %Calculate mean traces
    mean_L23=mean(tomean_L23.Vm,1);
    mean_L5=mean(tomean_L5.Vm,1);
    std_L23=std(tomean_L23.Vm,0,1);
    std_L5=std(tomean_L5.Vm,0,1);
    sem_L23=std_L23/sqrt(size(tomean_L23.Vm,1));
    sem_L5=std_L5/sqrt(size(tomean_L5.Vm,1));
    %% Make non parametric comparisons of average traces
    
    statssavedir='\\lexport\iss01.charpier\analyses\wod\Antoine\stats\Intra';
    
 
    %make wilcoxon test for each sample
    for icol=1:size(tomean_L23.Vm,2)
        p(icol)=ranksum(tomean_L23.Vm(:,icol),tomean_L5.Vm(:,icol));
    end
    
    %correct p-values
    [~,~,~,adj_p]=fdr_bh(p);
    
    %make table with p-values
    pValuestable= table.empty;
    pValuestable.original_p=p';
    pValuestable.adjusted_p=adj_p';
    
    fname_table=fullfile(statssavedir,'Average_traces_pval');
    writetable(pValuestable,fname_table,'FileType','Spreadsheet');
    
    %prepare line to plot on average traces
    idx_signif= p<0.05;
    idx_nnsignif= p>0.05;
    
    pVal_plot=ones(1,size(p,2));
    pVal_plot=pVal_plot.*idx_signif;
    pVal_plot=pVal_plot*60;
    pVal_plot(idx_nnsignif)=nan;
    %% Plot Overdraw with average trace and mean +/- SEM
    fig_L23=figure;hold
    for ineuron= 1:size(tomean_L23.Vm,1)
        plot(tomean_L23.time(ineuron,:), tomean_L23.Vm(ineuron,:),'Color',[0.8500 0.3250 0.0980 0.4])
    end
    plot(tomean_L23.time(1,:),mean_L23,'Color','r','LineWidth',1);
    plot(tomean_L5.time(1,:),pVal_plot,'Color','r','LineWidth',2);
    
    averagetracepath=fullfile(config{1}.imagesavedir,'average_trace');
    
    if ~isfolder(averagetracepath)
        mkdir(averagetracepath);
    end
    
    fname_23=fullfile(averagetracepath,'Intra_superficial');
    dtx_savefigure(fig_L23,fname_23,'png','pdf','close');
    
    
    fig_L5=figure;hold
    for ineuron= 1:size(tomean_L23.Vm,1)
        plot(tomean_L5.time(ineuron,:), tomean_L5.Vm(ineuron,:),'Color',[0 0 0 0.4])
    end
    plot(tomean_L5.time(1,:),mean_L5,'Color','r','LineWidth',1);
    plot(tomean_L5.time(1,:),pVal_plot,'Color','r','LineWidth',2);
    
    
    fname_5=fullfile(averagetracepath,'Intra_deep');
    dtx_savefigure(fig_L5,fname_5,'png','pdf','close');
    
    fig_L23=figure;hold;
    plot(tomean_L23.time(1,:),mean_L23,'Color',[0.8500 0.3250 0.0980],'LineWidth',1);
    patch( [tomean_L23.time(1,:),tomean_L23.time(1,end:-1:1)],[mean_L23- sem_L23, mean_L23(end:-1:1)+ sem_L23(end:-1:1)], [0.8500 0.3250 0.0980], 'facealpha', 0.3, 'edgecolor', 'none');
    plot(tomean_L5.time(1,:),pVal_plot,'Color','r','LineWidth',2);
    
    ylim([-10 60]);
    
    fname_23=fullfile(averagetracepath,'Intra_superficial_mean_sem');
    dtx_savefigure(fig_L23,fname_23,'png','pdf','close');
    
    fig_L5=figure;hold;
    plot(tomean_L5.time(1,:),mean_L5,'Color','k','LineWidth',1);
    patch( [tomean_L5.time(1,:),tomean_L5.time(1,end:-1:1)],[mean_L5- sem_L5, mean_L5(end:-1:1)+ sem_L5(end:-1:1)], 'k', 'facealpha', 0.3, 'edgecolor', 'none');
    plot(tomean_L5.time(1,:),pVal_plot,'Color','r','LineWidth',2);
    
    ylim([-10 60]);
    
    fname_5=fullfile(averagetracepath,'Intra_deep_mean_sem');
    dtx_savefigure(fig_L5,fname_5,'png','pdf','close');
    %% GET DATA measured data
    temp=load(fullfile(calculated_datapath,'data_allneurons.mat'));
    stats_intra=temp.stats_intra;
    clear temp
    %% Make table with all neuron data
    table_csv= readtable(fullfile(config{1}.datasavedir,'table_neurons.csv'));
    
    neuron_table=table_csv;
    irow=0;
    for ineuron=1:size(config,2)
        
        if config{ineuron}.Intra.rename{1}=='Intra_dep'
            intraclass=1;
        elseif config{ineuron}.Intra.rename{1}=='Intra_sup'
            intraclass=2;
        end
        irow=irow+1;
        neuron_table.depthclass(irow)=intraclass;
        neuron_table.AD(irow)= stats_intra{ineuron}.depol;
        neuron_table.depth(irow)= config{ineuron}.Intra.dep;
        neuron_table.AD_time(irow)=stats_intra{ineuron}.AD_time;
        neuron_table.AD_slope(irow)= stats_intra{ineuron}.AD_slope;
        neuron_table.ISO_time(irow)=stats_intra{ineuron}.iso_delay;
    end %ineuron
    
    
    tablename=fullfile(config{1}.datasavedir,'AllNeuronsTable');
    
    writetable(neuron_table,tablename,'FileType','Spreadsheet');
    
    
    if ~isfolder(statssavedir)
        mkdir(statssavedir);
    end
    
    %% Non-parametric comparisons of means
    
    statstable=table.empty;
    for idata= ["Vm" "SDVm" "FR" "CV2" "Path" "Rm" "Tm" "AD" "AD_time" "AD_slope" "ISO_time"]
        
        %between layers
        sel=neuron_table.depthclass==1;
        data_deep=neuron_table.(idata)(sel,:);
        
        sel_2=neuron_table.depthclass==2;
        data_sup=neuron_table.(idata)(sel_2,:);
        
        p=ranksum(data_deep,data_sup);
        clear sel sel_2
        
        %store in table
        statstable.(idata)=p;
        
    end %idata
    fnamestattable= fullfile(statssavedir,'Intra_parameters_ranksum');
    writetable(statstable,fnamestattable,'Filetype','spreadsheet');
    %% Make boxplots of comparisons
    boxplotpath=fullfile(config{1}.imagesavedir,'boxplots');
    
    if ~isfolder(boxplotpath)
        mkdir(boxplotpath);
    end
    
    %for origin time and depth
    for  idata= ["Vm" "SDVm" "FR" "CV2" "Path" "Rm" "Tm" "AD" "AD_time" "AD_slope" "ISO_time"]
        
        %separate by depth class
        boxplot(neuron_table.(idata),neuron_table.depthclass,'Colors','k','Symbol','o');
        
        fig=gcf;
        fname_boxplot_depth=fullfile(boxplotpath,sprintf('%s_per_layer',idata));
        dtx_savefigure(fig,fname_boxplot_depth,'pdf','png','close');
        
    end %idata
    %% Make timefrequency analysis of APless Vm and plot
    
    for ineuron= [13,14,16]
        
        t_voff= Events{ineuron}.markers.VentOff.synctime;
        t_wod= Events{ineuron}.markers.WoD.synctime;
        t_start= t_voff-30;
        t_stop= stats_intra{ineuron}.iso_time;
        %downsample Vm_APless
        cfgtemp=[];
        cfgtemp.resamplefs=Vm_APless{ineuron}.fsample/4;
        Vm_resampled{ineuron}=ft_resampledata(cfgtemp,Vm_APless{ineuron});
        
        %select Vent-Off-30 to WoD period
        cfgtemp=[];
        cfgtemp.latency=[t_start t_stop];
        WOD{ineuron}= ft_selectdata(cfgtemp,Vm_resampled{ineuron});
        
        %select baseline
        cfgtemp=[];
        cfgtemp.latency=[t_voff-120 t_voff];
        Baseline{ineuron}= ft_selectdata(cfgtemp,Vm_resampled{ineuron});
        
        %highpass filter Vm
        cfgtemp=[];
        cfgtemp.hpfilter='yes';
        cfgtemp.hpfreq=1;
        cfgtemp.hpfilttype= 'fir';
        WOD_filt{ineuron}=ft_preprocessing(cfgtemp,WOD{ineuron});
        
        cfgtemp=[];
        cfgtemp.hpfilter='yes';
        cfgtemp.hpfreq=1;
        cfgtemp.hpfilttype= 'fir';
        Baseline_filt{ineuron}=ft_preprocessing(cfgtemp,Baseline{ineuron});
        
        %Make time frequency analysis of WOD period
        cfgtemp=[];
        cfgtemp.channel                 = 'all';
        cfgtemp.method='mtmconvol';
        cfgtemp.output                  = 'pow';
        cfgtemp.taper                   = 'dpss'; %default = dpss
        cfgtemp.pad                     = 'nextpow2';
        cfgtemp.foi=config{ineuron}.timefreq.foi;
        cfgtemp.tapsmofrq=config{ineuron}.timefreq.tapsmofrq;
        cfgtemp.toi='all';
        cfgtemp.t_ftimwin= ones(size(cfgtemp.foi))*config{ineuron}.timefreq.t_ftimwin;
        timefreq_wod= ft_freqanalysis(cfgtemp,WOD_filt{ineuron});
        
        %Make time frequency analysis of baseline period
        cfgtemp=[];
        cfgtemp.channel                 = 'all';
        cfgtemp.method='mtmconvol';
        cfgtemp.output                  = 'pow';
        cfgtemp.taper                   = 'dpss'; %default = dpss
        cfgtemp.pad                     = 'nextpow2';
        cfgtemp.foi=config{ineuron}.timefreq.foi;
        cfgtemp.tapsmofrq=config{ineuron}.timefreq.tapsmofrq;
        cfgtemp.toi='all';
        cfgtemp.t_ftimwin= ones(size(cfgtemp.foi))*config{ineuron}.timefreq.t_ftimwin;
        timefreq_baseline= ft_freqanalysis(cfgtemp,Baseline_filt{ineuron});
        
        %Make baseline correction
        timefreq_wod_blcorrected=timefreq_wod;
        for ifreq= 1:size(timefreq_wod.freq,2)
            baseline_ifreq = nanmean(squeeze(timefreq_baseline.powspctrm(1,ifreq,:))); %baseline
            timefreq_wod_blcorrected.powspctrm(1,ifreq,:)=timefreq_wod.powspctrm(1,ifreq,:) ./ baseline_ifreq;
        end
        
        %make Vent Off t0
        timefreq_wod.time=timefreq_wod.time-t_voff;
        timefreq_wod_blcorrected.time=timefreq_wod_blcorrected.time-t_voff;
        WOD{ineuron}.time{1}=WOD{ineuron}.time{1}-t_voff;
        EEG_plot{ineuron}=EEG{ineuron};
        EEG_plot{ineuron}.time{1}= EEG_plot{ineuron}.time{1}-t_voff;
        
        %filter DC LFP
        cfgtemp=[];
        cfgtemp.hpfilter='yes';
        cfgtemp.hpfreq=0.1;
        cfgtemp.hpfilttype='fir';
        EEG_plot{ineuron}=ft_preprocessing(cfgtemp,EEG_plot{ineuron});
        
        
        
        
        timefreqpath=fullfile(config{ineuron}.datasavedir,'timefreq');
        
        if ~isfolder(timefreqpath)
            mkdir(timefreqpath);
        end
        
        fname_timefreq_wod=fullfile(timefreqpath,sprintf('%s_%s_timefreq_wod',config{ineuron}.Intra.rename{1},config{ineuron}.prefix));
        fname_timefreq_wod_blcorrected=fullfile(timefreqpath,sprintf('%s_%s_timefreq_wod_blcorrected',config{ineuron}.Intra.rename{1},config{ineuron}.prefix));
        fname_timefreq_baseline=fullfile(timefreqpath,sprintf('%s_%s_timefreq_baseline',config{ineuron}.Intra.rename{1},config{ineuron}.prefix));
        
        %save timefreq data
        save(fname_timefreq_wod,'timefreq_wod');
        save(fname_timefreq_wod_blcorrected,'timefreq_wod_blcorrected');
        save(fname_timefreq_baseline,'timefreq_baseline');
        
    end %ineuron
    %% Plotting Time frequency
    for ineuron=[13,14,16]
        
        temp= load(fullfile(timefreqpath,sprintf('%s_%s_timefreq_wod_blcorrected',config{ineuron}.Intra.rename{1},config{ineuron}.prefix)));
        timefreq_wod_blcorrected=temp.timefreq_wod_blcorrected;
        clear temp
        
        t_voff= Events{ineuron}.markers.VentOff.synctime;
        t_wod= Events{ineuron}.markers.WoD.synctime;
        t_start= t_voff-30;
        t_stop= stats_intra{ineuron}.iso_time;
        
        fig=figure;
        subplot(3,1,3)
        %plot TFR
        %voir les paramètres optionnels dans le descriptifs de la fonction pour
        %modifier l'aspect du TFR. Avec les paramètres par défaut :
        cfgtemp         = [];
        cfgtemp.channel = 'all';
        cfgtemp.interactive = 'no';
        cfgtemp.colormap= 'jet';
        cfgtemp.fontsize = 12;
        cfgtemp.ylim= [1 100];
        cfgtemp.masknans    = 'yes';
        ft_singleplotTFR(cfgtemp, timefreq_wod_blcorrected);
        ft_pimpplot(fig, jet(5000))
        caxis([1 2])
        axis tight
        
        xlim([-20 t_stop-t_voff])
        subplot(3,1,2)
        plot(WOD{ineuron}.time{1},WOD{ineuron}.trial{1},'Color','k','LineWidth',0.5);
        axis tight
        
        xlim([-20 t_stop-t_voff])
        
        subplot(3,1,1)
        plot(EEG_plot{ineuron}.time{1},EEG_plot{ineuron}.trial{1},'Color','k','LineWidth',0.5);
        axis tight
        xlim([-20 t_stop-t_voff])
        
        
        fig_imagepath=fullfile(config{ineuron}.imagesavedir,'timefreq');
        
        if ~isfolder(fig_imagepath)
            mkdir(fig_imagepath);
        end
        
        fname_fig=fullfile(fig_imagepath,sprintf('%s_%s',config{ineuron}.Intra.rename{1},config{ineuron}.prefix));
        dtx_savefigure(fig,fname_fig,'png','pdf','close');
        
    end %ineuron
    %% Make V-off aligned Vm and ECoG
    
    for ineuron= 1:size(Vm,2)
        %get Vent Off timing
        t_voff= Events{ineuron}.markers.VentOff.synctime;
        
        %express data with this timing
        Vm{ineuron}.time{1}=Vm{ineuron}.time{1}-t_voff;
        EEG{ineuron}.time{1}=EEG{ineuron}.time{1}-t_voff;
        t_von=  Events{ineuron}.markers.VentOn.synctime-t_voff;
        
        
        fig=figure;
        subplot(2,1,2)
        plot(Vm{ineuron}.time{1},Vm{ineuron}.trial{1},'k','LineWidth',0.5);
        xlim([-10 t_von+60])
        
        
        subplot(2,1,1)
        plot(EEG{ineuron}.time{1},EEG{ineuron}.trial{1},'k','LineWidth',0.5);
        xlim([-10 t_von+60])
        
        
        fname=fullfile(config{ineuron}.imagesavedir,'EEG_Intra',sprintf('%s_V_off_aligned',config{ineuron}.prefix));
        dtx_savefigure(fig,fname,'png','pdf','close');
        
    end %ineuron
    %% Make cross-correlogram of LFP and Vm and plot
    
    for ineuron= [13,14,16]
        
        %read LF events
        Events{ineuron}=readCEDevents(fullfile(config{ineuron}.rawdir,config{ineuron}.directorylist{1}{1}));
        
        %filter DC recording
        if ineuron==16
            cfgtemp=[];
            cfgtemp.hpfilter='yes';
            cfgtemp.hpfreq=0.1;
            cfgtemp.hpfilttype='fir';
            EEG{ineuron}=ft_preprocessing(cfgtemp,EEG{ineuron});
        end
        
        %define period to do xcorr on
        t_1= Events{ineuron}.markers.LF.synctime(1);
        t_2= Events{ineuron}.markers.LF.synctime(2);
        t_sel=[t_1 t_2];
        
        cfgtemp=[];
        cfgtemp.latency=t_sel;
        LF_Vm{ineuron}=ft_selectdata(cfgtemp,Vm_APless{ineuron});
        LF_EEG{ineuron}=ft_selectdata(cfgtemp,EEG{ineuron});
        
        if length(LF_Vm{ineuron}.trial{1})>length(LF_EEG{ineuron}.trial{1})
            LF_Vm{ineuron}.trial{1}(end)=[];
        end
        
        clear t_1 t_2 t_sel
        
        %make xcorr
        [r,lags]=xcorr(LF_EEG{ineuron}.trial{1},LF_Vm{ineuron}.trial{1},'normalized');
        
        %plot xcorr
        fig_xcorr=figure;
        plot(lags/LF_Vm{ineuron}.fsample,r)
        xline(0);
        ylim([-0.1,0.1])
        xcorrpath=fullfile(config{ineuron}.imagesavedir,'cross_corr');
        
        if~isfolder(xcorrpath)
            mkdir(xcorrpath);
        end
        
        fname_xcorr=fullfile(xcorrpath,sprintf('%s_Low_freq',config{ineuron}.prefix));
        dtx_savefigure(fig_xcorr,fname_xcorr,'pdf','png','close');
        clear fname_xcorr
        %plot LFP and Vm
        
        %select part
        t_1= Events{ineuron}.markers.VentOff.synctime;
        t_2= stats_intra{ineuron}.iso_time;
        t_sel=[t_1 t_2];
        
        cfgtemp= [];
        cfgtemp.latency=t_sel;
        toplot_Vm=ft_selectdata(cfgtemp,Vm{ineuron});
        toplot_LFP=ft_selectdata(cfgtemp,EEG{ineuron});
        
        fig=figure;
        plot(toplot_Vm.time{1},toplot_Vm.trial{1})
        hold
        plot(toplot_LFP.time{1},toplot_LFP.trial{1}*10,'Color','r')
        fname=fullfile(config{ineuron}.imagesavedir,'EEG_Intra','till_Iso',sprintf('%s_til_Iso',config{ineuron}.prefix));
        
        if ~isfolder(fullfile(config{ineuron}.imagesavedir,'EEG_Intra','till_Iso'))
            mkdir(fullfile(config{ineuron}.imagesavedir,'EEG_Intra','till_Iso'));
        end
        
        dtx_savefigure(fig,fname,'pdf','png','close');
        
    end %ineuron
end %slurm

end %function