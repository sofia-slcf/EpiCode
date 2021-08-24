starttrial              = LFP_cleaned.trialinfo.begsample / LFP_cleaned.fsample;
endtrial                = LFP_cleaned.trialinfo.endsample / LFP_cleaned.fsample;
offsettrial             = LFP_cleaned.trialinfo.offset / LFP_cleaned.fsample;


cfgtemp=[];
cfgtemp.lpfilter='yes';
cfgtemp.lpfreq=10;
LFP_cleaned=ft_preprocessing(cfgtemp,LFP_cleaned);

t=LFP_cleaned.time{itrial};
voff= MuseStruct{1}{1}.markers.Vent_Off.synctime(itrial);


t_1=t > min(t)/2; 
t_2=t < voff - starttrial(itrial) + offsettrial(itrial);
t_sel=t_1 & t_2;

plot(LFP_cleaned.time{itrial},LFP_cleaned.trial{itrial}(ichan,:))
[uplimit dnlimit]=envelope(LFP_cleaned.trial{itrial}(ichan,t_sel));
hold
plot(LFP_cleaned.time{itrial}(t_sel),uplimit)

Amp_baseline= abs(uplimit-dnlimit);
mean_Ampbaseline=mean(Amp_signal);
threshold=0.1*mean_Ampbaseline;

clear uplimit dnlimit


[up lo]=envelope(LFP_cleaned.trial{itrial}(ichan,:));
Amp_signal= up-lo;


downsample

plot(t,Amp_signal);
yline(threshold);

LFP_cleaned.fsample



for i=1:size(Amp_signal,2)
   
    if i==size(Amp_signal,2)
        change(i)=1;
    end
    change(i)=(Amp_signal(i+1)-mean_Ampbaseline)/mean_Ampbaseline;

end

decrease=change==-0.9;

plot(t,change);

