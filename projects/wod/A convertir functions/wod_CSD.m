





itrial=1;
starttrial              = LFP.trialinfo.begsample / LFP.fsample;
endtrial                = LFP.trialinfo.endsample / LFP.fsample;
offsettrial             = LFP.trialinfo.offset / LFP.fsample;


%prepare
t_1= freq_data;
t_end= t_voff+10;


cfgtemp=[];
cfgtemp.latency=[t_voff t_end];
LFP_select= ft_selectdata(cfgtemp,LFP);


cfgtemp=[];
cfgtemp.resamplefs=100;
LFP_select=ft_resampledata(cfgtemp,LFP_select);

cfgtemp=[];

cfgtemp.hpfilter = 'yes';
cfgtemp.hpfreq= 40;
cfgtemp.hpfilttype= 'fir';
LFP_hpfilt=ft_preprocessing(cfgtemp,LFP_select);



clear LFP_CSD
for ichan=1:size(LFP_hpfilt.label,1)
    LFP_CSD(:,ichan)=LFP_hpfilt.trial{1}(ichan,:);
    
end %ichan
LFP_CSD(:,:)= LFP_CSD(:,:)*1E-6;



CSDOutput= CSD(LFP_CSD,LFP_hpfilt.fsample,1E-4);
CSDOutput=CSDOutput(:,2:end-1);

%% Interp


data = squeeze(CSDOutput'); %freq x time
%remove nans
data(:, any(isnan(data),1)) = [];

datax = nan(size(data));
datay = nan(size(data));
for itime = 1:size(data, 2)
    datax(:, itime) = LFP_hpfilt.time{1}(itime);
end
for ichan = 1:size(data, 1)
    data(ichan, :) = movmean(data(ichan, :), 50); %smooth data of each channel : adapter la valeur de smoothing
    datay(ichan, :) = config32{6}.LFP.chan_depth{ichan};
end

dat = reshape(data, [], 1); %colonnes concaténées
chanX = reshape(datax, [], 1); %colonnes concaténées
chanY = reshape(datay, [], 1); %colonnes concaténées

%limits of the interpolation
hlim = [min(chanX) max(chanX)];
vlim = [min(chanY) max(chanY)];

% values of the interpolation
xi = linspace(hlim(1), hlim(2), size(datax,2)); %no interpolation on X data : same nr of points than the original data
yi = linspace(vlim(1), vlim(2), 1000); %Nr of Y samples. 1xN       

%interpolate
[Xi, Yi, Zi] = griddata(chanX', chanY, dat, xi', yi, 'v4'); 

%plot
figure; hold on;
h = surf(Xi, Yi, zeros(size(Zi)), Zi, 'EdgeColor', 'none', 'FaceColor', 'flat');
view(0,90);
colormap jet
contour(Xi, Yi, Zi, 'k');
% contour(Xi, Yi, Zi, 4, 'k'); %si besoin d'ajuster le nombre de lignes du contour

figure;
imagesc(Zi);







% ft_progress('init', 'text');
% for isample=1:size(CSDOutput,1)
%     ft_progress(0, '%d\n',  isample);
%     smoothed(isample, :)=ksdensity(CSDOutput(isample, :),'Bandwidth',0.1, 'NumPoints', 1000); %'Kernel','triangle');
% end
% ft_progress('close');
% 
% 
% fig = figure; hold on;
% imagesc(smoothed');
% colormap jet
% colorbar
% caxis([0 1]);
% min(min(-smoothed))

% smoothed_temp = {};
% ft_progress('init', 'text');
% for isample=1:1000%size(CSDOutput,1)
%     ft_progress(0, '%d\n',  isample);
%     evalc('smoothed_temp{isample} = pchip(1:size(CSDOutput,2),CSDOutput(isample,:),1:1000)');
% end
% ft_progress('close');
% smoothed = cat(1,smoothed_temp{:});
% 
% fig = figure; hold on;
% imagesc(smoothed');
% colormap jet


fig = figure; hold on;
imagesc(-CSDOutput');
contourf(-CSDOutput');
contourf(-CSDOutput','LineStyle','none');
colormap jet
% ft_pimpplot(fig, jet, true);

xticklabels(floor(xticks/LFP_hpfilt.fsample));