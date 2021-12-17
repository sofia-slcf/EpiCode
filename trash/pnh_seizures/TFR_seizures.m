function [TFR] = TFR_seizures(cfg, Trialdata, force)

% This file is part of EpiCode, see
% http://www.github.com/stephenwhitmarsh/EpiCode for documentation and details.
%
%    EpiCode is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.
%
%    EpiCode is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with EpiCode. If not, see <http://www.gnu.org/licenses/>.

fname_out = fullfile(cfg.datasavedir,[cfg.prefix,'TFRtrials.mat']);

if exist(fname_out,'file') && force == false
    fprintf('************************************\n');
    fprintf('*** Loading precomputed TFR data ***\n');
    fprintf('************************************\n\n');
    load(fname_out,'TFR');
    return
end

fprintf('********************************\n');
fprintf('*** (re-) computing TFR data ***\n');
fprintf('********************************\n\n');

for ipart = 1 : size(Trialdata,2)
    
    for markername = string(fields(Trialdata{ipart}))'
        
        cfgtemp                         = [];
        cfgtemp.channel                 = 'all'; 
        cfgtemp.method                  = 'mtmconvol';
        cfgtemp.output                  = 'pow';
        cfgtemp.taper                   = 'dpss';
        cfgtemp.pad                     = 'nextpow2'; 
        cfgtemp.keeptrials              = 'yes';
        cfgtemp.foi                     = 10:2:200;
        cfgtemp.t_ftimwin               = ones(size(cfgtemp.foi));
        cfgtemp.tapsmofrq               = 5;
        
        cfgtemp.toi                     = cfg.epoch.toi.(markername)(1) : 0.10 : cfg.epoch.toi.(markername)(2);
%         cfgtemp.feedback = 'off';
        TFR{ipart}.(markername)         = ft_freqanalysis(cfgtemp,Trialdata{ipart}.(markername));
        
    end
end % ipart
    

save(fname_out,'TFR','-v7.3');

