rootdir = '/Volumes/Backstaff/field/gb_proc/';

% mooring = '1076';
% dep = 'b';
% mooring = '1077';
% dep = 'b';
mooring = '1078';
dep = 'a';
% mooring = '1078';
% dep = 'b';
% mooring = '1079';
% dep = 'a';

% mooring = '1080';
% dep = 'b';


if strcmp(mooring, '1080')
    height = '1';
else
    height = '3';
end

infile = [rootdir mooring dep '/' mooring height upper(dep) 'dw/' mooring height upper(dep) 'dwb-cal.nc'];

initial_instrument_height = ncreadatt(infile, '/', 'initial_instrument_height');
fs = 1/ncreadatt(infile, '/', 'sample_interval');
nperburst = ncreadatt(infile, '/', 'burst_length');
pres = ncread(infile, 'P_1ac'); % this assumes atmospherically corrected pressure

%%
addpath /Users/dnowacki/Documents/matlabdjn/diwasp_1_1GD
%%

ID.fs = fs;
SM.nperburst = nperburst;
SM.nsegs=16;
SM.nfft = 2^(nextpow2(SM.nperburst/SM.nsegs));
SM.iter = 100; 
SM.dres=180;
SM.nfreqs=SM.nfft/2;
SM.freqs = ID.fs/SM.nfft:ID.fs/SM.nfft:ID.fs/2;
SM.dirs = -180:360/SM.dres:180;
SM.xaxisdir = 90;
EP.method = 'IMLM';

for burst = 1:size(pres, 2)
    ID.data = pres(:,burst);
    ID.depth = mean(ID.data) + initial_instrument_height;
    
    ID.layout = [0
                 0
                 initial_instrument_height];
    ID.datatypes = {'pres'};
    
    
%     ID.fs = fs;
    
%     SM.freqs = 1/256:1/128:ID.fs/2-1/256;
%     SM.dirs = 22.5:45:360-22.5;
%     SM.xaxisdir = 90;
%     SM.funit = 'Hz';
%     SM.dunit = 'naut';
%     
%     EP.method = 'IMLM';
    % 'DFTM' Direct Fourier transform method
    % 'EMLM' Extended maximum likelihood method
    % 'IMLM' Iterated maximum likelihood method
    % 'EMEP' Extended maximum entropy principle
    % 'BDM' Bayesian direct method
    
    [diwasp.S(burst), diwasp.E(burst)] = dirspec(ID, SM, EP, {'MESSAGE', 0, 'PLOTTYPE', 0});
    [diwasp.H(burst),diwasp.HsConf(burst,:),diwasp.Tp(burst),diwasp.DTp(burst),diwasp.Dp(burst)] = infospec(diwasp.S(burst));
%     [diwasp.Hs(burst), diwasp.Tp(burst), diwasp.Dtp(burst), diwasp.Dp(burst)] = infospec(diwasp.S(burst));
    disp(num2str(burst))
end
%%
dw.wh_4061 = diwasp.H;
dw.wp_peak = diwasp.Tp;
dw.frequency = diwasp.S(1).freqs;
dw.direction = diwasp.S(1).dirs;
dw.wvdir = diwasp.DTp;
dw.dwvdir = diwasp.Dp;
for n = 1:length(diwasp.S)
    dw.pspec(:,n) = sum(diwasp.S(n).S, 2) * diff(diwasp.S(1).dirs(1:2));
    m0 = sum(sum(diwasp.S(n).S)) * diff(diwasp.S(1).dirs(1:2)) * diff(diwasp.S(1).freqs(1:2));
    m1 = sum(sum(repmat(diwasp.S(1).freqs', 1, 181) .* diwasp.S(n).S)) * diff(diwasp.S(1).dirs(1:2)) * diff(diwasp.S(1).freqs(1:2));
    dw.wp_4060(n) = m0/m1;
end
%%




% load([rootdir mooring dep '/' mooring height upper(dep) 'dw/' mooring '2' upper(dep) 'dws-a.mat'])
outfile = [rootdir mooring dep '/' mooring height upper(dep) 'dw/' mooring height upper(dep) 'diwasp.nc'];

nccreate(outfile, 'wh_4061', 'dimensions', {'time', size(dw.wh_4061, 2)});
ncwrite(outfile, 'wh_4061', dw.wh_4061);

nccreate(outfile, 'wp_peak', 'dimensions', {'time', size(dw.wp_peak, 2)});
ncwrite(outfile, 'wp_peak', dw.wp_peak);

nccreate(outfile, 'wp_4060', 'dimensions', {'time', size(dw.wp_4060, 2)});
ncwrite(outfile, 'wp_4060', dw.wp_4060);

nccreate(outfile, 'frequency', 'dimensions', {'frequency', size(dw.frequency, 2)});
ncwrite(outfile, 'frequency', dw.frequency);

nccreate(outfile, 'pspec', 'dimensions', {'frequency', size(dw.pspec, 1), 'time', size(dw.pspec, 2)});
ncwrite(outfile, 'pspec', dw.pspec);

save([outfile(1:end-2) 'mat'], 'diwasp')