infile = '/Volumes/Backstaff/field/gb_proc/1076b/10763Bdw/10763Bdwb-cal.nc';

initial_instrument_height = ncreadatt(infile, '/', 'initial_instrument_height');
fs = 1/ncreadatt(infile, '/', 'sample_interval');
pres = ncread(infile, 'P_1ac'); % this assumes atmospherically corrected pressure

%%
for burst = 1:size(pres, 2)
    ID.data = pres(:,burst);
    ID.layout = [0
                 0
                 initial_instrument_height];
    ID.datatypes = {'pres'};
    
    ID.depth = mean(ID.data);
    ID.fs = fs;
    
    SM.freqs = 1/256:1/128:ID.fs/2-1/256;
    SM.dirs = 22.5:45:360-22.5;
    SM.xaxisdir = 90;
    SM.funit = 'Hz';
    SM.dunit = 'naut';
    
    EP.method = 'IMLM';
    % 'DFTM' Direct Fourier transform method
    % 'EMLM' Extended maximum likelihood method
    % 'IMLM' Iterated maximum likelihood method
    % 'EMEP' Extended maximum entropy principle
    % 'BDM' Bayesian direct method
    
    [diwasp.S(burst), diwasp.E(burst)] = dirspec(ID, SM, EP, {'MESSAGE', 0, 'PLOTTYPE', 0});
    [diwasp.Hs(burst), diwasp.Tp(burst), diwasp.Dtp(burst), diwasp.Dp(burst)] = infospec(diwasp.S(burst));
    disp(num2str(burst))
end
%%
dw.wh_4061 = diwasp.Hs;
dw.wp_peak = diwasp.Tp;
dw.frequency = diwasp.S(1).freqs;
for n = 1:length(diwasp.S)
    dw.pspec(:,n) = sum(diwasp.S(n).S, 2) * diff(diwasp.S(1).dirs(1:2));
    m0 = sum(sum(diwasp.S(n).S)) * diff(diwasp.S(1).dirs(1:2)) * diff(diwasp.S(1).freqs(1:2));
    m1 = sum(sum(repmat(diwasp.S(1).freqs', 1, 8) .* diwasp.S(n).S)) * diff(diwasp.S(1).dirs(1:2)) * diff(diwasp.S(1).freqs(1:2));
    dw.wp_4060(n) = m0/m1;
end
%%
rootdir = '/Volumes/Backstaff/field/gb_proc/';
% % mooring = '1076';
% % dep = 'b';
% mooring = '1077';
% dep = 'b';
% mooring = '1078';
% dep = 'a';
% mooring = '1078';
% dep = 'b';
% mooring = '1079';
% dep = 'a';

mooring = '1080';
dep = 'b';

if strcmp(mooring, '1080')
    height = '1';
else
    height = '3';
end

load([rootdir mooring dep '/' mooring height upper(dep) 'dw/' mooring '2' upper(dep) 'dws-a.mat'])
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