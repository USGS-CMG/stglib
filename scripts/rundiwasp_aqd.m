rootdir = '/Volumes/Backstaff/field/gb_proc/';
% mooring = '1076';
% dep = 'a';
% mooring = '1076';
% dep = 'b';
% mooring = '1078';
% dep = 'a';
% mooring = '1078';
% dep = 'b';
mooring = '1079';
dep = 'b';

height = '1';

infile = [rootdir mooring dep '/' mooring height upper(dep) 'aqd/' mooring height upper(dep) 'aqdwvsb-cal.nc'];

initial_instrument_height = ncreadatt(infile, '/', 'initial_instrument_height');
nperburst = ncreadatt(infile, '/', 'WaveNumberOfSamples');
fs = ncreadatt(infile, '/', 'WaveSampleRate');
fs = str2num(fs(1:strfind(fs, ' ')-1));
pres = ncread(infile, 'P_1ac'); % this assumes atmospherically corrected pressure
vel1= ncread(infile, 'vel1_1277'); % this assumes atmospherically corrected pressure
vel2 = ncread(infile, 'vel2_1278'); % this assumes atmospherically corrected pressure
vel3 = ncread(infile, 'vel3_1279'); % this assumes atmospherically corrected pressure
cellpos = ncread(infile, 'cellpos'); % this assumes atmospherically corrected pressure
adcpheight = ncreadatt(infile, '/', 'initial_instrument_height');
heading = ncread(infile, 'Hdg_1215');
pitch = ncread(infile, 'Ptch_1216');
roll = ncread(infile, 'Roll_1217'); 

%%
addpath /Users/dnowacki/Documents/matlabdjn/diwasp_1_1GD
%% Process AQD wave data with DIWASP

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

for burst = 1:size(pres,2)
    ID.depth = mean(pres(:,burst)) + adcpheight;
    ID.data = [pres(:,burst) vel1(:,burst)/1000 vel2(:,burst)/1000 vel3(:,burst)/1000];
    ID.layout = make_xyzpos(0, heading(burst), pitch(burst), roll(burst), cellpos(burst), adcpheight)'; % magvar has already been applied
    ID.datatypes={'pres' 'radial' 'radial' 'radial'};
    
%     SM.freqs = 1/256:1/128:ID.fs/2-1/256;
%     SM.dirs = 5:10:360-5;

%     SM.funit = 'Hz';
%     SM.dunit = 'naut';
    
    
    % 'DFTM' Direct Fourier transform method
    % 'EMLM' Extended maximum likelihood method
    % 'IMLM' Iterated maximum likelihood method
    % 'EMEP' Extended maximum entropy principle
    % 'BDM' Bayesian direct method
    %EP.nfft = 8;
    %EP.dres = 10;
    %EP.smooth = 'off';
    
    [diwasp.S(burst), diwasp.E(burst)] = dirspec(ID, SM, EP, {'MESSAGE', 0, 'PLOTTYPE', 0});
%     [diwasp.Hs(burst), diwasp.Tp(burst), diwasp.Dtp(burst), diwasp.Dp(burst)] = infospec(diwasp.S(burst));
    [diwasp.H(burst),diwasp.HsConf(burst,:),diwasp.Tp(burst),diwasp.DTp(burst),diwasp.Dp(burst)] = infospec(diwasp.S(burst));
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
    dw.dspec(:,:,n) = diwasp.S(n).S;
    dw.pspec(:,n) = sum(diwasp.S(n).S, 2) * diff(diwasp.S(1).dirs(1:2));
    m0 = sum(sum(diwasp.S(n).S)) * diff(diwasp.S(1).dirs(1:2)) * diff(diwasp.S(1).freqs(1:2));
    m1 = sum(sum(repmat(diwasp.S(1).freqs', 1, 181) .* diwasp.S(n).S)) * diff(diwasp.S(1).dirs(1:2)) * diff(diwasp.S(1).freqs(1:2));
    dw.wp_4060(n) = m0/m1;
end
%%

outfile = [rootdir mooring dep '/' mooring height upper(dep) 'aqd/' mooring height upper(dep) 'aqdwvs-diwasp.nc'];

nccreate(outfile, 'wh_4061', 'dimensions', {'time', size(dw.wh_4061, 2)});
ncwrite(outfile, 'wh_4061', dw.wh_4061);

nccreate(outfile, 'wp_peak', 'dimensions', {'time', size(dw.wp_peak, 2)});
ncwrite(outfile, 'wp_peak', dw.wp_peak);

nccreate(outfile, 'wp_4060', 'dimensions', {'time', size(dw.wp_4060, 2)});
ncwrite(outfile, 'wp_4060', dw.wp_4060);

nccreate(outfile, 'wvdir', 'dimensions', {'time', size(dw.wvdir, 2)});
ncwrite(outfile, 'wvdir', dw.wvdir);

nccreate(outfile, 'dwvdir', 'dimensions', {'time', size(dw.dwvdir, 2)});
ncwrite(outfile, 'dwvdir', dw.dwvdir);

% TODO: need to do wd_4062

nccreate(outfile, 'frequency', 'dimensions', {'frequency', size(dw.frequency, 2)});
ncwrite(outfile, 'frequency', dw.frequency);

nccreate(outfile, 'direction', 'dimensions', {'direction', size(dw.direction, 2)});
ncwrite(outfile, 'direction', dw.direction);

nccreate(outfile, 'pspec', 'dimensions', {'frequency', size(dw.pspec, 1), 'time', size(dw.pspec, 2)});
ncwrite(outfile, 'pspec', dw.pspec);

nccreate(outfile, 'dspec', 'dimensions', {'frequency', size(dw.pspec, 1), 'direction', size(dw.direction, 2), 'time', size(dw.pspec, 2)});
ncwrite(outfile, 'dspec', dw.dspec);

save([outfile(1:end-2) 'mat'], 'diwasp')

%% 
function xyzpositions = make_xyzpos(magvar, heading, pitch, roll, height, adcpheight)

xyzpos=ones(3,3);
pos=height*tand(25);

% as x,  y , z
% meas 1
% meas 2
% meas 3
xyzpos(:,1)=[0,pos,height];
xyzpos(:,2)=[pos*cosd(30),-pos*.5,height];
xyzpos(:,3)=[-pos*cosd(30),-pos*.5,height];

% set up the new coordinate transformation matrix
CH = cosd(heading+magvar);
SH = sind(heading+magvar);
CP = cosd(pitch);
SP = sind(pitch);
CR = cosd(-roll);
SR = sind(-roll);

%  let the matrix elements be ( a b c; d e f; g h j);
a = CH.*CR - SH.*SP.*SR;  b = SH.*CP; c = -CH.*SR - SH.*SP.*CR;
d = -SH.*CR - CH.*SP.*SR; e = CH.*CP; f = SH.*SR - CH.*SP.*CR;
g = CP.*SR;              h = SP;     j = CP.*CR;

%transform the original x,y,z positions to the new positions accounting for
%heading, pitch and roll... we also add adcpheight back in

new_xyzpos(1,:)=xyzpos(1,:)*a+xyzpos(2,:)*b+xyzpos(3,:)*c;
new_xyzpos(2,:)=xyzpos(1,:)*d+xyzpos(2,:)*e+xyzpos(3,:)*f;
new_xyzpos(3,:)=xyzpos(1,:)*g+xyzpos(2,:)*h+xyzpos(3,:)*j+adcpheight;
new_xyzpos(4,:)=[0,0,adcpheight];

xyzpositions=new_xyzpos;
xyzpositions = xyzpositions([4, 1:3], :);

end