function wvf=load_wavelet(wavelet,normw)
%Loads definition and properties of a wavelet filter
%wvf=load_wavelet(wavelet,normw)
%
%Input: 
% wavelet - wavelet identification string
% normw - specifies DC and Nyquist normalisation of the wavelet; or can
%         be used to specify type of normalisation:
%         'E' - equal gains of the synthesis filters, with the constraint 
%         that the low-pass analysis filter gain is sqrt(2)
%         'Eu' -  equal gains of the synthesis filters, equal to 1 
%         'V' - equalises error between even and odd reconstructed samples
%               (e,o) subsampling lattice assumed
%         If omitted, the default values are used - sqrt(2)
%
%Output: 
% wvf - structure containing properties of a wavelet:
%  'id' - identifier, i.e. a string descriptor of a wavelet 
%  'wvf_type' - symmetry type of a wavelet, can be 'symmetric_odd',
%   'symmetric_even', 'non_symmetric'. 'odd' specifies an odd number of
%   coefficients in a filter, and 'even' an even number
%  'lift_coeff' - lifting coefficients. If wavelet is not symmetric, there
%   can be two different coefficients per lifting step.
%  'lift_cnct' - which of the two neighbouring pixels are used in the 
%   lifting step.
%  'lift_norm' - normalisation factors used in lifting to achieve the 
%   requested (or default of sqrt(2)) DC and Nyquist gains.
%  'filt_H0','filt_H0_delay' - analysis low-pass filter, with delays
%  'filt_H1','filt_H1_delay' - analysis high-pass filter, with delays
%  'filt_G0','filt_G0_delay' - synthesis low-pass filter, with delays
%  'filt_G1','filt_G1_delay' - synthesis high-pass filter, with delays
%
%Note:
% Decomposition <=> Analysis
% Reconstruction <=> Synthesis
% It is assumed that DC and Nyquist gain factor of the wavelet being 
% retrieved equals sqrt(2). Lifting and filterbank coefficients must 
% correspond in that regard.
% Wavelets specifications are stored in the directory .\Wavelets, each in 
% a human-readable .wvf file. The syntax of a .wvf file is as follows:
%   
% %multiple-line comment�
% [empty line]
% %comment(wvf_type)  
%  wavelet symmetry type
% %comment(lift_coef)
%  1.prediction step = left lifting coecfficient, right lifting coefficient
%  1.update step = left lift. coecff., right lift. coeff.
%  2.prediction step = left lift. coecff., right lift. coeff.
%  2.update step = left lift. coecff., right lift. coeff.
%  ...
% %comment(filt_H0)
%  coefficients of the low-pass analysis filter, one per line - delay,value��
% %comment(filt_H1)
%  coefficients of the high-pass analysis filter, one per line - delay,value��
% %comment(filt_G0 and filt_G1 are not specified, derived from the analysis pair))
%
% The zero-delay refers for the low-pass filter to the even samples, while 
% for the high-pass filter to the odd samples. In other words, 0 delay for
% filt_H1 corresponds to the sample one after the 0-delay sample for the 
% filt_H0. 
% The commonly used (e,o) downsamling lattice is employed.
%
%Example:
% wvf = load_wavelet('Haar');
% wvf = load_wavelet('CDF_9x7',[1 1]);
% wvf = load_wavelet('LeGall_5x3','E');

wvf = struct('id',{},'wvf_type',{},'lift_coeff',{},'lift_cnct',{},'lift_norm',{},...
    'filt_H0',{},'filt_H0_delay',{},'filt_H1',{},'filt_H1_delay',{},...
    'filt_G0',{},'filt_G0_delay',{},'filt_G1',{},'filt_G1_delay',{});
wvf(1).id = wavelet;

LoD_F = []; LoD_F_delay = [];
HiD_F = []; HiD_F_delay = [];
%LoR_F = []; LoR_F_delay = [];
%HiR_F = []; HiR_F_delay = [];

%sets the path where the wavelets are stored
path0=fileparts(which(mfilename));
waveletfile=[path0 '\Wavelets_Data\' wavelet '.wvf'];
fid=fopen(waveletfile,'r');
if fid~=-1 %read the wavelet from .wvf file
    numpar=0;cnt=1;
    while 1 %first skip the starting comments
        tline=fgetl(fid);
        if (isempty(tline) || (tline(1)~='%')) 
            break;
        end;
    end;
    while 1 %the rest are the parameters
        tline=fgetl(fid);
        if ~ischar(tline) 
            break;
        end; %end of file
        if (~isempty(tline)) && (tline(1)~='%')
            switch numpar
                case 1
                  wvf.wvf_type = lower(strtrim(tline));
                case 2
                  tmp = sscanf(tline, '%f,%f');
                  wvf.lift_coeff(cnt,1) = tmp(1);  
                  wvf.lift_coeff(cnt,2) = tmp(2);
                  leftcnct = (wvf.lift_coeff(cnt,1) ~= 0.0);
                  rghtcnct = (wvf.lift_coeff(cnt,2) ~= 0.0);                  
                  wvf.lift_cnct(cnt,:) = sprintf('%c%c',leftcnct+48, rghtcnct+48);                
                case 3
                  wvf.lift_norm(cnt) = str2double(tline);                
                case 4
                  tmp = sscanf(tline, '%f,%f');
                  LoD_F_delay(cnt) = tmp(1);
                  LoD_F(cnt) = tmp(2);
                case 5
                  tmp = sscanf(tline, '%f,%f');
                  HiD_F_delay(cnt) = tmp(1);
                  HiD_F(cnt) = tmp(2);
            end;
            cnt=cnt+1;
        else 
            cnt=1;
            numpar=numpar+1;
        end;
    end;
    fclose(fid);
else        
    error('The specified wavelet cannot be found!');
end;

% if (wvf.lift_coeff(1) == 0)
%     warning('Lifting coefficients not specified!');
% end;
%generate synthesis (reconstruction) from analysis (decomposition) filters
%compute H0(z) and H0(-z)
H0delay = LoD_F_delay;
H0z=LoD_F;
H0nz=LoD_F;
H0odd = mod(H0delay,2) == 1;
H0nz(H0odd) = -H0nz(H0odd);
%compute H1(z) and H1(-z)
H1delay = HiD_F_delay;
H1z=HiD_F;
H1nz=HiD_F;
H1odd = mod(H1delay,2) == 1;
H1nz(H1odd) = -H1nz(H1odd);
%compute P0(z)
[P0z,P0zdelay] = signal_mult(H0z,H0delay,H1nz,H1delay);
[P0nz,P0nzdelay] = signal_mult(H0nz,H0delay,H1z,H1delay);
[PR,PRdelay] = signal_add(P0z,P0zdelay,-1*P0nz,P0nzdelay);

%INSTRUCTION: change PR_tolerance depending on how precise the coefficients have to be
PR_tolerance = 10^(-5);
PR(abs(PR) < PR_tolerance) = 0;
[PR,PRdelay] = signal_add(PR,PRdelay,0,0); %to get rid of extra zeroes

if (length(PR) > 1)
    error('Perfect reconstruction not satisfied!');
end;
if (PR > 0)
    %high-pass synthesis
    HiR_F = LoD_F;
    HiR_F_delay = LoD_F_delay - PRdelay;
    HiReven = mod(LoD_F_delay,2) == 0;
    HiR_F(HiReven) = -HiR_F(HiReven);
    %low-pass synthesis
    LoR_F = HiD_F;
    LoR_F_delay = HiD_F_delay - PRdelay;
    LoRodd = mod(HiD_F_delay,2) == 1;
    LoR_F(LoRodd) = -LoR_F(LoRodd);
else
    %high-pass synthesis
    HiR_F = LoD_F;
    HiR_F_delay = LoD_F_delay - PRdelay + 1; %+1 is to compensate for the way the transform is implemented (high pass samples taken at odd positions);
    HiRodd = mod(LoD_F_delay,2) == 1;
    HiR_F(HiRodd) = -HiR_F(HiRodd);
    %low-pass synthesis
    LoR_F = HiD_F;
    LoR_F_delay = HiD_F_delay - PRdelay;
    LoReven = mod(HiD_F_delay,2) == 0;
    LoR_F(LoReven) = -LoR_F(LoReven);    
end;
%INSTRUCTION: change norm_tolerance depending on how precise the coefficients have to be
norm_tolerance = 10^(-6);
% if ~isempty(LoD_F)
%     if abs(sum(LoD_F) + sqrt(2)) < norm_tolerance
%         %to avoid sign reversal, make sure that the central coefficients are positive!
%         warning('Low-pass analysis wavelet is inverted (sign reversed)!');
%     else
%         if abs(sum(LoD_F) - sqrt(2)) > norm_tolerance
%             warning('The DC norm of LoD_F differs from sqrt(2) more than the predefined tolerance allows!');
%         end;
%     end;
% end;
% 
% if ~isempty(LoR_F)
%     if abs(sum(LoR_F) + sqrt(2)) < norm_tolerance
%         %to avoid sign reversal, make sure that the central coefficients are positive!
%         warning('Low-pass synthesis wavelet is inverted (sign reversed)!');
%     else
%         if abs(sum(LoR_F) - sqrt(2)) > norm_tolerance
%             warning('The DC norm of LoR_F differs from sqrt(2) more than the predefined tolerance allows!');
%         end;
%     end;
% end;

%Normalisation specification
%sqrt(2) for DC and Nyquist normalisation, by default
lpen = sqrt(2)/sum(LoD_F);
hpen = sqrt(2)/sum(LoR_F);
if nargin == 2
    if (upper(normw) == 'E') %error energy equalised normalisation
        %gains adapted so the analysis DC gain is sqrt(2)
        lpen = sqrt(2)/sum(LoD_F);
        hpen = sqrt(sum((lpen*HiR_F).^2)/sum(LoR_F.^2));
    elseif strcmp(upper(normw),'EU') %gains are equal to 1
        hpen = 1/sqrt(sum(LoR_F.^2));
        lpen = 1/sqrt(sum(HiR_F.^2));
    elseif all(upper(normw) == 'V') %equalises error on odd/even samples
        %even samples
        emse_Lo = sum(LoR_F(~LoRodd).^2);
        emse_Hi = sum(HiR_F(HiRodd).^2);
        %odd samples
        omse_Lo = sum(LoR_F(LoRodd).^2);
        omse_Hi = sum(HiR_F(~HiRodd).^2);
        Lodiff = emse_Lo - omse_Lo;
        Hidiff = omse_Hi - emse_Hi;
        %if the condition below is not satisfied the required normalisation
        %has been already applied, or it cannot be applied
        % no warning message produced!
        if ~((abs(Hidiff) <  norm_tolerance) || (abs(Lodiff) <  norm_tolerance))
            %gains adapted so the analysis DC gain is sqrt(2)
            r =  sqrt(Hidiff/Lodiff);
            lpen = 1;
            hpen = r;
            %minimisation with respect to min(2 - (Lomse + Himse))
            %r =  Hidiff/Lodiff;
            %Lomse = emse_Lo + omse_Lo;
            %Himse = omse_Hi + emse_Hi;
            %B = (r * Lomse + Himse) / (r^2 * Lomse^2 + Himse^2);
            %A = r * B;
            %lpen = 1/sqrt(A);
            %hpen = 1/sqrt(B);
            %ad-hoc
            %r =  sqrt(Hidiff/Lodiff);
            %lpen = 1/sqrt(r);
            %hpen = sqrt(r);
        end;
    elseif ~ischar(normw)
        if (length(normw) ~= 2)
            error('Normalisation not specified for both low-pass and high-pass filters!');
        end;
        lpen = normw(1)/sum(LoD_F);
        hpen = normw(2)/sum(LoR_F);
    end;
end;
LoD_F = LoD_F * lpen;
HiR_F = HiR_F * lpen;
HiD_F = HiD_F * hpen;
LoR_F = LoR_F * hpen;
wvf(1).lift_norm(1) = wvf(1).lift_norm(1) * lpen;
wvf(1).lift_norm(2) = wvf(1).lift_norm(2) * hpen;
        
wvf(1).filt_H0 = LoD_F;
wvf(1).filt_H0_delay = LoD_F_delay;
wvf(1).filt_H1 = HiD_F;
wvf(1).filt_H1_delay = HiD_F_delay;
wvf(1).filt_G0 = LoR_F;
wvf(1).filt_G0_delay = LoR_F_delay;
wvf(1).filt_G1 = HiR_F;
wvf(1).filt_G1_delay = HiR_F_delay;

function [signal,delay]=signal_add(s1,d1,s2,d2)
if (any(s2))
    mind = min(d1(1),d2(1));
    maxd = max(d1(end),d2(end));
    dtot = mind:maxd;
    totlen = maxd - mind + 1;
    s1ext = zeros(1,totlen);
    s2ext = zeros(1,totlen);
    s1ext((d1(1) - mind + 1):(d1(end) - mind + 1)) = s1;
    s2ext((d2(1) - mind + 1):(d2(end) - mind + 1)) = s2;
    sigcum = s1ext + s2ext;
else
    mind = d1(1);
    sigcum = s1;
end;
nz = find(sigcum ~= 0);
signal = sigcum(nz(1):nz(end));
delay = (mind + nz(1) - 1):(mind + nz(end) - 1);

function [signal,delay]=signal_mult(s1,d1,s2,d2)
mind = min(d1(1),d2(1));
maxd = max(d1(end),d2(end));
dtot = mind:maxd;
totlen = maxd - mind + 1;
s1ext = zeros(1,totlen);
s2ext = zeros(1,totlen);
s1ext((d1(1) - mind + 1):(d1(end) - mind + 1)) = s1;
s2ext((d2(1) - mind + 1):(d2(end) - mind + 1)) = s2;

mindelay = 2 * mind; %d1(1) + d2(1);
maxdelay = 2 * maxd; %d1(end) + d2(end);
sigcum = zeros(1, maxdelay - mindelay + 1);
for i=1:totlen
    scomp = s1ext * s2ext(i);
    sigind = (1:totlen) + i - 1;
    sigcum(sigind) = sigcum(sigind) + scomp;
end;
nz = find(sigcum ~= 0);
signal = sigcum(nz(1):nz(end));
delay = (mindelay + nz(1) - 1):(mindelay + nz(end) - 1);
