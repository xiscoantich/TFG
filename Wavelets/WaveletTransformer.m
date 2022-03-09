classdef WaveletTransformer < handle

    properties (Access = private)
        scale
        paramout
        k
    end

    methods (Access = public)

        function [wave] = directTransform(obj,cParams)
            data = cParams.data;
            signal = data.signal;
            mother = data.motherwave;
            dt = data.dt;
            if data.dim == 1
                [wave,~,obj.scale,~,~,obj.paramout, obj.k] = obj.contwt(signal,dt,[],[],[],[],mother,[]);
            else
                [wave,~,obj.scale,~,~, obj.paramout, obj.k] = obj.contwt2(signal,dt,[],[],[],[],mother,[]);
            end
        end
        
        function signal = inverseTransform(obj,cParams)
            data = cParams.data;
            wave = data.wave;
            mother = data.motherwave;
           
            if data.dim == 1
                signal = obj.invcwt(wave, mother, obj.scale, obj.paramout, obj.k);
            else
                signal = obj.invcwt2(wave, mother, obj.scale, obj.paramout, obj.k);
            end
        end
    end
    
    methods (Access = private)
        
        function [wave,period,scale,coi, dj, paramout, k] = contwt(obj,Y,dt,pad,dj,s0,J1,mother,param)
            
            if (nargin < 8) |isempty(param), param = -1; end
            if (nargin < 7) |isempty(mother) mother = -1; end
            if (nargin < 6) |isempty(J1), J1 = -1; end
            if (nargin < 5) | isempty(s0), s0 = -1; end
            if (nargin < 4) | isempty(dj), dj = -1; end
            if (nargin < 3) | isempty(pad), pad = 0; end
            if (nargin < 2)
                error('Must input a vector Y and sampling time DT')
            end
            
            n1 = length(Y);
            
            if (s0 == -1), s0=2*dt;, end
            if (dj == -1), dj = 1./4.;, end
            if (J1 == -1), J1=ceil((log(n1*dt/s0)/log(2))/dj);, end  %changed fix to ceil(), JE Oct 12 2014
            if (mother == -1), mother = 'MORLET'; end
            
            %....construct time series to analyze, pad if necessary
            x(1:n1) = Y - mean(Y);
            %x(1:n1) = Y;
            
            if (pad == 1)
                x = obj.makepowerof2(x);
            end
            n = length(x);
            
            %....construct wavenumber array used in transform [Eqn(5)]
            k = [1:fix(n/2)];
            k = k.*((2.*pi)/(n*dt));
            k = [0., k, -k(fix((n-1)/2):-1:1)];
            
            %....compute FFT of the (padded) time series
            
            f = fft(x); % [Eqn(3)]
            
            %....construct SCALE array & empty PERIOD & WAVE arrays
            scale = s0*2.^((0:J1)*dj);
            period = scale;
            wave = zeros(J1+1,n);  % define the wavelet array
            wave = wave + i*wave;  % make it complex
            
            % loop through all scales and compute transform
            for a1 = 1:J1+1
                [daughter,fourier_factor,coi,dofmin, paramout]=wave_bases(mother,k,scale(a1),param);
                wave(a1,:) = ifft(f.*daughter);  % wavelet transform[Eqn(4)]
            end
            
            period = fourier_factor*scale;
            coi = coi*dt*[1E-5,1:((n1+1)/2-1),fliplr((1:(n1/2-1))),1E-5];  % COI [Sec.3g]
            wave = wave(:,1:n1);  % get rid of padding before returning
            return
        end
        
        function Xrec = invcwt(obj,wvcfs, mother, scale, param, k)
            
            % take the real part of the wavelet transform.
            Wr = real(wvcfs);
            N = size(Wr, 2);
            %compute the sum
            scale = scale(:); %make a column vector
            s = repmat(scale, [1, size(wvcfs,2)]);
            
            summand = sum(Wr./sqrt(s), 1);
            
            %compute the constant factor out front
            
            %compute the fourier spectrum at each scale [Eq(12)]
            for a1 = 1:length(scale)
                daughter=wave_bases(mother,k,scale(a1),param);
                %     'using conj'
                %     Wdelta(a1) = (1/N)*sum(conj(daughter));  %  [Eqn(12)]
                %
                Wdelta(a1) = (1/N)*sum(daughter);  %  [Eqn(12)]
            end
            % whos Wdelta
            RealWdelta = real(Wdelta);
            RealWdelta = RealWdelta(:); %make a column vector
            
            C = sum(RealWdelta./sqrt(scale));
            
            Xrec = (1/C)*summand;
        end
        
        function [y] = makepowerof2(x)
            N = length(x);
            y = x;
            while mod(log(N)/log(2),1)~=0
                y(N+1) = 0;
                N = N+1;
            end
        end

    end
end