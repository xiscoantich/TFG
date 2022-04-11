classdef WaveletTransformer < handle

    properties (Access = public)
        type_wt
        N
        scale
        paramout
        k
        l
        
        packet
        packet_stream
        s
    end

    methods (Access = public)

        function [wave] = directTransform(obj,cParams)
            data = cParams.data;
            %obj.type_wt = data.type.wt;
            if data.dim == 1
                wave = obj.wt1d(data);
            elseif data.dim == 2
                wave = obj.wt2d(data);
            end
        end
        
        function signal = inverseTransform(obj,cParams)
            data = cParams.data;
%             obj.type_wt = data.type.wt;
%             obj.N = data.wave_info.N;
%             obj.scale = data.wave_info.scale;
%             obj.paramout = data.wave_info.paramout;
%             obj.k = data.wave_info.k;
           
            if data.dim == 1
                signal = obj.iwt1d(data);
            elseif data.dim == 2
                signal = obj.iwt2d(data);     
            end
        end
    end
    
    methods (Access = private)
        
        function wave = wt1d(obj,data)
            switch data.type.wt
                case 'cwt'
                    [wave,~,obj.scale,~,~,obj.paramout, obj.k] = obj.contwt(data.signal,data.dt,[],[],[],[],data.motherwave,[]);
                    obj.plotSpectogram(wave)
                case 'convolution'
                    [wave(:,1),wave(:,2)] = obj.dwt_conv1D(data.signal,data.motherwave);
                case 'lifting'
                    [wave(:,1),wave(:,2)] = obj.dwt_lifting1D(data.signal,data.motherwave);
                case 'dyadic_decomp' %Esta no se si la deberia eliminar
                    [wave,obj.N] = obj.dwt_dyadic_decomp(data.signal,data.motherwave,obj.N);
                case 'multilevel'
                    [wave,obj.l] = obj.mlwavelet(data.signal,data.level,data.motherwave);
                case 'matlab'
                    [wave,obj.l] = wavedec(data.signal,data.level,data.motherwave);
            end
        end
        
        function wave = wt2d(obj,data)
            switch data.type.wt
                case 'cwt'
                    [wave,~,obj.scale,~,~, obj.paramout, obj.k] = obj.contwt2(data.signal,data.dt,[],[],[],[],data.motherwave,[]);
                case 'dwt'
                    [A,H,V,D] = obj.dwt_2D(data.signal,data.motherwave);
                    wave(:,:,1) = A;
                    wave(:,:,2) = H;
                    wave(:,:,3) = V;
                    wave(:,:,4) = D;
                case 'multilevel'
                    [wave,obj.l] = obj.mlwavelet2(data.signal,data.level,data.motherwave); %Deberia guardar los datos de otra manera 
                case 'matlab'
                    [wave,obj.l] = wavedec2(data.signal,data.level,data.motherwave);
                case 'packet'
                    [wave,data.wave_info.packet_stream,obj.s,E]=obj.decomp_packets2D(data.signal,data.wave_info.par,data.wave_info.ent_par);
                    obj.draw_packets(wave,data.wave_info.par.N,data.wave_info.par.pdep,obj.s,obj.packet_stream);
            end
        end
        
        function signal = iwt1d(obj,data)
            switch data.type.wt
                case 'cwt'
                    signal = obj.invcwt(data.wave, data.motherwave, data.wave_info.scale, data.wave_info.paramout, data.wave_info.k);
                case 'convolution'
                    signal = obj.idwt_conv1D(data.wave,data.wave_info.d,data.motherwave);
                case 'lifting'
                    signal = obj.idwt_lifting1D(data.wave,data.wave_info.d,data.motherwave);
                case 'dyadic_decomp'
                    signal = obj.idwt_dyadic_recon(data.wave,data.motherwave,data.wave_info.N);
                case 'multilevel'
                    signal = waverec(data.wave,data.wave_info.l,data.motherwave); %Esta funcion es de matlab
                case 'matlab'
                    signal = waverec(data.wave,data.wave_info.l,data.motherwave);
            end
        end
        
        function signal = iwt2d(obj,data)
            switch data.type.wt
                case 'cwt'
                    signal = obj.invcwt2(data.wave, data.motherwave, obj.scale, obj.paramout, obj.k);
                case 'dwt'
                    A = data.wave(:,:,1);
                    H = data.wave(:,:,2);
                    V = data.wave(:,:,3);
                    D = data.wave(:,:,4);
                    signal = obj.idwt_2D(A,H,V,D,data.motherwave);
                case 'multilevel'
                    signal = waverec2(data.wave,data.wave_info.l,data.motherwave); %Esta funcion es de matlab y no funciona bien
                case 'matlab'
                    signal = waverec2(data.wave,data.wave_info.l,data.motherwave);
                case 'packet'
                    signal = recon_packets2D(data.wave,data.wave_info.par,data.wave_info.packet_stream);
                    %signal = obj.recon_packets2D(data.wave,data.wave_info);
            end
        end
        
        function [wave,period,scale,coi, dj, paramout, k] = contwt(obj,Y,dt,pad,dj,s0,J1,mother,param)
            
            if (nargin < 8) || isempty(param), param = -1; end
            if (nargin < 7) || isempty(mother) mother = -1; end
            if (nargin < 6) || isempty(J1), J1 = -1; end
            if (nargin < 5) || isempty(s0), s0 = -1; end
            if (nargin < 4) || isempty(dj), dj = -1; end
            if (nargin < 3) || isempty(pad), pad = 0; end
            if (nargin < 2)
                error('Must input a vector Y and sampling time DT')
            end
            
            n1 = length(Y);
            
            if (s0 == -1), s0=2*dt; end
            if (dj == -1), dj = 1./4.; end
            if (J1 == -1), J1=ceil((log(n1*dt/s0)/log(2))/dj); end  %changed fix to ceil(), JE Oct 12 2014
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
        
        function [wave,period,scale,coi, dj, paramout, k] = contwt2(obj,Y,dt,pad,dj,s0,J1,mother,param)
            m=size(x,2); %columnas
            for i=1:m
                [wave(:,:,i),period,scale,coi, dj, paramout, k] = contwt(obj,Y(:,i),dt,pad,dj,s0,J1,mother,param)
            end
            n=size(x1,1); %filas
            for i=1:n
                y(i,:)=obj.FFT(x1(i,:)); %FFT por filas
            end
        end
        
        function Xrec = invcwt(obj,wvcfs, mother, scale, param, k)
            
            if isempty(mother)
                mother = 'MORLET';
            end
            
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
        
        function [c,l] = mlwavelet(obj,Y,n,wavelet)
            %load the wavelet here
            if ischar(wavelet)
                wvf = obj.load_wavelet(wavelet,'E');
            else
                wvf = wavelet;
            end
            
            if isempty(n) 
                n = wmaxlev(length(Y),wavelet);
            end
            
            Y=double(Y);
            [Dcols]=size(Y,1);
            subc=Dcols/(2^n);%dimensions of the lowest subband
            %check the number of decompositions
            if (round(subc) ~= subc)
                %at the moment, only powers of two supported
                error('Illegal number of decompositions for a given matrix!');
            end
             
            % Initialization.
            s = size(Y);
            Y = Y(:).'; % row vector
            c = [];
            l = zeros(1,n+2,'like',real(Y([])));
            
            l(end) = length(Y);
            for k = 1:n
                [Y,d] = obj.dwt_lifting1D(Y,wvf); % decomposition
                c     = [d c];            % store detail
                l(n+2-k) = length(d);     % store length
            end
            
            % Last approximation.
            c = [Y c];
            l(1) = length(Y);
            
            if s(1)>1
                c = c.';
                l = l';
            end
        end
        
        function [c,s] = mlwavelet2(obj,x,n,wavelet)
            
            %load the wavelet here
            if ischar(wavelet)
                wvf = obj.load_wavelet(wavelet,'E');
            else
                wvf = wavelet;
            end
            
            if isempty(n)
               n = wmaxlev(length(x),wavelet);
            end
            % Initialization.
            c = [];
            sx =  size(x);
            s = zeros(n+2,length(sx));
            if isempty(x)
                c = cast(c,"like",x);
                return;
            end
            
            s(end,:) = size(x);
            for i=1:n  
                [x,h,v,d] = obj.dwt_2D(x,wvf); % decomposition
                c = [h(:)' v(:)' d(:)' c];     % store details
                s(n+2-i,:) = size(x);          % store size
            end
            
            % Last approximation.
            c = [x(:)' c];
            s(1,:) = size(x);
        end

        function [a,d] = dwt_conv1D(obj,x,wvf)
            if ischar(wvf)
                wvf = obj.load_wavelet(wvf,'E');
            end
            sym_ext = false;
            if (strcmp(wvf.wvf_type,'symmetric_even') || strcmp(wvf.wvf_type,'symmetric_odd'))
                sym_ext = true;
            end
            
            %low-pass filtering
            Lle = -1 * wvf.filt_H0_delay(1); %left extension length
            Lre = wvf.filt_H0_delay(end); %right extension length
            if (sym_ext)
                %setup the indices by using the symmetric extension
                I = [(Lle+1):-1:2 1:length(x) (length(x)-1):-1:(length(x) - Lre)];
            else
                %setup the indices by using the periodic extension
                I = [(length(x) - Lle):(length(x)-1) 1:length(x) 1:Lre];
            end
            h0 = fliplr(wvf.filt_H0);
            %oversampled approximation signal
            ao = conv2(x(I),h0','valid');
            %downsample
            a = ao(1:2:end);
            
            %high-pass filtering
            Hle = -1 * wvf.filt_H1_delay(1); %left extension
            Hre = wvf.filt_H1_delay(end); %right extension
            if (sym_ext)
                %setup the indices by using the symmetric extension
                I = [(Hle+1):-1:2 1:length(x) (length(x)-1):-1:(length(x) - Hre)];
            else
                %setup the indices by using the periodic extension
                I = [(length(x) - Hle):(length(x)-1) 1:length(x) 1:Hre];
            end
            h1 = fliplr(wvf.filt_H1);
            %oversampled approximation signal
            do = conv2(x(I),h1','valid');
            %downsample
            d = do(1:2:end);
        end
        
        function y = idwt_conv1D(obj,a,d,wvf)
            if ischar(wvf)
                wvf = obj.load_wavelet(wvf,'E');
            end
            sym_ext = false;
            if (strcmp(wvf.wvf_type,'symmetric_even') || strcmp(wvf.wvf_type,'symmetric_odd'))
                sym_ext = true;
            end
            
            %%low-pass filtering%%
            Lle = -1 * wvf.filt_G0_delay(1); %left extension
            Lre = wvf.filt_G0_delay(end); %right extension
            %upsample the approximation signal
            au = zeros(1, 2*length(a));
            au(1:2:end) = a;
            if (sym_ext)
                %setup the indices by using the symmetric extension
                I = [(Lle+1):-1:2 1:length(au) (length(au)-1):-1:(length(au) - Lre)];
            else
                %setup the indices by using the periodic extension
                I = [(length(au) - Lle):(length(au)-1) 1:length(au) 1:Lre];
            end
            g0 = fliplr(wvf.filt_G0);
            %convolution
            ao = conv2(au(I),g0,'valid');
            
            %%high-pass filtering%%
            Hle = -1 * wvf.filt_G1_delay(1); %left extension
            Hre = wvf.filt_G1_delay(end); %right extension
            %upsample the detail signal
            du = zeros(1, 2*length(d));
            du(2:2:end) = d; %subsampling defined so it starts on odd pixels
            %Note that symmetric extension is here different since the point of
            %symmetry is around odd numbered pixel!
            if (sym_ext)
                %setup the indices by using the symmetric extension
                I = [(Hle+1):-1:2 1:length(du) (length(du)-1):-1:(length(du) - Hre)];
            else
                %setup the indices by using the periodic extension
                I = [(length(du) - Hle):(length(du)-1) 1:length(du) 1:Hre];
            end
            g1 = fliplr(wvf.filt_G1);
            %convolution
            do = conv2(du(I),g1,'valid');
            
            %%combine%%
            y = ao + do;
        end
        
        function [Y,N]=dwt_dyadic_decomp(obj,X,wavelet,N)
            %Dyadic wavelet decomposition of a multidimensional signal
            %[Y,N]=dwt_dyadic_decomp(X,wavelet,N)
            %
            %Input:
            % X - matrix to be transformed containing the input signal or image's
            %     filename (in that case the transform is 2D)
            % wavelet - wavelet identification string
            % N - [optional, default = max] specifies the number of levels of decomposition
            %     default value is the maximum possible, e.g. down to the smallest possible
            %     LL subband
            %
            %Output:
            % Y - matrix of wavelet coefficients
            % N - number of actually performed levels of decomposition
            %
            %Note:
            % Performs dyadic decomposition, i.e. the dimensions of low-pass subband
            % are half of the original subband.
            %
            %Uses:
            % submatrix.m
            % dwt_dim.m
            % subband_dim.m
            %
            %Example:
            % Y=dwt_dyadic_decomp(X,'CDF_9x7',6);
            % [Y,N]=dwt_dyadic_decomp('Lena512.png','Haar');
            
            if ischar(X)
                X=imread(X);
            end
            Y=double(X);
            
            transposed = 0; %by default do not transpose
            if isvector(X)
                n = 1; %number of dimensions
                if (size(X,1) == 1) %if one-row vector, needs to be transposed
                    Y = Y';
                    transposed = 1; %remember to transpose it back later
                end
            else
                n = ndims(Y);
            end
            Xsiz=size(Y);
            
            %to find Nmax - the maximum number of decompositions possible
            [~,hd,Nmax]= obj.subband_dim(Xsiz, Inf);
            if isempty(N)
                %if not specified, the number of decomposition is set to Nmax
                N = Nmax;
            elseif (N > Nmax)
                warning(['Specified number of decompositions exceeds the maximum. N is set to Nmax = ' num2str(Nmax)]);
                N = Nmax;
            end
            
            Lsiz = Xsiz; %low-pass subband dimensions
            for i=1:N
                [Li,Lind] = obj.submatrix(Y,Lsiz);
                for j=1:n %transform in the j-th dimension
                    Li = obj.dwt_dim(Li,j,wavelet);
                end
                Y(Lind{:}) = Li;
                Lsiz = ceil(Lsiz/2); %i.e. low-pass of signal of 3 samples contains 2 samples
            end
            
            if transposed
                Y = Y';
            end
        end
        
        function X=dwt_dim(obj,X,d,wavelet)
            %DWT in specific dimension of an n-dimensional matrix
            %X=dwt_dim(X,d,wavelet)
            %
            %Input:
            % X - matrix to be transformed containing the input signal
            % d - dimension in which the transform will take place
            % wavelet - wavelet identification string, or wavelet data structure
            %
            %Output:
            % X - matrix of wavelet coefficients
            %
            %Note:
            % Filters across the d-th dimension. For instance, if X is 2D then it first
            % filters across columns, as size(X,1) is number of rows (column length).
            %
            %Uses:
            % load_wavelet.m
            % dwt_lifting1D.m
            % subband_dim.m
            %
            %Example:
            % Y = dwt_dim(X,1,'CDF_9x7');
            
            %load the wavelet here
            if ischar(wavelet)
                wvf = obj.load_wavelet(wavelet);
            else
                wvf = wavelet;
            end
            
            if ~isa(X,'double')
                X = double(X);
            end
            N = ndims(X);
            dimprod = numel(X);
            X = shiftdim(X,d-1); %rotates the order of dimensions
            sv = size(X); %matrix size before reshaping
            sizdim = sv(1); %size of the first dimension
            if sizdim > 1 %if non-singleton dimension
                sizcol = dimprod/sizdim; %product of other dimensions
                X = reshape(X,sizdim,sizcol); %reshape into 2D sizdim x sizcol matrix
                lpasiz = obj.subband_dim(sizdim, 1);
                for j=1:sizcol
                    [X(1:lpasiz,j),X(lpasiz+1:sizdim,j)] = obj.dwt_lifting1D(X(:,j),wvf);
                end
                X = reshape(X,sv);
            end;
            X = shiftdim(X,N-d+1); %rotates the order of dimensions forward to the original
        end
        
        function [a,d] = dwt_lifting1D(obj,x,wvf)
            %DWT of a 1D signal in lifting implementation
            %[a,d]=dwt_lifting1D(x,wvf)
            %
            %Input:
            % x - signal to be transformed
            % wvf - wavelet identification string, or wavelet data structure
            %
            %Output:
            % a - approximation (low-pass) signal
            % d - detail (high-pass) signal
            %
            %Uses:
            % load_wavelet.m
            %
            %Example:
            % [a,d] = dwt_lifting1D(x,wvf);
            % [a,d] = dwt_lifting1D(x,'CDF_9x7');
            
            if ischar(wvf)
                wvf = obj.load_wavelet(wvf);
            end;
            s = wvf.lift_coeff;
            K = wvf.lift_norm;
            cn = wvf.lift_cnct;
            
            %xe - for 1-pixel extended signal
            xe = zeros(1,length(x)+2);
            xe(2:end-1) = x;
            if (strcmp(cn,'00'))
                warning('Lifting not available!');
            else
                for i=1:size(s,1)
                    xe(1) = xe(3); %extension on the left
                    xe(end) = xe(end-2); %extension on the right
                    start = rem(i,2); %determines if it is prediction or update step, 1 - prediction, 0 - update
                    lind = 1+start:2:length(xe)-2;
                    cind = lind + 1;
                    rind = cind + 1;
                    if (cn(i,1) == '1') %left connection present
                        xe(cind) = xe(cind) + s(i,1)*xe(lind); %left pixel lifting
                    end;
                    if (cn(i,2) == '1') %right connection present
                        xe(cind) = xe(cind) + s(i,2)*xe(rind); %right pixel lifting
                    end;
                end;
            end;
            %normalisation
            a = xe(2:2:end-1) * K(1);
            d = xe(3:2:end-1) * K(2);
        end
            
        function Xr = idwt_dyadic_recon(obj,Y,wavelet,N)
            %Dyadic wavelet reconstruction of a multidimensional signal
            %X=idwt_dyadic_recon(Y,wavelet,N)
            %
            %Input:
            % X - matrix of wavelet coefficients
            % wavelet - wavelet identification string
            % N - specifies the number of levels of reconstruction (inverse DWT)
            %
            %Output:
            % Xr - reconstructed matrix
            %
            %Uses:
            % submatrix.m
            % idwt_dim.m
            %
            %Example:
            % Xr=idwt_dyadic_recon(Y,'CDF_9x7',6);
            
            Xr=double(Y);
            
            transposed = 0; %by default, do not transpose
            if (size(Y,1) == 1) %if one-row vector, needs to be transposed
                Xr = Xr';
                transposed = 1; %remember to transpose it back later
            end
            Xsiz=size(Xr);
            if isvector(Xr)
                n = 1;
            else
                n = ndims(Xr);
            end
            
            for i=N-1:-1:0
                Lsiz = ceil(Xsiz / 2^i); %low-pass subband dimensions
                [Li,Lind] = obj.submatrix(Xr,Lsiz);
                for j=n:-1:1 %inverse transform in j-th dimension
                    Li = obj.idwt_dim(Li,j,wavelet);
                end
                Xr(Lind{:}) = Li;
            end
            
            if transposed
                Xr = Xr';
            end
        end
        
        function X = idwt_dim(obj,X,d,wavelet)
            %IDWT in specific dimension of an n-dimensional matrix
            %X=idwt_dim(X,d,wavelet)
            %
            %Input:
            % X - matrix of wavelet coeffcients
            % d - dimension in which the transform will take place
            % wavelet - wavelet identification string, or wavelet data structure
            %
            %Output:
            % X - reconstructed matrix
            %
            %Uses:
            % load_wavelet.m
            % idwt_lifting1D.m
            % subband_dim.m
            %
            %Example:
            % X = idwt_dim(Y,1,'Haar');
            
            %load the wavelet here
            if ischar(wavelet)
                wvf = obj.load_wavelet(wavelet,'E');
            else
                wvf = wavelet;
            end;
            
            N = ndims(X);
            dimprod = numel(X);
            X = shiftdim(X,d-1); %rotates the order of dimensions
            sv = size(X); %matrix size before reshaping
            sizdim = sv(1); %size of the first dimension
            if sizdim > 1 %if non-singleton dimension
                sizcol = dimprod/sizdim; %product of other dimensions
                X = reshape(X,sizdim,sizcol); %reshape into 2D sizdim x sizcol matrix
                [lpasiz,hpasiz] = obj.subband_dim(sizdim, 1);
                for j=1:sizcol
                    X(:,j) = obj.idwt_lifting1D(X(1:lpasiz,j),X(lpasiz+1:sizdim,j),wvf);
                end;
                X = reshape(X,sv);
            end;
            X = shiftdim(X,N-d+1); %rotates the order of dimensions forward to the original
        end
        
        function y = idwt_lifting1D(obj,a,d,wvf)
            %IDWT of a 1D signal in lifting implementation
            %y=idwt_lifting1D(a,d,wvf)
            %
            %Input:
            % a - approximation (low-pass) signal
            % d - detail (high-pass) signal
            % wvf - wavelet identification string, or wavelet data structure
            %
            %Output:
            % y - reconstructed signal
            %
            %Uses:
            % load_wavelet.m
            %
            %Example:
            % y = idwt_lifting1D(a,d,wvf);
            % y = idwt_lifting1D(a,d,'CDF_9x7');
            
            if ischar(wvf)
                wvf = obj.load_wavelet(wvf,'E');
            end
            s = wvf.lift_coeff;
            K = wvf.lift_norm;
            cn = wvf.lift_cnct;
            
            %xe - for 1-pixel extended signal
            xe = zeros(1,length(a)+length(d)+2);
            %undo the normalisation
            xe(2:2:end-1) = a / K(1);
            xe(3:2:end-1) = d / K(2);
            if (strcmp(cn,'00'))
                warning('Lifting not available!');
            else
                for i=size(s,1):-1:1
                    xe(1) = xe(3); %extension on the left
                    xe(end) = xe(end-2); %extension on the right
                    start = rem(i,2); %determines if it is prediction or update step, 1 - prediction, 0 - update
                    lind = 1+start:2:length(xe)-2;
                    cind = lind + 1;
                    rind = cind + 1;
                    if (cn(i,1) == '1') %left connection present
                        xe(cind) = xe(cind) - s(i,1)*xe(lind); %left pixel lifting
                    end
                    if (cn(i,2) == '1') %right connection present
                        xe(cind) = xe(cind) - s(i,2)*xe(rind); %right pixel lifting
                    end
                end
            end
            y = xe(2:end-1);
            
        end
        
        function [A,H,V,D] = dwt_2D(obj,X,wavelet)
            %Two-dimensional separable DWT
            %[A,H,V,D]=dwt_2D(X,wavelet)
            %
            %Input:
            % X - matrix to be transformed containing the input signal
            % wavelet - wavelet identification string, or wavelet data structure
            %
            %Output:
            % A,H,V,D - approximation signal, horizontal, vertical and diagonal details
            %           signal
            %
            %Uses:
            % load_wavelet.m
            % dwt_dim.m
            % subband.m
            %
            %Example:
            % [A,H,V,D] = dwt_2D(X,'CDF_9x7');
            
            %load the wavelet here
            if ischar(wavelet)
                wvf = obj.load_wavelet(wavelet);
            else
                wvf = wavelet;
            end
            
            X = obj.dwt_dim(X,1,wvf); %columns
            X = obj.dwt_dim(X,2,wvf); %rows
            
            A = obj.subband(X,1,'ll');
            H = obj.subband(X,1,'hl');
            V = obj.subband(X,1,'lh');
            D = obj.subband(X,1,'hh');
        end
        
        function [S,Sind,Sdim]=subband(obj,D,N,band)
            %[S,Sind,Sdim]=subband(D,N,band)
            %Version: 3.01, Date: 2005/01/01, author: Nikola Sprljan
            %
            %Input:
            % D - array of wavelet coefficients
            % N - specifies the decomposition level
            % band - specifies the subband ('ll', 'hl', 'lh' or 'hh')
            %
            %Output:
            % S - array of the subband wavelet coefficients
            % Sind - indices of elements of the selected subband
            % Sdim - dimensions of the selected subband
            %
            %Note:
            % ('ll', 'hl', 'lh', 'hh') corresponds to ('a', 'h', 'v' ,'d')
            %
            %Uses:
            % subband_dim.m
            %
            %Example:
            % D=dwt_dyadic_decomp(A,'CDF_9x7',4);
            % S=subband(D,4,'hl');
            % S=subband(D,6,'ll');
            
            if (ndims(D) ~= 2)
                error('Wavelet coefficients array of other dimensions than 2D!');
            end;
            
            [sizrow,sizcol] = size(D);
            [ldr,hdr]= obj.subband_dim(sizrow, N);
            [ldc,hdc]= obj.subband_dim(sizcol, N);
            Sind = [];
            switch band
                case 'll'
                    Sind{1} = 1:ldr;
                    Sind{2} = 1:ldc;
                    Sdim = [ldr ldc];
                case 'hl'
                    Sind{1} = 1:ldr;
                    Sind{2} = ldc+1:ldc+hdc;
                    Sdim = [ldr hdc];
                case 'lh'
                    Sind{1} = ldr+1:ldr+hdr;
                    Sind{2} = 1:ldc;
                    Sdim = [hdr ldc];
                case 'hh'
                    Sind{1} = ldr+1:ldr+hdr;
                    Sind{2} = ldc+1:ldc+hdc;
                    Sdim = [hdr hdc];
            end;
            S = D(Sind{:});
        end

        function Y = idwt_2D(obj,A,H,V,D,wavelet)
            %Two-dimensional separable IDWT
            %Y=idwt_2D(A,H,V,D,wavelet)
            %
            %Input:
            % A,H,V,D - approximation signal, horizontal, vertical and diagonal details
            %           signal
            % wavelet - wavelet identification string, or wavelet data structure
            %
            %Output:
            % Y - reconstructed matrix
            %
            %Uses:
            % load_wavelet.m
            % idwt_dim.m
            %
            %Example:
            % Y = idwt_2D(A,H,V,D,'CDF_9x7');
            
            %load the wavelet here
            if ischar(wavelet)
                wvf = obj.load_wavelet(wavelet,'E');
            else
                wvf = wavelet;
            end
            Y = [A H;V D];
            Y = obj.idwt_dim(Y,2,wvf); %rows
            Y = obj.idwt_dim(Y,1,wvf); %columns
        end
        
        function [D,packet_stream,s,E]=decomp_packets2D(obj,Y,param,entp)
            %2D wavelet packets decomposition with entropy-based subband splitting
            %[D,packet_stream,s,E]=decomp_packets2D(Y,param,entp)
            %
            %Input:
            % Y - array to be transformed
            % param - structure containing decomposition parameters
            %         N: maximum depth of decomposition
            %         pdep: packet decompositon depth
            %         wvf: structure containing properties of a wavelet (see load_wavelet.m)
            %         dec: 'greedy' or otherwise 'full'
            % entp - parameters for splitting criterion based on entropy
            %
            %Output:
            % D - array of wavelet coefficients
            % packet_stream - stream of bits representing information on splitting decisions
            %                 of the wavelet packets decomposition
            % s - structure containing info on parent-children relationship between subbands
            %     given by wavelet packets decomposition
            %     s.scale: decomposition level of the subband
            %     s.parent: subband's parent
            %     s.children: subband's children
            %     s.band_abs: absolute position and size [x y dx dy] (in pixels)
            %     s.ban_rel: position and size relative to the original image size
            % E - entropy of the resulting decomposition
            %
            %Uses:
            % dwt_2D.m
            %
            %Example:
            % par=struct('N',5,'pdep',2,'wvf',load_wavelet('CDF_9x7'),'dec','greedy');
            % ent_par=struct('ent','shannon','opt',0);
            % [D,packet_stream,s,E]=decomp_packets2D('lena256.png',par,ent_par);
            % draw_packets(D,par.N,par.pdep,s,packet_stream); %displays the result of performed decomposition
            
            if ischar(Y)
                Y=imread(Y);
            end
            Y=double(Y);
            [Drows,Dcols]=size(Y);
            subr=Drows/(2^param.N); %dimensions of the lowest subband
            subc=Dcols/(2^param.N);
            D=zeros(Drows,Dcols);
            %check the number of decompositions
            if (round(subr) ~= subr) || (round(subc) ~= subc)
                %at the moment, only powers of two supported
                %error('Illegal number of decompositions for a given matrix!');
                cprintf('err', 'Wavelet - Illegal number of decompositions for a given matrix! \n');
                cprintf('err', 'Matrix has been turned power of for levels \n');
                Y = makematrixdivisible(Y,param.N);
                
            end
            %initialize the packet tree structure
            s=init_packettree(param.N,subr,subc);
            %reserving memory for D (coefficients matrix), and Dpom (helping matrix)
            packet_stream=[];
            E=0;
            %performing the wavelet decomposition of N levels
            N=param.N;
            param.N=param.N-param.pdep;
            for i=1:N
                [DA,DH,DV,DD]=obj.dwt_2D(Y,param.wvf);
                siz=size(DA);
                if (i < N) && (param.pdep > 0) %checking the subbands only if it's not the lowest level of decomposition
                    subbindex=(N-i+1)*3-1; %index of subband H in array structure 's'
                    [DH,Eh,packet_stream,s]=decompose_subband(subbindex,DH,param,entp,packet_stream,s);
                    subbindex=(N-i+1)*3;   %index of subband V in array structure 's'
                    [DV,Ev,packet_stream,s]=decompose_subband(subbindex,DV,param,entp,packet_stream,s);
                    subbindex=(N-i+1)*3+1; %index of subband D in array structure 's'
                    [DD,Ed,packet_stream,s]=decompose_subband(subbindex,DD,param,entp,packet_stream,s);
                else
                    Eh=subb_entropy(DH,entp.ent,entp.opt);
                    Ev=subb_entropy(DV,entp.ent,entp.opt);
                    Ed=subb_entropy(DD,entp.ent,entp.opt);
                end;
                E=E+Eh+Ev+Ed;
                A=[DA DH; DV DD];
                D(1:2*siz(1),1:2*siz(2))=A;
                Y=DA;
                param.pdep=param.pdep-1;
            end;
            E=E+subb_entropy(DA,entp.ent,entp.opt)+subb_entropy(DV,entp.ent,entp.opt)+...
                subb_entropy(DH,entp.ent,entp.opt)+subb_entropy(DH,entp.ent,entp.opt);
            
            function [band,E,p_stream,s]=decompose_subband(subbindex,band,param,entp,p_stream,s)
                %Decompose subband a bit further
                s0=struct('scale',0,'parent',0,'children',0,'band_abs',0,'band_rel',0);
                coord=struct('abs',0,'rel',0);
                parentindex=s(subbindex).parent;
                s0.parent=parentindex;
                coord.abs=s(subbindex).band_abs;
                coord.rel=s(subbindex).band_rel;
                if strcmp(param.dec,'greedy')
                    [band,E,p_stream,sub_list]=dec_greedy(coord,band,param,entp,p_stream,s0);
                else
                    [band,E,p_stream,sub_list]=dec_full(coord,band,param,entp,p_stream,s0);
                end;
                s=resolve_conflicts(subbindex,parentindex,s,sub_list);
            end
            function [band,E,p_stream,sub_list]=dec_full(obj,coord,band,param,entp,p_stream,sub_list)
                %Function for full packet decomposition and best basis selection. Provides optimal performance
                E=subb_entropy(band,entp.ent,entp.opt);
                if param.pdep>0
                    cntold=size(p_stream,2)+1;
                    [banda,bandh,bandv,bandd]=obj.dwt_2D(band,param.wvf);
                    siz=size(banda);
                    A_coord=struct('abs',[coord.abs(1:2) siz(2) siz(1)],'rel',[coord.rel(1:2) coord.rel(3:4)/2]);
                    H_coord.abs=A_coord.abs+[siz(2) 0 0 0];
                    H_coord.rel=A_coord.rel+[A_coord.rel(3) 0 0 0];
                    V_coord.abs=A_coord.abs+[0 siz(1) 0 0];
                    V_coord.rel=A_coord.rel+[0 A_coord.rel(4) 0 0];
                    D_coord.abs=A_coord.abs+[siz(2) siz(1) 0 0];
                    D_coord.rel=A_coord.rel+[A_coord.rel(3) A_coord.rel(4) 0 0];
                    param.pdep=param.pdep-1;
                    [banda,Ea,p_stream,list_A]=dec_full(A_coord,banda,param,entp,p_stream,sub_list);
                    [bandh,Eh,p_stream,list_H]=dec_full(H_coord,bandh,param,entp,p_stream,sub_list);
                    [bandv,Ev,p_stream,list_V]=dec_full(V_coord,bandv,param,entp,p_stream,sub_list);
                    [bandd,Ed,p_stream,list_D]=dec_full(D_coord,bandd,param,entp,p_stream,sub_list);
                    if Ea+Eh+Ev+Ed<E
                        p_stream(size(p_stream,2)+1)=1;
                        band=[banda bandh;bandv bandd];
                        E=Ea+Eh+Ev+Ed;
                        sub_list=[list_A list_H list_V list_D];
                    else
                        p_stream(cntold)=0;
                        p_stream=p_stream(1:cntold);
                        sub_list=add_to_list(sub_list,param.N+param.pdep+1,coord,sub_list(1).parent,0);
                    end;
                else
                    sub_list=add_to_list(sub_list,param.N+param.pdep,coord,sub_list(1).parent,0);
                end;
            end
            function [band,E,p_stream,sub_list]=dec_greedy(obj,coord,band,param,entp,p_stream,sub_list)
                %Function for greedy packet decomposition and best basis selection. Provides faster preformance
                [banda,bandh,bandv,bandd]=obj.dwt_2D(band,param.wvf);
                siz=size(banda);
                E =subb_entropy(band,entp.ent,entp.opt);
                Ea=subb_entropy(banda,entp.ent,entp.opt);
                Eh=subb_entropy(bandh,entp.ent,entp.opt);
                Ev=subb_entropy(bandv,entp.ent,entp.opt);
                Ed=subb_entropy(bandd,entp.ent,entp.opt);
                if Ea+Eh+Ev+Ed<E
                    A_coord=struct('abs',[coord.abs(1:2) siz(2) siz(1)],'rel',[coord.rel(1:2) coord.rel(3:4)/2]);
                    H_coord.abs=A_coord.abs+[siz(2) 0 0 0];
                    H_coord.rel=A_coord.rel+[A_coord.rel(3) 0 0 0];
                    V_coord.abs=A_coord.abs+[0 siz(1) 0 0];
                    V_coord.rel=A_coord.rel+[0 A_coord.rel(4) 0 0];
                    D_coord.abs=A_coord.abs+[siz(2) siz(1) 0 0];
                    D_coord.rel=A_coord.rel+[A_coord.rel(3) A_coord.rel(4) 0 0];
                    param.pdep=param.pdep-1;
                    if param.pdep>0
                        [banda,Ea,p_stream,sub_list]=dec_greedy(A_coord,banda,param,entp,p_stream,sub_list);
                        [bandh,Eh,p_stream,sub_list]=dec_greedy(H_coord,bandh,param,entp,p_stream,sub_list);
                        [bandv,Ev,p_stream,sub_list]=dec_greedy(V_coord,bandv,param,entp,p_stream,sub_list);
                        [bandd,Ed,p_stream,sub_list]=dec_greedy(D_coord,bandd,param,entp,p_stream,sub_list);
                    else
                        sub_list=add_to_list(sub_list,param.N+param.pdep,A_coord,sub_list(1).parent,0);
                        sub_list=add_to_list(sub_list,param.N+param.pdep,H_coord,sub_list(1).parent,0);
                        sub_list=add_to_list(sub_list,param.N+param.pdep,V_coord,sub_list(1).parent,0);
                        sub_list=add_to_list(sub_list,param.N+param.pdep,D_coord,sub_list(1).parent,0);
                    end;
                    p_stream(size(p_stream,2)+1)=1;
                    band=[banda bandh;bandv bandd];
                    E=Ea+Eh+Ev+Ed;
                else
                    p_stream(size(p_stream,2)+1)=0;
                    sub_list=add_to_list(sub_list,param.N+param.pdep,coord,sub_list(1).parent,0);
                end;
            end
            function sub_list=add_to_list(sub_list,scale,coord,parent,children)
                ind=size(sub_list,2);
                if sub_list(ind).scale
                    ind=ind+1;
                end;
                sub_list(ind).scale=scale;
                sub_list(ind).band_abs=coord.abs;
                sub_list(ind).band_rel=coord.rel;
                sub_list(ind).parent=parent;
                sub_list(ind).children=children;
            end
            function s=resolve_conflicts(subbindex,parentindex,s,sub_list)
                siz_new_list=size(sub_list,2);
                if siz_new_list>1 %if subband was further decomposed
                    siz_list=size(s,2);
                    new_children=siz_list+1:siz_list+siz_new_list; %indices of new subbands in 's'
                    s(new_children)=sub_list(1:siz_new_list); %then add new subbands in 's'
                else %else subband H stays as child of it's parent
                    new_children=subbindex;
                end;
                move_children_up=[]; %children that will be moved upwards to low parent (or ancestor) in subband tree structure
                prev_children=s(subbindex).children; %children assigned to subband in previous level of decomposition
                %resolve parenting conflicts between current and previous level of decomposition
                if prev_children(1)>0 %if it's not the subband from the highest level
                    siz_prev=size(prev_children,2);
                    coord_prev=reshape([s(prev_children).band_rel],4,siz_prev)'; %coordinates of subbands
                    siz_new=size(new_children,2);
                    for k=1:siz_new
                        band_scale=s(new_children(k)).scale;
                        int_areas_overlap=rectint(s(new_children(k)).band_rel,coord_prev)>0;
                        overlapping_bands=prev_children(int_areas_overlap);
                        scale_diff=[s(overlapping_bands).scale]-band_scale; %scale diference of the overlapping bands
                        s(new_children(k)).children=overlapping_bands(scale_diff>-1); %link subband to children
                        if isempty(s(new_children(k)).children)
                            s(new_children(k)).children=0;
                        end;
                        if ~isempty(find(scale_diff>1,1))
                            disp(['Heavy parenting conflict resolved on subband ' num2str(overlapping_bands(scale_diff>1))]);
                        end;
                        move_children_up=[move_children_up overlapping_bands(scale_diff<=-1)]; %move up in subband tree
                        indices=~int_areas_overlap; %just the rest will be tested again
                        prev_children=prev_children(indices);
                        coord_prev=coord_prev(indices,:);
                    end;
                end;
                s(parentindex).children(1)=[]; %remove first child (from initial structure)
                s(parentindex).children=[s(parentindex).children new_children move_children_up]; %and assign new children to parent
                if isempty(s(parentindex).children)
                    s(parentindex).children=0;
                end;
            end
            function s=init_packettree(N,subr,subc)
                %subband ordering is H,V,D
                bands=N*3+1;
                s=repmat(struct('scale',0,'parent',0,'children',0,'band_abs',0,'band_rel',0),1,bands);
                s(1).scale=1;
                s(1).parent=0;
                s(1).children=[2 3 4];
                s(1).band_abs=[0 0 subc subr]; %coordinates that define subband as a rectangle
                s(1).band_rel=[0 0 1 1];
                for i=2:bands
                    mul=2^floor((i-2)/3);
                    subrs=subr*mul;
                    subcs=subc*mul;
                    s(i).scale=floor((i+1)/3);
                    s(i).parent=i-3;
                    s(i).children=i+3;
                    s(i).band_abs=[(mod(i,3)>0)*subcs (mod(i+1,3)>0)*subrs subcs subrs];
                    s(i).band_rel=[0 0 1 1];
                end;
                s(2).parent=1;
                s(3).parent=1;
                s(4).parent=1;
                s(bands-2).children=0;
                s(bands-1).children=0;
                s(bands).children=0;
            end
            function ent = subb_entropy(A,type,par)
                switch type
                    case 'shannon'
                        A = A(find(A)).^2;
                        ent = -sum(A.*log(A));
                    case 'threshold'     % par is the threshold.
                        ent = sum(abs(A) > par);
                    case 'logenergy'     % in3 not used.
                        A = A(find(A)).^2;
                        ent = sum(log(A));
                end;
            end
            function  Y = makematrixdivisible(Y,level)
                %Make image divisible by 2 in the levels
                if mod(size(Y,1),2) ~= 0
                    Y(end+1,:)=Y(end,:);
                end
                if mod(size(Y,2),2) ~= 0
                    Y(:,end+1)=Y(:,end);
                end
                extension = 1;
                while mod(size(Y,1),2^level) ~= 0
                    Y = wextend('addrow','sym',Y,extension);
                end
                while mod(size(Y,2),2^level) ~= 0
                    Y = wextend('addcol','sym',Y,extension);
                end
            end
        end
        
        function A=recon_packets2D(obj,D,wavelet_info)
            %2D wavelet packets reconstruction
            %A=recon_packets2D(D,param,packet_stream)
            %
            %Input:
            % D - array of wavelet coefficients
            % param - structure containing decomposition parameters (see in
            %         decomp_packets.m)
            % packet_stream - stream of bits representing information on splitting decisions
            %                 of wavelet packets decomposition
            %
            %Output:
            % A - reconstructed array
            %
            %Uses:
            % idwt_2D.m
            %
            %Example:
            % [D,packet_stream]=decomp_packets2D(Y,par,ent_par);%see decomp_packets2D.m
            % A=recon_packets2D(D,par,packet_stream);
            
            param = wavelet_info.par;
            packet_stream = wavelet_info.packet_stream;
            
            param.pdep=param.pdep-param.N+1; %e.g. if N=5 and pdep=2 -> param.pdep=-2
            if size(packet_stream,2)>1
                packet_stream=fliplr(packet_stream);
                entropy=1;
            else
                entropy=0;
            end
            cnt=0;
            siz=size(D)/(2^param.N);
            for i=1:param.N
                D1=D(1:siz(1),1:siz(2));
                DH1=D(1:siz(1),siz(2)+1:2*siz(2));
                DV1=D(siz(1)+1:2*siz(1),1:siz(2));
                DD1=D(siz(1)+1:2*siz(1),siz(2)+1:2*siz(2));
                if (i > 1) && (param.pdep > 0) && entropy
                    [DD1,cnt]=recon_entropy(DD1,param,packet_stream,cnt);
                    [DV1,cnt]=recon_entropy(DV1,param,packet_stream,cnt);
                    [DH1,cnt]=recon_entropy(DH1,param,packet_stream,cnt);
                end
                A=obj.idwt_2D(D1,DH1,DV1,DD1,param.wvf);
                siz=siz*2;
                D(1:siz(1),1:siz(2))=A;
                param.pdep=param.pdep+1;
            end
            
            function [band,cnt]=recon_entropy(band,param,packet_stream,cnt)
                cnt=cnt+1;
                if packet_stream(cnt)
                    param.pdep=param.pdep-1;
                    siz=size(band)/2;
                    D1=band(1:siz(1),1:siz(2));
                    DH1=band(1:siz(1),siz(2)+1:2*siz(2));
                    DV1=band(siz(1)+1:2*siz(1),1:siz(2));
                    DD1=band(siz(1)+1:2*siz(1),siz(2)+1:2*siz(2));
                    if param.pdep>0
                        [DD1,cnt]=recon_entropy(DD1,param,packet_stream,cnt);
                        [DV1,cnt]=recon_entropy(DV1,param,packet_stream,cnt);
                        [DH1,cnt]=recon_entropy(DH1,param,packet_stream,cnt);
                        [D1,cnt]=recon_entropy(D1,param,packet_stream,cnt);
                    end
                    band=idwt_2D(D1,DH1,DV1,DD1,param.wvf);
                end
            end
        end
        
    end
    methods (Access = private, Static)
        
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
            % %multiple-line comment 
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
            %  coefficients of the low-pass analysis filter, one per line - delay,value  
            % %comment(filt_H1)
            %  coefficients of the high-pass analysis filter, one per line - delay,value  
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
            end
            
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
            end
        end
        
        function [S,Sind] = submatrix(X,Ssiz,Soff)
            %Extracts submatrix from a multidimensional matrix
            %[S,Sind]=submatrix(X,Ssiz,Soff)
            %
            %Input:
            % X - matrix of n dimensions
            % Ssiz - submatrix size to be extracted
            % Soff - submatrix starts at this offset
            %
            %Output:
            % S - submatrix
            % Sind - indices of the extracted submatrix in the original submatrix
            %
            %Note:
            % To get the extracted indices use X(Sind{:})
            %
            %Example:
            % Suppose I have an N-dimensional matrix X, and want to extract the
            % submatrix of which each dimension is half than in X (with
            % the offset at the first element of X). Then, the function would be called
            % with: S = submatrix(X,ceil(size(X)/2));
            
            if isvector(X)
                n = 1;
            else
                n = ndims(X);
            end;
            if (nargin < 3)
                Soff = ones(n,1); %sets the default offset
            end;
            Sind = repmat({':'},n,1);
            for i = 1:n
                Sind{i} = Soff(i):(Soff(i) + Ssiz(i) - 1);
            end;
            S = X(Sind{:});
        end
        
        function [ld,hd,N] = subband_dim(sdim, N)
            %Computes the subband dimensions for a specified number of decompositions
            %[ld,hd,N]=subband_dim(sdim, N)
            %
            %Input:
            % sdim - vector of lengths of the original signal
            % N - number of signal decompositions
            %
            %Output:
            % ld - size of the low-pass signal after N levels of decomposition
            % hd - size of the high-pass signal after N levels of decomposition
            % N - number of actually performed number of decompositions
            %
            %Note:
            % If N is too large, the function will return ld = 1, and N will equal
            % the number of allowed decompositions.
            % In the smallest non-singleton dimension direction the low-pass subband
            % can have only one coefficient.
            %
            %Example:
            % [ld,hd,N]= subband_dim(9, 3); %if upper_limit = 1 -> ld = 2, hd = 1
            % [ld,hd,N]= subband_dim(100, Inf); %for max. number of decompositions
            
            sdim = double(sdim);
            ld = zeros(size(sdim));
            hd = zeros(size(sdim));
            n = zeros(size(sdim));
            
            for d=1:length(sdim)
                if (sdim(d) == 1) %singleton dimension, skip it!
                    ld(d) = 1;
                    hd(d) = 1;
                    n(d) = N;
                else
                    n(d) = sb_dim(sdim(d),N);
                end;
            end;
            
            N = min(n);
            for d=1:length(sdim)
                [n(d),ld(d),hd(d)] = sb_dim(sdim(d),N);
            end;
            
            function [i,ld,hd] = sb_dim(sdim,N)
                ld = sdim;
                hd = 0;
                for i=1:N
                    ldold = ld;
                    ld = ceil(ldold / 2);
                    hd = ldold - ld;
                    if (ld == 1)
                        break;
                    end;
                end
            end
        end
        
        function [y] = makepowerof2(x)
            N = length(x);
            y = x;
            while mod(log(N)/log(2),1)~=0
                y(N+1) = 0;
                N = N+1;
            end
        end
        
        function draw_packets(D,N,pdep,s,packet_stream)
            %Visualises the wavelet packets decomposition
            %draw_packets(D,N,pdep,s,packet_stream)
            %
            %Input:
            % D - array of wavelet coefficients
            % N - number of dyadic (pre-packet) decompositions
            % pdep - "packet decomposition depth"
            % s - structure containing info on parent-children relationship between subbands
            %     (see in decomp_packets.m)
            % packet_stream - stream of bits representing information on splitting decisions
            %                 of wavelet packets decomposition
            %
            %Note:
            % Draws 4 plots.
            %
            %Example:
            % par=struct('N',5,'pdep',2,'wvf',load_wavelet('CDF_9x7'),'dec','greedy');
            % ent_par=struct('ent','shannon','opt',0);
            % [D,packet_stream,s,E]=decomp_packets2D('lena256.png',par,ent_par);
            % draw_packets(D,par.N,par.pdep,s,packet_stream);
            
            scrsz = get(0,'ScreenSize');
            figure('Name','Wavelet packet structure drawn by bit information');
            pdep=pdep-N+1;
            if size(packet_stream,2)>1
                packet_stream=fliplr(packet_stream);
            end;
            cnt=0;
            set(gca,'YDir','reverse','PlotBoxAspectRatio',[1 1 1],'XTick',[],'YTick',[]);
            siz=1/(2^N);
            for i=1:N
                rectangle('Position',[0,0,siz,siz]);
                rectangle('Position',[0,siz,siz,siz]);
                rectangle('Position',[siz,0,siz,siz]);
                rectangle('Position',[siz,siz,siz,siz]);
                if pdep>0
                    cnt=draw_subband(siz,siz,packet_stream,cnt,siz,pdep);
                    cnt=draw_subband(0,siz,packet_stream,cnt,siz,pdep);
                    cnt=draw_subband(siz,0,packet_stream,cnt,siz,pdep);
                end;
                siz=siz*2;
                pdep=pdep+1;
            end;
            %second plot - drawn by subband structure s (small - convinient for copy/paste)
            figure('Position',[300 200 200 200],'Name','Wavelet packet structure drawn by ''s'' subband structure');
            set(gca,'YDir','reverse','PlotBoxAspectRatio',[1 1 1],'XTick',[],'YTick',[],'Position',[0.005 0.005 0.99 0.99]);
            axis([0 size(D,2) 0 size(D,1)]);
            wavelet_dec_bands=N*3+1;
            for i=1:size(s,2)
                if i<=wavelet_dec_bands
                    lw=2;
                else
                    lw=1;
                end;
                rectangle('Position',s(i).band_abs,'LineWidth',lw);
            end;
            %third plot - drawn by subband structure s, with linked subbands
            %figure('Position',[1 scrsz(4)/2 scrsz(3)/2 scrsz(4)/2])
            figure('Position',[scrsz(3)/2-3*scrsz(4)/8 scrsz(4)/2-3*scrsz(4)/8 6*scrsz(4)/8 6*scrsz(4)/8],...
                'Name','Wavelet packet structure drawn by s - subabnd structure');
            set(gca,'YDir','reverse','PlotBoxAspectRatio',[1 1 1],'XTick',[],'YTick',[],'Position',[0.005 0.005 0.99 0.99]);
            axis([0 size(D,2) 0 size(D,1)]);
            wavelet_dec_bands=N*3+1;
            for i=1:size(s,2)
                if i<=wavelet_dec_bands
                    lw=2;
                else
                    lw=1;
                end;
                rectangle('Position',s(i).band_abs,'LineWidth',lw);
            end;
            cmap=colormap('prism');
            lv=1;
            link_subbands(s(1),s,cmap,lv);
            %and finally draw decomposed subbands
            figure('Name','Wavelet packet transform coefficients');
            set(gca,'YDir','reverse','PlotBoxAspectRatio',[1 1 1],'XTick',[],'YTick',[]);
            axis([0 size(D,2) 0 size(D,1)]);
            image('CData',100*log10(abs(D)));colormap(gray(256));
            for i=1:size(s,2)
                if i<=wavelet_dec_bands
                    lw=2;
                else
                    lw=1;
                end;
                rectangle('Position',s(i).band_abs,'LineWidth',lw,'EdgeColor','white');
            end;
            
            function cnt=draw_subband(sx,sy,packet_stream,cnt,siz,pdep)
                cnt=cnt+1;
                if packet_stream(cnt)
                    pdep=pdep-1;
                    siz=siz/2;
                    rectangle('Position',[sx,sy,siz,siz]);
                    rectangle('Position',[sx+siz,sy,siz,siz]);
                    rectangle('Position',[sx,sy+siz,siz,siz]);
                    rectangle('Position',[sx+siz,sy+siz,siz,siz]);
                    if pdep>0
                        cnt=draw_subband(sx+siz,sy+siz,packet_stream,cnt,siz,pdep);
                        cnt=draw_subband(sx,sy+siz,packet_stream,cnt,siz,pdep);
                        cnt=draw_subband(sx+siz,sy,packet_stream,cnt,siz,pdep);
                        cnt=draw_subband(sx,sy,packet_stream,cnt,siz,pdep);
                    end;
                    %siz=siz*2;
                end;
            end
            function link_subbands(node,s,cmap,lv)
                n=node.children;
                if n(1)>0
                    for i=1:size(n,2)
                        child=s(n(i));
                        link_subbands(child,s,cmap,lv+1);
                        %krc0=floor(lv/2)/lv;
                        %krc1=floor((lv+1)/2)/(lv+1);
                        childcoordx=child.band_abs(1)+child.band_abs(3)/4;
                        childcoordy=child.band_abs(2)+child.band_abs(4)/4;
                        parentcoordx=node.band_abs(1)+node.band_abs(3)/2;
                        parentcoordy=node.band_abs(2)+node.band_abs(4)/2;
                        line([parentcoordx childcoordx],[parentcoordy childcoordy],'Color',cmap(lv,:));
                        text(childcoordx,childcoordy,int2str(n(i)));
                        %next line for use with arrow.m by F. Golnaraghi et al.
                        %arrow([parentcoordx parentcoordy],[childcoordx childcoordy],'Length',10);
                    end;
                end;
            end
        end
        
        function plotSpectogram(wave)
            figure
            imagesc(abs(wave))
            xlabel('Time (integer index)')
            ylabel('Scale')
        end
    end
end