function example_dwt
folder = fileparts(which(mfilename)); 
addpath(genpath(folder));
d.type = 'TEMPORAL';
d.typesignal = 'AUDIO';
d.name = 'sinus';
d.motherwave = 'Haar';
d.type_ft = 'matlab';
d.type_wt = 'multilevel';
a = Data(d); 

[c,l] = wavedec(a.signal,1,'Haar');
%[cA,cD] = dwt(a.signal,'Haar');

%Compare results


end


%%  Errores o cosas a mejorar
%1. Cuando se hace la stft tambien quiere hacer despues una inversa despues y no deberia
%2. guardar las reconstrucciones porque no se como hacer para que no se sobre escriban las reconstrucciones ft y la wave
% Esto tengo que crear un nuevo struct (clase) para cada uno de los metodos de
% compresion donde se guarde su recuperacion. No guardarlos todos en la
% misma clase. (Pasar de guardar en rec_f a signal)