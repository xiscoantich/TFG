function trans = Transformer(data,cParams)
switch cParams.transtype
    case 'Fourier'
        ft = FourierTransformer(cParams);
        trans = ft.directTransform(data);
        
    case 'Wavelet'
        wt = WaveletTransformer(cParams);
        trans = wt.directTransform(data);
    case 'PCA'
        cParams.transmethod = 'nothing';
        pca = PCATransformer(cParams);
        trans = pca.directTransform(data);
end
end