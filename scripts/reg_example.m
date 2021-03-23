h = waitbar(0,'please wait...','Name','xy registration...');
for i = 1:dim3
    tmp = double(Fim_raw(:,:,i));
    [output, Greg] = dftregistration(fft2(reg_temp),fft2(tmp));
    Fim_reg(:,:,i) = abs(ifft2(Greg));
    waitbar(i/dim3)
end