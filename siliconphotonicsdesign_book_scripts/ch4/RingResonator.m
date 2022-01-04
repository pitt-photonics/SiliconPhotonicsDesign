%% Ring resonator spectrum
clear all

lambda = linspace(1.5, 1.6, 5000);    % wavelength in um
R = 3.1;                              % radius in um
k = 0.31;
Filter_type = 'add-drop';
[Ethru, Edrop] = RingSpectrum(lambda, R, k, Filter_type);

subplot(211)
plot(lambda, abs(Ethru).^2, lambda, abs(Edrop).^2)
xlim([min(lambda),max(lambda)])
ylim([0, 1])
xlabel('Wavelength (um)')
ylabel('Transmission')
legend({'Through Port','Drop Port'})

subplot(212)
semilogy(lambda, abs(Ethru).^2, lambda, abs(Edrop).^2)
xlim([min(lambda),max(lambda)])
ylim([1e-3, 1])
xlabel('Wavelength (um)')
ylabel('Transmission')

function [Ethru, Edrop] = RingSpectrum(lambda, R, k, Filter_type)
    % ring parameters
    neff = 1.91;                           % neff for SOI TE mode
    ng = 4.63;                              % group index
    L_rt = 2*pi*R;                          % round trip length
    phi_rt = (2*pi*neff./lambda)*L_rt;      % round trip phase
    alpha_wg_dB=5;                          % optical loss of intrinsic optical waveguide, in dB/cm
    alpha_wg=-log(10^(-alpha_wg_dB/10));    % converted to /cm
    A = exp(-alpha_wg*L_rt/1e4);            % optical loss from one round trip

    % coupling coefficients
    if (Filter_type=='all-pass')
        t=sqrt(1-k^2);
        Ethru=(-sqrt(A)+t*exp(-1i*phi_rt))./(-sqrt(A)*conj(t)+exp(-1i*phi_rt));
        Edrop=zeros(size(lambda));
        Qc=-(pi*L_rt*ng)/(mean(lambda)*log(abs(t)));
    elseif (Filter_type=='add-drop')
        k1=k; k2=k1;
        t1=sqrt(1-k1^2);  t2=sqrt(1-k2^2);
        Ethru=(t1-conj(t2)*sqrt(A)*exp(1i*phi_rt))./(1-sqrt(A)*conj(t1)*conj(t2)*exp(1i*phi_rt));
        Edrop=-conj(k1)*k2*sqrt(sqrt(A))*exp(1i*phi_rt/2)./(1-sqrt(A)*conj(t1)*conj(t2)*exp(1i*phi_rt));
        Qc1=-(pi*L_rt*ng)/(mean(lambda)*log(abs(t1)));
        Qc2=-(pi*L_rt*ng)/(mean(lambda)*log(abs(t2)));
        Qc=1/(1/Qc1+1/Qc2);
    else
        error(1, 'The''Filter_type'' has to be ''all-pass'' or ''add-drop''.\n');
    end
end