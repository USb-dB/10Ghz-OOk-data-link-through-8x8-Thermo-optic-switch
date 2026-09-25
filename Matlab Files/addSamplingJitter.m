function y = addSamplingJitter(x, sigma_UI, Ns)
    n = (1:length(x));
    jitter = sigma_UI * Ns * randn(size(n));
    y = interp1(n, x, n + jitter, 'linear', 'extrap');
end
