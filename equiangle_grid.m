function [thcomp, phcomp] = equiangle_grid(bw)
%bw = bandwidth for sampling
%returns equiangular grid 2bw by 2bw

    theta = linspace(pi/(4*bw), pi*(4*bw - 1)/(4*bw), 2*bw);
    phi = linspace(0, pi*(2*bw - 1)/bw, 2*bw);
    [th, ph] = meshgrid(theta, phi);
    thcomp = th(:);
    phcomp = ph(:);
end

