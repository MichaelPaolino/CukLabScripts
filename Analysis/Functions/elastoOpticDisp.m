function [delta_epsilon, p_norm, p] = elastoOpticDisp(l, n, l0, D, K, Ed, p_norm_val, l_norm_val)
%Calculate the elasto-optic tensor element as a function of wavelength.
%Elasto-optic tensor relates a change in strain to a change in inverse
%permittivity, also known as impermeability. Code based on the paper:
%Wemple, S.H.; DiDomenico, M., Jr., Theory of the ELasto-Optic Effect in
%Nonmetallic Crystals. Phys. Rev. B. 1970, 1(1), 193.
%
%Inputs:
%   l: wavelength array in nm (1 x n)
%   n: complex index of refraction of elasto-optic material (1 x n)
%   l0: Effective Sellmiere bandgap of material in nm
%   D: Phenomenological deformation potential in eV
%   K: Elasto-optic dispersion parameter
%   Ed: "Dispersion energy" in eV, e.g. 23.7 eV for STO*
%   p_norm_val: experimental value for elasto-optic effect if available
%   l_norm_val: corresponding value at which p_norm_val was measured in nm
%
%Outputs:
%   delta_epsilon: change in permittivity for a unit change in strain
%   p_norm: normalized elasto-optic tensor element (1 x n)
%   p: unnormalized elasto-optic tensor element (1 x n)
%
%Todo:
%   Include usage of complex index of refraction to calculate kappa
%
%*For Ed see: Wemple, S.H.; DiDomenico, M., Jr., Optical Dispersion and the
% Structure of Solids. Phys. Rev. Lett. 1969, 23(20), 1156.
%**For piezo-optic constants for STO see: Levin, S.B.; Field, N. J.; Plock,
% F.M.; Merker, L.; Some Optical Properties of Strontium Titanate Crystal. J
% Opt. Soc. Am. 1955, 45(9), 737.

p = (1-1./real(n.^2)).^2.*(2*D/Ed).*(1+K*(1-l0^2./l.^2));
[~,l_norm_ind] = min(abs(l-l_norm_val));
%this normalizes fit results to experimental measured value, e.g. at 540 nm**
p_norm = p/(p(l_norm_ind)/p_norm_val);    
delta_epsilon = -p_norm.*(n.^4);