function [rs,ts,rp,tp,theta_t] = fresnel(n1,n2,theta_i)
% FRESNEL calculates the complex reflection and refraction coefficients
% at an interface between two absorbing media. All inputs can be complex
% valued.
%
% [rs,ts,rp,tp,theta_t] = fresnel(n1,n2,theta_i)
%   Given complex index of refraction vectors n1, n2 [np,1] and complex 
%   angle of incidince scalar or vector theta_i in radians [1,1] or [np,1] 
%   calculates complex fresnel factors:
%       rs: field reflection amplitude for s-pol
%       ts: field transmission amplitude for s-pol
%       rp: field reflection amplitude for p-pol
%       tp: field transmission amplitude for p-pol
%       theta_t: complex transmission angle in radians
%   The angle of incidence theta_i is defined in medium 1 (n1) as the angle
%   between the interface normal and the input field's k-vector. The s- and
%   p-polarizations are defined with respect to the incidence plane, the
%   plane formed between the incident and reflected k-vectors as well as
%   the interface normal vector. s- is E-field perpendicular to the plane,
%   p- is E-field parallel to the plane.

% Snell's law
theta_t = asin(n1./n2.*sin(theta_i));

% Fresnel amplitude factors
rs = (n1.*cos(theta_i)-n2.*cos(theta_t))./(n1.*cos(theta_i)+n2.*cos(theta_t));
ts = 2*n1.*cos(theta_i)./(n1.*cos(theta_i)+n2.*cos(theta_t));

rp = (n2.*cos(theta_i)-n1.*cos(theta_t))./(n2.*cos(theta_i)+n1.*cos(theta_t));
tp = 2*n1.*cos(theta_i)./(n2.*cos(theta_i)+n1.*cos(theta_t));

