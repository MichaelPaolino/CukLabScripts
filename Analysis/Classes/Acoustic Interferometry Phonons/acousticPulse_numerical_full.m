%% Set initial parameters (currently set to unitless)
%Useful constants to convert to real units
R = 0.928; %acoustic reflection coefficient between STO and H2O
B = 170.57; %Bulk modulus in GPa
G = 108.46; %Shear modulus in GPa
rho = 5.11; %density in g/cm2
c_optical = 299792.458; %speed of light in nm/ps
c_acoustic = sqrt((B+4/3*G)*1e10/rho)*1e-5; %acoustic velocity in nm/ps calcualted from elasticity constants, should be close to 7.9 nm/ps

lambda = 15; %penetration depth in nm
dt = 0.2; %in ps
dx = 0.5; %in nm
t = 0:dt:30; %in ps
x = 0:dx:(max(t+0.8)*c_acoustic); %in nm
t0 = 0.1; %pulse duration in ps

%tau = ([9.9 10 10.1]*7.85/15)+4; %for directly inputting tau
%xi = 0:0.1:(1.2*max(tau));

tau = t*c_acoustic/lambda; %dimensionless time
xi = x/lambda;  %dimensionless space

tau0 = t0*c_acoustic/lambda;%t0/2.355*c/lambda; %dimensionless pulse duration
offset = 0; theta = @(xipt,tau,tau0) exp(-xipt).*(1-exp(-tau/tau0));
theta_matrix = @(xipt,tau,tau0) exp(-xipt')*(1-exp(-tau/tau0));
% offset = 4; offset_tmp = erf(offset/sqrt(2)); 
% theta = @(xipt,tau,tau0) 0.5*exp(-xipt).*(offset_tmp+erf((tau/tau0-offset)/sqrt(2)));
% theta_matrix = @(xipt,tau,tau0) 0.5*exp(-xipt')*(offset_tmp+erf((tau/tau0-offset)/sqrt(2)));
tau = tau + offset*tau0; %in ps

% alpha_therm = 3.23e-5; %1/K
% Cp = 2.730; %J/cm^3K
% R_pump = 0.8932*0.9989*0.9648;
% Q = R_pump*0.6*1e-3*(1-3.25/4.64360); %J/cm^2
% eta0 = alpha_therm/3*Q/(lambda*1e-7*Cp);

eta0 = -0.06;
U_conv = eta0*3*lambda*B*1e10/(rho*(c_acoustic*1e5)^2); %converts between unitless and unit displacement in nm
eta_conv = eta0*3*B*1e10/(rho*(c_acoustic*1e5)^2); %converts between unitless and unit strain
sigma_conv = eta_conv*(B+4*G/3);

%% This generates a 2D grid for theta(xi,tau), with xi as rows and tau as cols
%This runs really slow because it needs to calculate multiple values of t
%and x, run the next one below for just the acoustic pulse 
% theta = exp(-xi)'*g_int;

% Calcualte V as a function of x and t
V = zeros(length(xi),length(tau));

for ii = 1:length(tau)
    tic;
    for jj = 1:length(xi)
        if tau(ii)<=xi(jj)
            int_1 = integral(@(y) theta(y,(tau(ii)-xi(jj)+y),tau0), xi(jj)-tau(ii), xi(jj));
            int_2 = integral(@(y) theta(y,(tau(ii)+xi(jj)-y),tau0), xi(jj), xi(jj)+tau(ii));
            V(jj,ii) = 1/2*(int_1-int_2);
        else
            int_1 = integral(@(y) theta(y,(tau(ii)-xi(jj)+y),tau0), 0, xi(jj));
            int_2 = integral(@(y) theta(y,(tau(ii)+xi(jj)-y),tau0), xi(jj), xi(jj)+tau(ii));
            int_3 = integral(@(y) theta(y,(tau(ii)-xi(jj)-y),tau0), 0, tau(ii)-xi(jj)); 
            V(jj,ii) = 1/2*(int_1-int_2-R*int_3);
        end
    end
    toc
end

%%
% figure; plot(t,g,t,g_int);
% figure; contour(x,t,theta');
video = false;
real_units = true;
if video
    writerObj = VideoWriter('CC_disp_strain_scale.avi','Motion JPEG AVI');
    writerObj.Quality = 80;
    writerObj.FrameRate = 15;
    open(writerObj);
end

Vtrue = V;
Strain = diff(Vtrue,1,1)/(xi(2)-xi(1));
Stress = diff(Vtrue,1,1)/(xi(2)-xi(1))-theta_matrix(xi(2:end)-0.5*(xi(2)-xi(1)),tau,tau0);
Velocity = diff(Vtrue,1,2)/(tau(2)-tau(1)); Velocity = [zeros(size(Velocity,1),1) interp1_2D(tau(1:end-1)+0.5*(tau(2)-tau(1)),Velocity,tau(2:end-1))];
Acceleration = diff(Vtrue,2,2)/(tau(2)-tau(1))^2; Acceleration = [zeros(size(Acceleration,1),1) Acceleration];
xip = xi(2:end)-0.5*(xi(2)-xi(1));

if real_units
   xi_tmp = xi*lambda;
   xip_tmp = xip*lambda;
   
   Vtrue = 1000*Vtrue*U_conv;
   Strain = Strain*eta_conv;
   Stress = Stress*sigma_conv;
   Velocity = 1000*Velocity*U_conv*lambda/c_acoustic;
   Acceleration = 1000*Acceleration*U_conv*(lambda/c_acoustic)^2;
   
   xlabval = 'Depth (nm)';
else
   xi_tmp = xi;
   xip_tmp = xip;
   xlabval = '\xi';
end

figure('Position',[230 250 590 325]);
for ii = 1:1:(length(tau)-2)
    subplot(2,1,1);
        plot([0, max(xi_tmp)],[0,0],'--k');    
        hold on; plot(xi_tmp,Vtrue(:,ii)); hold off;
        %set(gca, 'xdir', 'reverse');
        if real_units
            title(['\lambda = ' num2str(lambda) ' nm, t_0 = ' num2str(tau0*lambda/c_acoustic) ' ps, t = ' num2str((tau(ii)-offset*tau0)*lambda/c_acoustic,'%.1f') ' ps']); %-offset*tau0 for OC
        else
            title(['\tau_0 = ' num2str(tau0) ', \tau = ' num2str(tau(ii),'%.1f')]);
        end
        ylim(1.1*[min(min(Vtrue)), max(max(Vtrue))]);
        ylabel('Disp. (pm)');
        xlim([0 max(xi_tmp)])
%     subplot(5,1,2);
%         plot([0, max(xi_tmp)],[0,0],'--k'); 
%         hold on; plot(xi_tmp,Velocity(:,ii)); hold off;
%         %set(gca, 'xdir', 'reverse');
%         ylim(1.1*[min(min(Velocity)), max(max(Velocity))]);
%         ylabel('Vel. (pm/ps)');
%         xlim([0 max(xi_tmp)]);
%     subplot(5,1,3);
%         plot([0, max(xi_tmp)],[0,0],'--k');
%         hold on; plot(xi_tmp,Acceleration(:,ii)); hold off;
%         %set(gca, 'xdir', 'reverse');
%         ylim(1.1*[min(min(Acceleration)), max(max(Acceleration))]);
%         ylabel('Accel. (pm/ps^2)');
%         xlim([0 max(xi_tmp)]);
    subplot(2,1,2);
        plot([0, max(xi_tmp)],[0,0],'--k');
        hold on; plot(xip_tmp,Strain(:,ii)); hold off;
        %set(gca, 'xdir', 'reverse');
        ylim(1.1*[min(min(Strain)), max(max(Strain))]);
        %ylim([-0.006, 0.012]);
        ylabel('Strain (%)');
        xlim([0 max(xi_tmp)]);
        xlabel(xlabval);
        
        %inset:
        axes('Position',[.2 .135 .705 .22])
            box on;
            plot([0, max(xi_tmp)],[0,0],'--k');
            hold on;
            plot(xip_tmp,1000*Strain(:,ii));
            hold off;
            xlim([9 100]);
            ylim([-1, 3]);
            %ylabel('% x10^3');
            text(11, 2.5, '%, x10^{-3}');
            set(gca,'xtick',[])
%     subplot(5,1,5);
%         plot([0, max(xi_tmp)],[0,0],'--k');
%         hold on; plot(xip_tmp,Stress(:,ii)); hold off;
%         %set(gca, 'xdir', 'reverse');
%         ylim(1.1*[min(min(Stress)), max(max(Stress))]);
%         ylabel('Stress (GPa)');
%         xlim([0 max(xi_tmp)]);
%         xlabel(xlabval);
    pause(0.01);
    
    if video
        frame = getframe(gcf);
        writeVideo(writerObj,frame);
    end
end

if video
    close(writerObj);
end

%% Calculate scattered signal amplitude
l = 375;
n = 2.4;

kz = 2*pi*n./l;
k0 = 2*pi./l;
z = xip*lambda;
dz = z(2)-z(1);
rdeps = 1i*k0^2/2/kz*sum(Strain.*exp(2*1i*kz*z(:)),1)*dz;

tPlot = (tau-offset*tau0)*lambda/c_acoustic;
dR = abs(rdeps).*cos(angle(rdeps));
expGrowth = 0.5*(max(dR)+min(dR(floor(numel(dR)/2):end)))*(1-exp(-tPlot/t0)).*(1-exp(-tPlot/lambda*c_acoustic));

figure;
plot(tPlot,dR,tPlot,expGrowth);
xlabel('t (ps)');
ylabel('dR');
legend('Numerical','Guess');

figure;
a0 = 0.9; %norm factor for guess plot
phi - 180; %phase for guess plot
plot(tPlot,dR-expGrowth,tPlot,a0*expGrowth.*cos(4*pi*c_acoustic*n/l*tPlot+phi/180*pi));
xlabel('t (ps)');
ylabel('dR');
legend('Numerical','Guess');

%% Plot just one of the above
figure; plot(xip/0.5,Strain(:,75));