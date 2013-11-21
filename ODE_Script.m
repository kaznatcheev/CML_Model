clear all
close all

% non-cancer initial conditions

N_1h=0.25;
N_0h=0.25;
N_0l=0.25;
N_1l=0.25;
P=0;

% cancer initial conditions

Nc_1h=0.375; %all 0.375
Nc_0h=0.375;
Nc_0l=0.375;
Nc_1l=0.375;
Pc=0;

N_1=N_1h+N_1l+Nc_1h+Nc_1l;
N_0=N_0h+N_0l+Nc_0h+Nc_0l;

% parameters 

k1=0.0433;
k0=0.0918;
c=1/16; %upper bound of 1/48? tc-c1 = 16
cp=1/24; %\tilde{tau_c}=24
r=0.005;
d=1/240; %rdeg = 1/2.8;
dp=1/28; %\lambda_m=8 \lambda_p=20

amin = 0.002;
al = 0.25;
ah = 0.75;
amax = 1;

g=al/ah;
b=1/g;

% non-cancer transition 
% 1 -> 0 
nu1=1.9972;
nu2=0.0028;
nu3=8.7716;

% 0 -> 1
n1=1.6667;
n2=0.333;
n3=3.2189;

mu=ah/amax/(nu1+nu2*exp(nu3*N_0));

om=amin/ah/(n1+n2*exp(n3*N_1));

% cancer transition 
% 1 -> 0
nuc1=0.9992;
nuc2=0.00077068;
nuc3=9.9558;

% 0 -> 1
nc1=16.6667;
nc2=8.3333;
nc3=1.3868;
nc4=0.96;

muc=ah/amax/(nuc1+nuc2*exp(nuc3*N_0));

omc=amin/ah*(1/(nc1+nc2*exp(nc3*N_1))+nc4);

%

% treatment transition
% 0 -> 1
ni1=1666.7;
ni2=833.3333;
ni3=1.3863;
ni4=0.0496;

omi=amin/ah*(1/(ni1+ni2*exp(ni3*N_1))+ni4);

% treatment parameters



interesting_plot=0;
plot_step=5000;
dt=0.001;

for tidx=1:20000000 %1000000
    
    N_1h_new=N_1h + N_1h*(c-k1-d-mu)*dt+ N_0h*om*dt;
    
    N_0h_new=N_0h + N_1h*mu*dt - N_0h*om*dt + N_0l*k0*dt;
    
    N_0l_new=N_0l + N_1h*(g*mu)*dt - N_0l*b*om*dt - N_0l*k0*dt;
    
    N_1l_new=N_1l + N_1l*(c - d - r - g*mu)*dt + N_1h*k1*dt + b*om*N_0l*dt;
    
    P_new=P + r*N_1l*dt - dp*P*dt + cp*P*(1-P)*dt;
    
    N_1h = N_1h_new;
    N_0h = N_0h_new;
    N_0l = N_0l_new;
    N_1l = N_1l_new;
    P=min(P_new,10);  % REMOVE THIS????????????
    
    
    Nc_1h_new=Nc_1h + Nc_1h*(c-k1-d-muc)*dt+ Nc_0h*omc*dt;
    
    Nc_0h_new=Nc_0h + Nc_1h*muc*dt - Nc_0h*omc*dt + Nc_0l*k0*dt;
    
    Nc_0l_new=Nc_0l + Nc_1h*(g*muc)*dt - Nc_0l*b*omc*dt - Nc_0l*k0*dt;
    
    Nc_1l_new=Nc_1l + Nc_1l*(c - d - r - g*muc)*dt + Nc_1h*k1*dt + b*omc*Nc_0l*dt;
    
    Pc_new=Pc + r*Nc_1l*dt - dp*Pc*dt + cp*Pc*(1-Pc)*dt;
    
    Nc_1h = Nc_1h_new;
    Nc_0h = Nc_0h_new;
    Nc_0l = Nc_0l_new;
    Nc_1l = Nc_1l_new;
    Pc=min(Pc_new,10);
    
    
    N_1=N_1h+N_1l+Nc_1h+Nc_1l;
    N_0=N_0h+N_0l+Nc_0h+Nc_0l;
    
    mu=ah/amax/(nu1+nu2*exp(nu3*N_0));
    om=amin/ah/(n1+n2*exp(n3*N_1));
    
    muc=ah/amax/(nuc1+nuc2*exp(nuc3*N_0));
    omc=amin/ah*(1/(nc1+nc2*exp(nc3*N_1))+nc4);
    
    if mod(tidx,plot_step)==0
        
        if interesting_plot==0
                        plot(tidx/plot_step,N_1h,'.b','markersize',5)
                        hold on
                        plot(tidx/plot_step,N_0h,'.r','markersize',5)
                        plot(tidx/plot_step,N_0l,'.g','markersize',5)
                        plot(tidx/plot_step,N_1l,'.k','markersize',5)
                        plot(tidx/plot_step,P,'.m','markersize',5)
                        drawnow
            
%             plot(tidx/plot_step,Nc_1h,'ob','markersize',5)
%             hold on
%             plot(tidx/plot_step,Nc_0h,'or','markersize',5)
%             plot(tidx/plot_step,Nc_0l,'og','markersize',5)
%             plot(tidx/plot_step,Nc_1l,'ok','markersize',5)
%             plot(tidx/plot_step,Pc,'om','markersize',5)
%             set(gca,'fontsize',18)
%             drawnow
%             
%             title('Blue: N_1h, Red: N_0h, Green: N_0low, Black: N_1low')
            
            
%             plot(tidx/plot_step,Pc/(Pc+2*P),'ob','markersize', 5)
%             hold on
%             drawnow
%             
%             title('Pc/(Pc+2*P)')
        else
            
            
            
            
        end
        
    end
    
    
end
