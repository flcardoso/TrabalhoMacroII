% Trabalho Macroeconomia II - Small Open Economy - Carlos Viana Carvalho
% SOE_Model
% Fernando Cardoso
% PUC-Rio - 2H 2018

%----------------------------------------------------------------
%% Housekeeping
%----------------------------------------------------------------



%----------------------------------------------------------------
%% Defining Variables
%----------------------------------------------------------------

var pi $\hat\pi$ mu $\hat\mu$ mc $\widehat{mc}$ w $\hat w$ s $\hat s$ 
ps $\hat p^*$ p $\hat p$ a $a$ x $\hat x$ g $\hat g$ 
y $\hat y$ pi_m $\overline\pi$ l $\hat l$ z $\hat z$ phi $\hat \phi$ 
i $\hat i$ c $\hat c$ nu $\nu$ 
pis $\pi^*$ is $\hat i^*$
;

varexo err_mu $\varepsilon^\mu$ err_g $\varepsilon^g$ err_a $\varepsilon^a$
err_nu $\varepsilon^\nu$ err_phi $\varepsilon^\phi$ err_pi $\varepsilon^\pi$ 
err_is $\varepsilon^{i^*}$
;

parameters 
beta $\beta$ h sigma $\sigma$ varphi $\varphi$ rho $\rho$ i_ $\overline i$
delta $\delta$ lambda $\lambda$ gamma $\gamma$ theta $\theta$ 
rho_i $\rho_i$ phi_pi $\phi_\pi$ phi_c $\phi_c$ phi_dc $\phi_{\Delta c}$
phi_2 $\phi_2$ chi $\chi$

rho_g $\rho_g$      sigma_g $\sigma_g$      g_ $\overline g$
rho_a $\rho_a$      sigma_a $\sigma_a$      a_ $\overline a$
rho_pi $\rho_\pi$   sigma_pi $\sigma_\pi$ 
rho_phi $\rho_\phi$ sigma_phi $\sigma_\phi$ phi_ $\overline \phi$
rho_mu $\rho_\mu$   sigma_mu $\sigma_\mu$   mu_ $\overline\mu$
rho_is $\rho_{i^*}$ sigma_is $\sigma_{i^*}$ is_ $\overline i^*$
rho_nu $\rho_\nu$   sigma_nu $\sigma_\nu$   

vartheta_l $\vartheta_l$ vartheta_z $\vartheta_z$
varrho_pi $\varrho_\pi$ kappa_m $\kappa_m$
;

%----------------------------------------------------------------
%% Calibration
%----------------------------------------------------------------

rho     = 0.02/4;
h       = .25;      % Habit Persistence 0<=h<1
sigma   = 1;        % Inverse Intertemp Subs/Risk Aversion
varphi  = 5;        % Inverse Frisch
delta   = 1/3;      % Imports
lambda  = 1/4;      % Calvo
gamma   = 0.25;     % Inflation Inertia
rho_i   = .25;        % Interest Rate Inertia (MPR)
phi_pi  = 1.5;      % Inflation (MPR)
phi_c   = 0;        % Consumption (MPR)
phi_dc  = .25;        % Delta Consumption (MPR)
phi_2   = .2;       % Exchange Rate (MPR)
chi     = .001;     % Exchange Rate Imbalance Coeff (UIP)
theta   = 9 ;       % S-S Elasticity of Subst.

rho_g = .8; 
sigma_g = 0.0025;
g_ = 0;

rho_a  = .8;
sigma_a = 0.01;
a_ = 0;

rho_pi  = .8;
sigma_pi  = 0.2;

rho_phi   = .8;
sigma_phi  = 0.01;
phi_ = 0;

rho_mu   = .8;
sigma_mu  = 0.2;

rho_is   = .8;
sigma_is  = 0.0025;


rho_nu   = .3;
sigma_nu  = 0.0025;



% Functional Parameters
beta        = exp(-rho);
i_          = rho  ;
is_         = rho  ;
mu_         = log(theta/(theta-1)); % Steady State Markup

vartheta_l  = delta*(log((1-delta)/delta));
vartheta_z  = log(delta/(1-delta)) + vartheta_l;
varrho_pi   = 1/(1+beta*gamma);
kappa_m       = (lambda/(1-lambda))*(1-(1-lambda)*beta);


%----------------------------------------------------------------
%% Model
%----------------------------------------------------------------

model(linear);
    
    [name = 'Phillips Curve',type='endogenous']
    pi = varrho_pi * (gamma*pi(-1) + beta*pi(1) + kappa_m*((mu-mu_) + mc)) ;
    
    [name = '(Real) Marginal Cost',type='endogenous']
    mc = (1-delta)*w + delta*(s + ps) - p - a ;
    
    [name = 'Euler Equation',type='endogenous']
    x = x(1) - sigma^(-1)*((i-i_) - pi(1) + rho_g*(g-g_)) ;
    
    [name = 'Habit',type='endogenous']
    x = (1/(1-h))*(c - h*c(-1)) ;
    
    [name = 'Mkt Clearing',type='endogenous']
    y = c;
    
    [name = 'Monetary Policy Rule',type='endogenous']
    (i - i_) = rho_i*(i(-1) - i_)
        + (1-rho_i)*(phi_pi*(p(2)-p(-1)-pi_m) + phi_c*c + phi_dc*(c-c(-4)) 
        + phi_2*(s - s(-2))) + nu ;
        
    [name = 'Labor',type='endogenous']
    l = y - a + delta*(s + ps - w) ;
    
    [name = 'Wages',type='endogenous']
    w = p + varphi*l + sigma*x;
    
    [name = 'Imports',type='endogenous']
    z = y - a + (1-delta)*(w - (s + ps)) ;
    
    [name = 'UIP',type='endogenous']
    (i - i_) - (is - is_) = s(1) - s + phi - chi*(s + ps - p) ;
     
% Exogenous Processes
    
    [type='exogenous']
    mu-mu_ =  rho_mu*(mu(-1)-mu_) + sigma_mu*err_mu ;
    
    [type='exogenous']
    g-g_ = rho_g*(g(-1)-g_) + sigma_g*err_g ;
    
    [type='exogenous']
    a-a_ = rho_a*(a(-1)-a_) + sigma_a*err_a ;
    
    [type='exogenous']
    nu = rho_nu*nu(-1) + sigma_nu*err_nu ;
    
    [type='exogenous']
    phi-phi_ = rho_phi*(phi(-1)-phi_) + sigma_phi*err_phi ;
    
    [type='exogenous']
    pis = rho_pi*pis(-1) + sigma_pi*err_pi ;
    
    [type='exogenous']
    is-is_ = rho_is*(is(-1)-is_) + sigma_is*err_is ;
    
% Definitions

    [name = 'Inflation',type='endogenous']
    pi = p - p(-1) ;
    
    [name = 'Foreign Inflation',type='exogenous']
    pis = ps - ps(-1) ;
    
    [name = 'Inflation Target',type='exogenous']
    pi_m = 0 ;
    
    
end;

%----------------------------------------------------------------
%% Initializing Values
%----------------------------------------------------------------

initval;

pi  = 0;
mu  = mu_;
mc  = 0 ;
w   = 0 ;
s   = 0 ;
ps  = 0 ;
p   = 0 ;
a   = 0 ;
x   = 0 ;
g   = 0 ;
y   = 0 ;
l   = 0 ;
z   = 0 ;
phi = 0 ;
i   = rho ;
c   = 0 ;
nu  = 0 ;
pis = 0 ;
is  = rho ;

end;

%----------------------------------------------------------------
%% Shocks
%----------------------------------------------------------------

shocks;
var err_mu  = sigma_mu^2;
var err_g   = sigma_g^2;
var err_a   = sigma_a^2;
var err_nu  = sigma_nu^2;
var err_phi = sigma_phi^2;
var err_pi  = sigma_pi^2;
var err_is  = sigma_is^2;
end;


%----------------------------------------------------------------
%% Run
%----------------------------------------------------------------
write_latex_static_model;
write_latex_dynamic_model;
write_latex_original_model;
resid;
model_diagnostics;
steady;
check;
stoch_simul(order=1,periods=1000,irf=25,tex);
