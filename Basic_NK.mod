addpath c:\Users\mahon\dynare\4.6.4\matlab;

var pi 
    y_gap 
    y_nat
    y
    yhat
    r_nat
    r_real
    i
    n
    nu
    a
    z
    p
    c
    w
    w_real
    mu
    mu_hat
    r_real_ann
    i_ann
    r_nat_ann
    pi_ann
;

varexo eps_a
       eps_nu
       eps_z
;

parameters alpha
    beta
    rho_a
    rho_nu
    rho_z
    sigma
    varphi
    phi_y
    phi_pi
    eta
    epsilon
    theta
    chi
;

chi=1
sigma=1.5;
phi_pi=1;
phi_y=1.1;
theta=3/4;
rho_nu=0.5;
rho_z=0.5;
rho_a=0.5;
beta=0.99;
eta=5;
alpha=1/4;
epsilon=10;

model(linear);
//Composite parameters
#Omega=(1-alpha)/(1-alpha-alpha*epsilon);
#lambda=(1-theta)*(1-beta*theta)/theta*Omega;
#psi_n_ya=(1+varphi)/(sigma*(1-alpha)+varphi+alpha);
#kappa=lambda*(sigma+(varphi+alpha)/(1-alpha));

[name='NKPC']
pi=beta*pi(+1)-kappa*y_gap;

[name='IS Curve']
y_gap=-1/sigma*(i-pi(+1)-r_nat)+y_gap(+1);

[name='interest rate rule']
i=phi_pi*pi+phi_y*yhat+nu;

[name='Natural rate']
r_nat=-sigma*psi_n_ya*(1-rho_a)*a+(1-rho_z)*z;

[name='real rate']
r_real=i-pi(+1);

[name='Natty output']
y_nat=psi_n_ya*a;

[name='defn output gap']
y_gap=y-y_nat;

[name='MP shock']
nu=rho_nu*nu(-1)+eps_nu;

[name='TFP shock']
a=rho_a*a(-1)+eps_a;

[name='production fn']
y=a+(1-alpha)*n;

[name='pref shock']
z = rho_z*z(-1)-eps_z;

[name='Output deviation from steady state']
yhat=y-steady_state(y);

[name='Price level']
pi=p-p(-1);

[name='market clearing']
c=y;

[name='Labor supply FOC']
w-p=sigma*c+chi*n;

[name='defn real wage']
w_real=w-p;

[name='average price markup']
mu=-(sigma+(varphi+alpha)/(1-alpha))*y_nat+(1+varphi)/(1-alpha)*a;

[name='Average markup']
mu_hat=-(sigma+(varphi+alpha)/(1-alpha))*y_gap;

[name='Annualized nominal interest']
i_ann=4*i;

[name='Annualized real interest']
r_real_ann=4*r_real;

[name='Annualized nat interest']
r_nat_ann=4*r_nat;

[name='Annualized inflation']
pi_ann=4*pi;

end;

shocks;
var eps_nu=0.25^2;
end;

resid(1);
steady;
check;

stoch_simul(order = 1,irf=15) y_gap pi_ann y n w_real p i_ann r_real_ann nu;



















































