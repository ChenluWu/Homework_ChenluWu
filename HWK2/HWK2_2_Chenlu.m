clear all;
clc;

rng(1234);

n = 1000;
t=40;
beta = 0.99^(1/12);
sigma_u = sqrt(0.2);
sigma_eps = sqrt(0.2);

ln_eps = normrnd(0,sigma_eps,n,t);
eps = exp(ln_eps);

ln_u = normrnd(0,sigma_u,n,1);
u = exp(ln_u);

z = exp(-sigma_u^2/2) * u;

SeaComponent = [0.863,0.727,0.932;0.691,0.381,0.845;1.151,1.303,1.076;
                1.140,1.280,1.070;1.094,1.188,1.047;1.060,1.119,1.030;
                1.037,1.073,1.018;1.037,1.073,1.018;1.037,1.073,1.018;
                1.002,1.004,1.001;0.968,0.935,0.984;0.921,0.843,0.961];
            
StoSeaComponent = [0.085 0.171 0.043;0.068 0.137 0.034;0.290 0.580 0.145;
                   0.283 0.567 0.142;0.273 0.546 0.137;0.273 0.546 0.137;
                   0.239 0.478 0.119;0.205 0.410 0.102;0.188 0.376 0.094;
                   0.188 0.376 0.094;0.171 0.341 0.085;0.137 0.273 0.068];
       
ln_eps_m = normrnd(0,sqrt(StoSeaComponent));
eps_m = exp(ln_eps_m);
StoSeaComponent = exp(-StoSeaComponent/2) .* eps_m;

NonSeaComponent = exp(-sigma_eps^2/2) * eps;

theta = 0.3224;
kappa = ((1-theta)/0.3)/(28.5 * 30/7);
v = 1;
 
results_all = zeros(1000,3);
results_conrmv = zeros(1000,3);
results_labrmv = zeros(1000,3);

for degree = 1:3
   
    u_life = zeros(1000,1);
    u_life_conrmv = zeros(1000,1);
    u_life_labrmv = zeros(1000,1);
    
    for t = 1:40
       
        u_year = zeros(1000,1);
        u_year_conrmv = zeros(1000,1);
        u_year_labrmv = zeros(1000,1);
        
        for m = 1:12
            c_mt = z .* SeaComponent(m,degree).*StoSeaComponent(m,degree) .* NonSeaComponent(:,t);
            c_searmv = z .* NonSeaComponent(:,t);
            %c_nonsearmv = z .* SeaComponent(m,degree).*StoSeaComponent(m,degree);
          
            h_mt =  0.5 * (z .* SeaComponent(m,degree).*StoSeaComponent(m,degree) .* NonSeaComponent(:,t));
            h_searmv =  0.5 * (z .* NonSeaComponent(:,t));
            %h_nonsearmv = z .* SeaComponent(m,degree).*StoSeaComponent(m,degree);
            
            u_cmt = log(c_mt) - kappa*(h_mt.^(1+1/v)/(1+1/v));
            u_conrmv = log(c_searmv) - kappa*(h_mt.^(1+1/v)/(1+1/v));
            u_labrmv = log(c_mt) - kappa*(h_searmv.^(1+1/v)/(1+1/v));
            
            u_year = u_year + beta^(m-1) * u_cmt;
            u_year_conrmv = u_year_conrmv + beta^(m-1) * u_conrmv;
            u_year_labrmv = u_year_labrmv + beta^(m-1) * u_labrmv;
        end
        
        u_life = u_life + beta^(12*t) * u_year;
        u_life_conrmv = u_life_conrmv + beta^(12*t) * u_year_conrmv;
        u_life_labrmv = u_life_labrmv + beta^(12*t) * u_year_labrmv;
    end
    
    results_all(:,degree) = u_life;
    results_conrmv(:,degree) = u_life_conrmv;
    results_labrmv(:,degree) = u_life_labrmv;
end

A = 0;
B = 0;
for m = 1:12
    B = B + beta^(m-1);
end
for t = 1:40
    A = A + beta^(12*t)*B;
end

Gain_conrmv = exp((results_conrmv - results_all)/A)-1;
Gain_labrmv = exp((results_labrmv - results_all)/A)-1;

avg_conrmv = mean(Gain_conrmv)
avg_labrmv = mean(Gain_labrmv)
std_conrmv = std(Gain_conrmv)
std_labrmv = std(Gain_labrmv)


