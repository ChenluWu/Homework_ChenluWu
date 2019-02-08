clear all;
clc;

rng(1234);

eta = 4;

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
            
NonSeaComponent = exp(-sigma_eps^2/2) * eps;

results_all = zeros(1000,3);
results_searmv = zeros(1000,3);
results_nonsearmv = zeros(1000,3);

for degree = 1:3
   
    u_life = zeros(1000,1);
    u_life_searmv = zeros(1000,1);
    u_life_nonsearmv = zeros(1000,1);
    
    for t = 1:40
       
        u_year = zeros(1000,1);
        u_year_searmv = zeros(1000,1);
        u_year_nonsearmv = zeros(1000,1);
        
        for m = 1:12
            c_mt = z .* SeaComponent(m,degree) .* NonSeaComponent(:,t);
            c_searmv = z .* NonSeaComponent(:,t);
            c_nonsearmv = z .* SeaComponent(m,degree);
            
            u_cmt = c_mt.^(1-eta) / (1-eta);
            u_searmv = c_searmv.^(1-eta) / (1-eta);
            u_nonsearmv = c_nonsearmv.^(1-eta) / (1-eta);
%             u_cmt = log(c_mt);
%             u_searmv = log(c_searmv);
%             u_nonsearmv = log(c_nonsearmv);
            
            u_year = u_year + beta^(m-1) * u_cmt;
            u_year_searmv = u_year_searmv + beta^(m-1) * u_searmv;
            u_year_nonsearmv = u_year_nonsearmv + beta^(m-1) * u_nonsearmv;
        end
        
        u_life = u_life + beta^(12*t) * u_year;
        u_life_searmv = u_life_searmv + beta^(12*t) * u_year_searmv;
        u_life_nonsearmv = u_life_nonsearmv + beta^(12*t) * u_year_nonsearmv;
    end
    
    results_all(:,degree) = u_life;
    results_searmv(:,degree) = u_life_searmv;
    results_nonsearmv(:,degree) = u_life_nonsearmv;
end


Gain_searmv2 = (results_searmv./results_all).^(1/(1-eta))-1;
Gain_nonsearmv2 = (results_nonsearmv./results_all).^(1/(1-eta))-1;

avg_searmv2 = mean(Gain_searmv2)
avg_nonsearmv2 = mean(Gain_nonsearmv2)
std_searmv2 = std(Gain_searmv2)
std_nonsearmv2 = std(Gain_nonsearmv2)


histogram(Gain_nonsearmv2)
xlabel('g');ylabel('number of hh');
            
            