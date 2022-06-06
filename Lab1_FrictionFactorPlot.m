% Lab 1 Friction Factor Pipe Roughness
% Data collected 06 Jun 2022 in Upson I 112 Hydraulics Laboratory, UND.
%

% Last modified 06 Jun 2022 // MATLAB 9.9.0.1570001 (R2020b) Update 4 Win10
% David T. Wang (david.wang@und.edu) + Summer 2022 CE 423L Group XX Members

clear;

%% Set up friction factor calculation for laminar flow

ReiL = logspace(2,log10(4000),10);

fiL = 64./ReiL;

%% Set up Darcy-Weisbach friction factor calculation

rriT = [0e-6, 1e-6, 5e-6, 1e-5 5e-5, 1e-4, 5e-4, ... % relative roughness = e/D, lines to plot
    .001, .002, .005, .01, .015, .02, .03, .04, .05]; 

ReiT = logspace(log10(2100), 8, 100); % Reynolds number values to plot

% invsqrtf = zeros(length(ReiT), length(rriT)); 
fiT = zeros(length(ReiT), length(rriT)); 

for i = 1:length(ReiT),
    for j = 1:length(rriT),
        %invsqrtf(i,j) = -1.8 * log10( 6.9/ReiT(i) + (rriT(j)/3.7)^(10/9)
        %); % Haaland (1983) equation ... doesn't seem to work?
        %fiT(i,j) =0.0625 / ( log10( rriT(j)/3.7 + 5.74/(ReiT(i).^0.9) )^2 );
        fiT(i,j) = fzero(@(f) 1/sqrt(f)+2*log10(rriT(j)/3.71+2.51/(ReiT(i)*sqrt(f))),[eps,1]); % solve implicit Colebrook-White equation
%         fiT(i,j) = fzero(@(f) 1/sqrt(f)+4*log10(rriT(j)/3.71+1.255/(ReiT(i)*sqrt(f))),[eps,1]);
    end
end

%fiT = 1./sqrt(invsqrtf); % f for [Re e/D]

%fiT = fiT*4

%% Calculations on lab data

L = 0.524; %meters
D = 0.003; %meters
A = 0.785*D^2; %m2
nu = 1e-6; %m2/s

rhoHg = 13550; %kg/m3
rhoW = 1000; %kg/m3

data = readtable('Lab1_Data.xlsx'); % [Head_in, Time_sec, Vol_ml]

H = data.Head_in * 2.54/100; %m
Q = (data.Vol_ml ./ 1e6) ./ data.Time_sec ; %m3/s
V = Q/A; %m/s
Re = V*D/nu; %dimensionless

f = zeros(6,1);

for i=1:3,
    f(i) = 64./Re(i);
end

for i=4:6,
    f(i) = ( 2*9.81*D*H(i) )/( L*V(i)^2 ) * (rhoHg/rhoW - 1);
end

% Show results table
results = table(H, Re, f)

%% Best-fit line through lab data

rriF = linspace(0.001, 0.003, 11); %test e/D = 0.002, 0.0021, 0.0022, ...
fF = zeros(length(rriF), 6); %[e/D datapt]
for k=1:length(rriF), %k=test values of e/D
    for m=4:6, %m = columns of Reynolds number data
        fF(k,m) = fzero(@(f) 1/sqrt(f)+2*log10(rriF(k)/3.71+2.51/(Re(m)*sqrt(f))),[eps,1]);
    end
end

% Determine best-fit value by rastering across likely values of e/D
fF = fF(:,4:6);
bestfittests = [rriF' sum(bsxfun(@minus,f(4:6)',fF).^2, 2)*1e5] % [e/D, sumsq{resid}]
indexbest = find(bestfittests(:,2)==min(bestfittests(:,2))); %index of minimimum sumsq resdiuals
rri_best = bestfittests(indexbest,1);

fprintf(['Best-fit value of relative pipe roughness e/D = ' num2str(rri_best) '\n'])
fprintf(['Sum of squares of residuals in f @ best-fit value = ' num2str(bestfittests(indexbest,2)) '\n'])

% Make best-fit curve
f_best = zeros(100,1);
for i=1:100, %m = columns of Reynolds number data
    f_best(i) = fzero(@(f) 1/sqrt(f)+2*log10(rri_best/3.71+2.51/(ReiT(i)*sqrt(f))),[eps,1]); % solve implicit Colebrook-White equation
end

%% Plot Moody diagram with Fanning friction factor

figure(1); clf;

% Laminar Flow
loglog(ReiL, fiL, 'color', [0, 0.4470, 0.7410])        
text(2.9e3, 0.015, 'laminar flow', ...
    'FontSize', 8, ...
    'Color', [0, 0.4470, 0.7410])

% Transition Region
rectangle('Position', [2100, 0.0015, 4000, 1.10], ...
                'Curvature', 0.2, ...
                'FaceColor', [0, 0, 0, 0.1], ...
                'EdgeColor', [0, 0, 0, 0.0]);             
text(4.1e3, 0.021, 'transition region', ...
    'FontSize', 8, ...
    'Color', [0.5, 0.5, 0.5])      

hold on;

% Turbulent Flow
redMap = [linspace(0, 255, length(rriT))', zeros(length(rriT), 1), zeros(length(rriT), 1)]./255; %make red colormap for 15 lines
for j = 1:length(rriT),
    plot(ReiT, fiT(:,j), 'color', redMap(j,:))
    text(1e8, fiT(end,j), [' ' num2str(rriT(j),'%g')], ...
    'FontSize', 6, ...
    'Color', redMap(j,:), ...
    'VerticalAlign', 'middle', ...
    'HorizontalAlign', 'left')
end
text(1.5e6, 0.08, 'turbulent flow', ...
    'FontSize', 8, ...
    'Color', 'r')

% Best-fit e/D curve for turbulent flow data
plot(ReiT, f_best, 'k--')
text(1e8, f_best(end), [' ' num2str(rri_best,'%g')], ...
    'FontSize', 6, ...
    'Color', 'black', ...
    'FontWeight', 'bold', ...
    'VerticalAlign', 'middle', ...
    'HorizontalAlign', 'left')
text(1.1e7, 0.027, 'best-fit', ...
    'FontSize', 6, ...
    'FontWeight', 'bold', ...
    'Color', 'k')

%% Plot data

plot(Re, f, 'k.')
plot(Re, f, 'ko')

%% Format plot

ax = gca();
ax.YAxis.Exponent = 0;
% xtickformat('%.2f')
yticks([0.01, 0.015, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.1])
yticks = get(gca,'YTick');
set(gca,'YTickLabel',yticks);

ylim([0.006 0.1])
xlim([5e2 1e8])

grid on

xlabel('Reynolds number, $ Re = \frac{\rho\, VD}{\mu}$', 'interpreter', 'latex')
ylabel('Friction factor, $f$', 'interpreter', 'latex')

title('Moody diagram for Lab 1')

