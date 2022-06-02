%% Figures from "Asymptotic Behaviour of Optimal Vaccination Policies"


%% Default Parameter Values
beta1 = [1 2; 2 4];
beta2 = [0.5 1; 1 2];
beta3 = beta2;
beta4 = [0.25 0.5; 0.5 1];

N2 = 1;
mu1 = [1,1];
mu2 = [1,1];


p2 = 1;
Istar = 0.01;
A = 1;
%% Fig 1
close all
%Initialise Objective values
Hvals = zeros(100,2);
count = 1;
for eps = 0.01:0.01:1
    %Initial conditions
    y0 = [eps,0,0,0,0,0,0,1-Istar,Istar,0,0,0,0,0];
    %Set vulnerability and population size
    p1 = 1/eps;
    N1 = eps;
    % Define vaccination policies
    U = @(t) [0,A*(1 - heaviside(t-1/A))];
    Utilde = @(t) [A*(1-heaviside(t -eps/A)),A*(1-heaviside((t-eps/A)*(t-1/A)))]; 
    %Solve equations and log objective values
    [t,y] = ode23(@(t,y) SIRODE(y,U(t),beta1,beta2,beta3,beta4,mu1,mu2,[N1,N2]),[0,100],y0);
    Hvals(count,1) = y(length(y),3)*p1 + y(length(y),6)*p1 + p2*y(length(y),10)  + p2*y(length(y),13) ;
    [t2,y2] = ode23(@(t,y) SIRODE(y,Utilde(t),beta1,beta2,beta3,beta4,mu1,mu2,[N1,N2]),[0,100],y0);
    Hvals(count,2) =y2(length(y2),3)*p1 + y2(length(y2),6)*p1 +p2*y2(length(y2),10)  + p2*y2(length(y2),13) ;

    count = count + 1;
end
%Create figure
figure;
semilogx(0.01:0.01:1,Hvals,'linewidth',3)
xlabel('$\epsilon$','interpreter','latex')
ylabel('Objective value')
legend('$H(U)$','$H(\tilde{U})$','interpreter','latex','location','northwest')
set(gca,'fontsize',18)
print('-depsc','Fig1.eps')

%% Fig 2
%Initialise values
EpsstarvalsBonusN = zeros(101,3);
pnewvals = [0,0.5,1];
Inewvals = -6:0.05:-1;
for  i = 1:101
    for j = 1:3
        %Set Istar and p2
        Istar = 10^(Inewvals(i));
        p2 = pnewvals(j);
        %Set trigger value for finding epsilon star
        Found = 0;
        %Set initial epsilon value
        epsexp = -3;
        while Found == 0
            %Iterate through different epsilon values until objective
            %condition is satisfied
            eps = 10^(epsexp);
            y0 = [eps,0,0,0,0,0,0,1-Istar,Istar,0,0,0,0,0];
            p1 = 1/eps;
            N1 = eps;
            U = @(t) [0,A*(1 - heaviside(t-1/A))];
            Utilde = @(t) [A*(1-heaviside(t -eps/A)),A*(1-heaviside((t-eps/A)*(t-1/A)))]; 
            [t,y] = ode23(@(t,y) SIRODE(y,U(t),beta1,beta2,beta3,beta4,mu1,mu2,[N1,N2]),[0,100],y0);
            Val1 = y(length(y),3)*p1 + y(length(y),6)*p1 + p2*y(length(y),10)  + p2*y(length(y),13) ;
            [t2,y2] = ode23(@(t,y) SIRODE(y,Utilde(t),beta1,beta2,beta3,beta4,mu1,mu2,[N1,N2]),[0,100],y0);
            Val2 =y2(length(y2),3)*p1 + y2(length(y2),6)*p1 +p2*y2(length(y2),10)  + p2*y2(length(y2),13) ;
            if Val1 < Val2
                Found = 1;
                EpsstarvalsBonusN(i,j) = epsexp;
            else
                epsexp = epsexp + 0.01;
            end
            if epsexp > 0
                Found = 1;
                EpsstarvalsBonusN(i,j) = 0;
            end
        end
        disp(i)
        disp(j)
    end
end
%Create figure 
close all
figure;
plot(-6:0.05:-1,EpsstarvalsBonusN(:,1),'Linewidth',3)

hold on

plot(-6:0.05:-1,EpsstarvalsBonusN(:,2),'Linewidth',3)

plot(-6:0.05:-1,EpsstarvalsBonusN(:,3),'Linewidth',3)
xlabel('$\log_{10}(I^*)$','Interpreter','latex')
ylabel('$\log_{10}(\epsilon^*(I^*,p^*))$','Interpreter','latex')
legend('$p^* = 0$','$p^*=0.5$','$p^* = 1$','interpreter','latex','location','northwest')
set(gca,'fontsize',18)
print('-depsc','Fig2.eps')

%% Fig 3
%The process of creating Fig 2 is repeated but for different values of
%beta2, beta3 and beta4.
beta2n = zeros(2,2);
beta3n = zeros(2,2);
beta4n = zeros(2,2);
Inewnewvals = -5:0.05:1;
EpsstarvalsBonusNN = zeros(81,3);
for  i = 1:81
    for j = 1:3
        Istar = 10^(Inewnewvals(i));
        p2 = pnewvals(j);
        Found = 0;
        epsexp = -1;
        while Found == 0
            eps = 10^(epsexp);
            y0 = [eps,0,0,0,0,0,0,1-Istar,Istar,0,0,0,0,0];
            p1 = 1/eps;
            N1 = eps;
            U = @(t) [0,A*(1 - heaviside(t-1/A))];
            Utilde = @(t) [A*(1-heaviside(t -eps/A)),A*(1-heaviside((t-eps/A)*(t-1/A)))]; 
            [t,y] = ode23(@(t,y) SIRODE(y,U(t),0.0001*beta1,beta2n,beta3n,beta4n,mu1,mu2,[N1,N2]),[0,100],10000*y0);
            Val1 = y(length(y),3)*p1 + y(length(y),6)*p1 + p2*y(length(y),10)  + p2*y(length(y),13) ;
            [t2,y2] = ode23(@(t,y) SIRODE(y,Utilde(t),0.0001*beta1,beta2n,beta3n,beta4n,mu1,mu2,[N1,N2]),[0,100],10000*y0);
            Val2 =y2(length(y2),3)*p1 + y2(length(y2),6)*p1 +p2*y2(length(y2),10)  + p2*y2(length(y2),13) ;
            if Val1 < Val2
                Found = 1;
                EpsstarvalsBonusNN(i,j) = epsexp;
            else
                epsexp = epsexp + 0.01;
            end
            if epsexp > 0
                Found = 1;
                EpsstarvalsBonusNN(i,j) = 0;
            end
        end
        disp(i)
        disp(j)
    end
end
close all
figure;
plot(-5:0.05:-1,EpsstarvalsBonusNN(:,1),'Linewidth',3)

hold on

plot(-5:0.05:-1,EpsstarvalsBonusNN(:,2),'Linewidth',3)

plot(-5:0.05:-1,EpsstarvalsBonusNN(:,3),'Linewidth',3)
xlabel('$\log_{10}(I^*)$','Interpreter','latex')
ylabel('$\log_{10}(\epsilon^*(I^*,p^*))$','Interpreter','latex')
legend('$p^* = 0$','$p^*=0.5$','$p^* = 1$','interpreter','latex','location','northwest')
set(gca,'fontsize',18)
print('-depsc','Snapshot2.eps')
%% UK Data
%Input contact matrix and population sizes for UK
muUK = [1.7112383721011	0.80335911324295	0.388876091777961	0.28572486821271	0.33696343169784	0.628688406705816	0.725146694937456	0.91886117781732	0.459706790825267	0.174630433861908	0.194182456536066	0.131538050213491	0.074698377482648	0.043368114629815	0.037356498709477	0.00720339776988
0.681647122643733	4.03425249578426	0.714526756407179	0.209567190507402	0.100595729199919	0.426362993332467	0.739049366869117	0.815510692470327	0.650366840671183	0.290179894310514	0.227030412002046	0.073027448998374	0.092351479235614	0.058038016270078	0.039535667745856	0.013488696707224
0.317547244688793	2.00878867776471	4.54147098246679	0.817541122335809	0.141833625971166	0.2763371779291	0.452574419683193	0.880121048871909	0.865244809195991	0.340404267624021	0.185715793854421	0.094278921165269	0.078321777377439	0.032326892507348	0.040480391795631	0.052535091386635
0.17644454896498	0.496186583452507	2.27907535644973	7.87311346700604	1.32394506650135	0.632567776147492	0.401289024506748	0.984423203294234	0.900470374284428	0.952732012220013	0.600994592671502	0.219218070443359	0.087133096664032	0.06080025969978	0.015326110133108	0.00062601530237
0.320395511005378	0.247639067952491	0.270288567695945	2.02438212689977	2.42463164525645	1.29474541305459	0.897445389158323	1.08244889010628	0.691615914429445	0.977956249583577	0.361817345222361	0.387602972759624	0.120246542344494	0.037727348033907	0.022461150830349	0.05130199607221
0.841884650035674	0.456517309028613	0.14265645815947	0.758164265159102	1.73113538947572	2.40508711992341	1.4735225350913	0.970026692858644	1.28523343334207	0.932933999281656	0.916492247994984	0.296863368186688	0.149136506241853	0.055388604835446	0.017096676595794	6.42930854276393E-06
0.483288236190419	0.743776350691811	0.419156281817357	0.365889257682038	0.962718216887447	1.1548217423198	1.994987645559	1.63874229856943	0.862762983625699	0.921011769842121	0.502912299674096	0.501537238162951	0.142881958039452	0.058254188722895	0.009742403801571	0.024821296169974
0.543164548546523	1.09108499332856	0.752279269408771	0.653622161800438	0.63007926520706	0.923993645161487	0.98227639250514	1.88743812512276	1.35526325981955	0.961319152453284	0.694139106582885	0.411206840390126	0.23073157825577	0.118413244206554	0.111418751170672	0.038397813830407
0.236055712743079	0.86979767264152	0.920948050875685	0.790691075005073	0.769328643682645	0.773330971972483	1.03684752420129	1.42882015012761	1.74976440334236	1.27132404390175	0.612171130433804	0.454891280980261	0.147189038602815	0.035496980789485	0.061634354208815	0.00131900989626
0.14272927635623	0.153915964344865	0.643726812145509	1.32353818131216	0.933062511412995	0.972016047384488	0.981871920905164	1.15890480677256	1.45259702268274	2.18342259262626	0.917529049836053	0.38463354070245	0.21200436872059	0.099062676978984	0.062858793760666	0.066674510851858
0.150004662749013	0.130060304363535	0.277435285895216	0.64485755755201	0.729446082075889	0.973979370522516	0.779171735079313	0.975921162039635	1.38588749354254	1.17618021340199	1.0186201223853	0.615888029724108	0.266016570369211	0.143804724804949	0.12872166073862	0.012295309427609
0.172105230351602	0.207377040650624	0.223204341110905	0.473102962845615	0.915592368793147	1.11497708596801	0.954317923986381	0.728159820079553	0.936304871375322	0.660008299552242	0.851975025307012	1.24584809916084	0.508881854587708	0.207374089669334	0.12379844974119	0.109684868083495
0.021279441241425	0.223166220637787	0.113044110753442	0.059927343497398	0.422373830455612	0.392356673584836	0.471942054072193	0.583251196547724	0.500815638326449	0.374818282177818	0.388842697809227	0.560892008981676	0.739544544280317	0.156816206433602	0.128423947908426	0.081782931726588
0.0649151909962	0.245142970109999	0.216258172677734	0.213251941735953	0.123207866929172	0.501809975315919	0.344794734492308	0.593423629506972	0.633226483056339	0.457567371030771	0.405526646679515	0.598569703290443	0.37158315387733	0.726331606939168	0.191239431540953	0.056493947502016
0.072816898130484	0.048429811137548	0.254360693525395	0.487877319805465	0.514897333716445	0.439531837719742	0.242898620478756	0.138245804193342	0.785103900484758	0.441180313988881	0.371836287529357	0.467690646103474	0.54048193431852	0.694249409983356	0.398546702699878	0.339001888253868
0.020607384349943	0.000688316176792	0.077566053602813	0.170458170021847	0.087028501848366	0.038604883601546	0.170578859512259	0.215447768007773	0.268900161890028	0.461031245301225	0.429559106046982	0.000198892237717	0.466582938224367	0.000338156571627	0.302541461262816	0.738331961063434
];
PopUK = [  3924	  4120	  3956	  3686	  4075	  4484	  4707	  4588	  4308	  4296	  4635	  4539	  3905	  3382	  3388	  2442	  1737	  1078	   491	   130	   16
];
%Group together the oldest age groups
PopUK = [PopUK(1:15) sum(PopUK(16:21))];

%% Fig 4
%Create next generation matrix and scale so its leading eigenvalue is 4.
BetaUK = muUK./PopUK;
ScaleUK = 4/eigs((PopUK/sum(PopUK)).*BetaUK',1);
BetaUK = 4*(BetaUK/eigs((PopUK/sum(PopUK)).*BetaUK',1));

%Create figure
close all
figure;
ax = heatmap(BetaUK.*(PopUK/sum(PopUK)),'CellLabelColor','none');
labels = categorical({'0-4','5-9','10-14','15-19','20-24','25-29','30-34','35-39','40-44','45-49','50-54','55-59','60-64','65-69','70-74','75+'});
labels = reordercats(labels,{'0-4','5-9','10-14','15-19','20-24','25-29','30-34','35-39','40-44','45-49','50-54','55-59','60-64','65-69','70-74','75+'});


annotation('textarrow',[0.98,0.98],[0.5,0.5],'string','$R_{ij}$', ...
      'HeadStyle','none','LineStyle','none','HorizontalAlignment','center','TextRotation',90,'interpreter','latex','fontsize',18);
ax.XData = labels;
ax.YData = labels;
ax.XLabel = "Age Group (years)";
ax.YLabel = "Age Group (years)";
ax.NodeChildren(3).XAxis.Label.Interpreter = 'latex';
ax.NodeChildren(3).YAxis.Label.Interpreter = 'latex';

ax.NodeChildren(3).Title.Interpreter = 'latex';
set(gca,'fontsize',18)
print('-depsc','Fig4.eps')
%% Fig 5
close all
%Set epsilon exponent
exponent = -1;
%Initialise comparisons
yvals = zeros(length(0.01:0.01:5),1);
actualdecs = zeros(length(0.01:0.01:5),1);
count = 1;
for beta = 0.01:0.01:5
    %Solve for no vaccination
    U= @(t) 0;
    [t3,y3] = ode23(@(t,y) SIRODE(y,U(t),beta,0.5*beta,0.5*beta,0.25*beta,1,1,1),[0,10000],[(1-10^(-4)),10^(-4),0,0,0,0,0]);
    Rvals = [y3(length(y3),3),0];
    M = zeros(2,2);
    SIndices = [1,2];
    y0 = [1-10^(-4),0];
    %Solve for linear approximation to vaccination
    BetaDash = [beta, 0.5*beta;0.5*beta,0.25*beta];
    for i = 1:2
        for j = 1:2
            M(i,j) = (eq(i,j) - y0(SIndices(i)) * BetaDash(i,j)* exp(-dot(BetaDash(i,:),Rvals)) ) / (1 - exp(-dot(BetaDash(i,:),Rvals)));


        end
    end
    %disp(M)
    p = ones(2,1);
    x = M'\p;

    yUK = zeros(1,1);
    for i = 1:1
        yUK(i) = (1-10^(-4))*(x(i+1) - x(i));
    end
    %Create final linear approximation
    yvals(count) = yUK(1);
    ycand = (1-10^(-4))*((-1 + ((0.5*beta*y0(1)*exp(-beta*Rvals(1)))/(1-exp(-beta*Rvals(1))))*(1-exp(-0.5*beta*Rvals(1))))*((1-exp(-beta*Rvals(1)))/(1-beta*y0(1)*exp(-beta*Rvals(1)))) + (1-exp(-beta*0.5*Rvals(1))));
    disp(ycand)
    disp(yvals(count))
    disp('------------------------------')
    U = @(t) (1- heaviside(t-10^(exponent)));
    %Get real decrease in deaths
    [t4,y4] = ode89(@(t,y) SIRODE(y,U(t),beta,0.5*beta,0.5*beta,0.25*beta,1,1,1),[0,10000],[1-10^(-4),10^(-4),0,0,0,0,0]);
    actualdecs(count) = (y4(length(y4),3) + y4(length(y4),6) -y3(length(y3),3) -y3(length(y3),6));
    count = count + 1;
%Repeat for epsilon = 0.01 = 10^{-2}
endexponent = -2;
yvals2 = zeros(length(0.01:0.01:5),1);
actualdecs2 = zeros(length(0.01:0.01:5),1);
count = 1;
for beta = 0.01:0.01:5
    U= @(t) 0;
    [t3,y3] = ode23(@(t,y) SIRODE(y,U(t),beta,0.5*beta,0.5*beta,0.25*beta,1,1,1),[0,10000],[(1-10^(-4)),10^(-4),0,0,0,0,0]);
    Rvals = [y3(length(y3),3),0];
    M = zeros(2,2);
    SIndices = [1,2];
    y0 = [1-10^(-4),0];
    BetaDash = [beta, 0.5*beta;0.5*beta,0.25*beta];
    for i = 1:2
        for j = 1:2
            M(i,j) = (eq(i,j) - y0(SIndices(i)) * BetaDash(i,j)* exp(-dot(BetaDash(i,:),Rvals)) ) / (1 - exp(-dot(BetaDash(i,:),Rvals)));


        end
    end
    %disp(M)
    p = ones(2,1);
    x = M'\p;

    yUK = zeros(1,1);
    for i = 1:1
        yUK(i) = (1-10^(-4))*(x(i+1) - x(i));
    end
    yvals2(count) = yUK(1);
    ycand = (1-10^(-4))*((-1 + ((0.5*beta*y0(1)*exp(-beta*Rvals(1)))/(1-exp(-beta*Rvals(1))))*(1-exp(-0.5*beta*Rvals(1))))*((1-exp(-beta*Rvals(1)))/(1-beta*y0(1)*exp(-beta*Rvals(1)))) + (1-exp(-beta*0.5*Rvals(1))));
    disp(ycand)
    disp(yvals2(count))
    disp('------------------------------')
    U = @(t) (1- heaviside(t-10^(exponent)));
    [t4,y4] = ode89(@(t,y) SIRODE(y,U(t),beta,0.5*beta,0.5*beta,0.25*beta,1,1,1),[0,10000],[1-10^(-4),10^(-4),0,0,0,0,0]);
    actualdecs2(count) = (y4(length(y4),3) + y4(length(y4),6) -y3(length(y3),3) -y3(length(y3),6));
    count = count + 1;
end
%Create figure
figure;

tiledlayout(2,1)
nexttile
plot(0.01:0.01:5,yvals*10^(-1),'Linewidth',2)
hold on
plot(0.01:0.01:5,actualdecs,'Linewidth',2)
xticks(0:10)
legend('Predicted $\rho_1$','Actual $\rho_1$','interpreter','latex','location','southeast')
xlabel('$\beta$','interpreter','latex')
ylabel('$\rho_1$','Interpreter','latex')
title('(a): $\epsilon = 0.1$','interpreter','latex')
set(gca,'fontsize',18)
nexttile
plot(0.01:0.01:5,yvals2*10^(-2),'Linewidth',2)
hold on
plot(0.01:0.01:5,actualdecs2,'Linewidth',2)
xticks(0:10)
legend('Predicted $\rho_1$','Actual $\rho_1$','interpreter','latex','location','southeast')
xlabel('$\beta$','interpreter','latex')
ylabel('$\rho_1$','Interpreter','latex')
title('(b): $\epsilon = 0.01$','interpreter','latex')
set(gca,'fontsize',18)
print('-depsc','Fig5.eps')
%% Fig 6
%Initialise vector for initial conditions
y0 = zeros(16*7,1);
%These indices relate to the susceptible groups
Indices = 1:7:(7*15+1);
%Input initial conditions
y0(Indices) = (1-10^(-4))*(PopUK/sum(PopUK));
y0(Indices + 1) = (10^(-4))*(PopUK/sum(PopUK));
%Set no vaccination
U = @(t) zeros(16,1);
%Solve ODE system
[t3,y3] = ode23(@(t,y) SIRODE(y,U(t),BetaUK,BetaUK*0.5,BetaUK*0.5,BetaUK*0.25,ones(16,1),ones(16,1),PopUK/sum(PopUK)),[0,100],y0);
%Create vector of final R values
Rvals2 = zeros(1,32);
Rvals2(1:16) = y3(length(y3(:,1)),Indices+2);
Rvals2(17:32) = zeros(1,16);
%Create new matrix with 2n groups
rho = 0.5;
chi = 0.5;
BetaDashUK = [BetaUK BetaUK*rho; BetaUK*chi BetaUK*rho*chi];
%Solve for M
M = zeros(32,32);
SIndices = zeros(32,1);
SIndices(1:16) = Indices;
SIndices(17:32) = Indices + 3;
for i = 1:32
    for j = 1:32
        M(i,j) = (eq(i,j) - y0(SIndices(i)) * BetaDashUK(i,j) * exp(-dot(BetaDashUK(i,:),Rvals2)) ) / (1 - exp(-dot(BetaDashUK(i,:),Rvals2)));
    end
end

%Input case fatality rates
Unifp = ones(32,1)*0.1;
Unifp(17:32) = 0.01*ones(16,1);

CovidInit = [0.00299245340815
0.000598959227885
0.000938864114634
0.002535768420739
0.006230466251811
0.012372736200681
0.023995896845171
0.036734022963793
0.069577008713216
0.117082813689097
0.197024642405357
0.308834654916426
0.472778299511567
1.00791673834262
1.65644905016982
3.53138611819271
];
%Scale largest agegroup
Covidp = [CovidInit*0.01;CovidInit*0.001];
TotalOver80 =  1737	 + 1078	+   491	+   130	 +  16;
Covidp(16) = 0.01*8.88444160047717*(TotalOver80)/PopUK(16) + Covidp(15)*(PopUK(16) - TotalOver80)/PopUK(16);
Covidp(32) = 0.1*Covidp(16);

%Solve system to get linear approximation

x = M'\Unifp;

yUnif = zeros(16,1);
for i = 1:16
    yUnif(i) = (y0(Indices(i))*sum(PopUK)/PopUK(i))*(x(i+16) - x(i));
end

x = M'\Covidp;

yCovid = zeros(16,1);
for i = 1:16
    yCovid(i) = (y0(Indices(i))*sum(PopUK)/PopUK(i))*(x(i+16) - x(i));
end
%Create figure
close all
figure;
bar([yUnif/min(yUnif),yCovid/min(yCovid)])
xlabel('Age Group (years)','Interpreter','latex')
ylabel('Proportion of Optimal Effectiveness','Interpreter','latex')
legend('Uniform Case Fatality Ratio','COVID-19 Case Fatality Ratio','location','northoutside')
labels = categorical({'0-4','5-9','10-14','15-19','20-24','25-29','30-34','35-39','40-44','45-49','50-54','55-59','60-64','65-69','70-74','75+'});
labels = reordercats(labels,{'0-4','5-9','10-14','15-19','20-24','25-29','30-34','35-39','40-44','45-49','50-54','55-59','60-64','65-69','70-74','75+'});
xticks(1:16)
xticklabels(labels)
set(gca,'fontsize',18)
print('-depsc', 'Fig6.eps')

%% Fig 7
%Select different values of epsilon
epsvals = [0.01,0.0498];
y0 = zeros(16*7,1);
%Initial conditions
y0(Indices) = (1-10^(-4))*(PopUK/sum(PopUK));
y0(Indices + 1) = (10^(-4))*(PopUK/sum(PopUK));
ObjChanges = zeros(16,2);
%Simulate epidemic for COVID death rates
for i = 1:16
    for eps = 1:2
        epsilon = epsvals(eps);
        [t3,y3] = ode23(@(t,y) SIRODE(y,U_Disc(t,epsilon,1,i),BetaUK,BetaUK*0.5,BetaUK*0.5,BetaUK*0.25,ones(16,1),ones(16,1),PopUK/sum(PopUK)),[0,100],y0);
        RvalsnewU = y3(length(y3(:,1)),Indices + 2);
        RvalsnewV = y3(length(y3(:,1)),Indices + 5);
        ObjChanges(i,eps) = dot(RvalsnewU,Covidp(1:16)) + dot(RvalsnewV,Covidp(1:16))*0.1 - dot(Rvals2(1:16),Covidp(1:16));
        disp(eps)
    end
    disp(i)
    disp('----------------')
end
%Simulate epidemic for uniform death rates
ObjChanges2 = zeros(16,2);
for i = 1:16
    for eps = 1:2
        epsilon = epsvals(eps);
        [t3,y3] = ode23(@(t,y) SIRODE(y,U_Disc(t,epsilon,1,i),BetaUK,BetaUK*0.5,BetaUK*0.5,BetaUK*0.25,ones(16,1),ones(16,1),PopUK/sum(PopUK)),[0,1000],y0);
        RvalsnewU = y3(length(y3(:,1)),Indices + 2);
        RvalsnewV = y3(length(y3(:,1)),Indices + 5);
        ObjChanges2(i,eps) = dot(RvalsnewU,Unifp(1:16)) + dot(RvalsnewV,Unifp(1:16))*0.1 - dot(Rvals2(1:16),Unifp(1:16));
        disp(eps)
    end
    disp(i)
    disp('----------------')
end

%Create figure (note extra formatting steps)
close all


hold on
box on
bar([-ObjChanges(:,2),-yCovid*0.0498],'w','EdgeColor','none');

print('-depsc','test2.eps')
xlabel('Vaccinated Age Group','Interpreter','latex')
ylabel('Objective Decrease','Interpreter','latex')
legend('Predicted Values','Actual Values','location','northwest')
title('COVID-19 Case Fatality Ratio')
labels = categorical({'0-4','5-9','10-14','15-19','20-24','25-29','30-34','35-39','40-44','45-49','50-54','55-59','60-64','65-69','70-74','75+'});
labels = reordercats(labels,{'0-4','5-9','10-14','15-19','20-24','25-29','30-34','35-39','40-44','45-49','50-54','55-59','60-64','65-69','70-74','75+'});
xticks(1:16)
xticklabels(labels)
set(gca,'fontsize',18)
%print('-depsc', 'ObjCompCovid.eps')
figure;
tiledlayout(2,1)
nexttile
bar([-ObjChanges(:,2),-yCovid*0.0498],'w','EdgeColor','none');
%print('-depsc','testhere.eps')
b = bar([-ObjChanges(:,2),-yCovid*0.0498],'facecolor','flat');
b(1).CData = [ 0    0.4470    0.7410].*ones(16,3);
b(2).CData = [0.8500    0.3250    0.0980];
xlabel('Vaccinated Age Group','Interpreter','latex')
ylabel('Objective Decrease','Interpreter','latex')
legend([b(1),b(2)],["Predicted Values","Actual Values"],'location','northwest')
title('(a): COVID-19 Case Fatality Ratio')
print('-depsc', 'ObjCompCovid.eps')
labels = categorical({'0-4','5-9','10-14','15-19','20-24','25-29','30-34','35-39','40-44','45-49','50-54','55-59','60-64','65-69','70-74','75+'});
labels = reordercats(labels,{'0-4','5-9','10-14','15-19','20-24','25-29','30-34','35-39','40-44','45-49','50-54','55-59','60-64','65-69','70-74','75+'});
xticks(1:16)
xticklabels(labels)
nexttile
bar([-ObjChanges2(:,2),-yUnif*0.0498],'w','EdgeColor','none');
%print('-depsc','testhere.eps')
b = bar([-ObjChanges2(:,2),-yUnif*0.0498],'facecolor','flat');
b(1).CData = [ 0    0.4470    0.7410].*ones(16,3);
b(2).CData = [0.8500    0.3250    0.0980];
xlabel('Vaccinated Age Group','Interpreter','latex')
ylabel('Objective Decrease','Interpreter','latex')
legend([b(1),b(2)],["Predicted Values","Actual Values"],'location','northeast')
title('(b): Uniform Case Fatality Ratio')
labels = categorical({'0-4','5-9','10-14','15-19','20-24','25-29','30-34','35-39','40-44','45-49','50-54','55-59','60-64','65-69','70-74','75+'});
labels = reordercats(labels,{'0-4','5-9','10-14','15-19','20-24','25-29','30-34','35-39','40-44','45-49','50-54','55-59','60-64','65-69','70-74','75+'});
xticks(1:16)
xticklabels(labels)
print('-depsc', 'Fig7.eps')

