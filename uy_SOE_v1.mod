
var N x xf URlog u labstar ls  z zf lam lamf robs pinfobs dy dc dinve dwAWE  ewma epinfma  zcapf rkf kf pkf    cf invef yf labf wf rrf mc zcap rk k pk    c inve yH lab pinf w r a  b g qs  ms  spinf sw kpf kp winfAWE winf ygap;    
var rstar xi xCF xCH pH pF xIF xIH pI xHstar rer y_star bstar piS zeta zeta2 tb m pistar dxHstar Rstar_obs piS_obs pistar_obs xHstarf rstarf pinfh piHobs dy_star xi_obs; 
var auxi;
varexo ea eb eg  eqs  em  epinf ew  els eAWE eps_zeta eps_zeta2 eps_y_star eps_pistar eps_Rstar eps_auxi;  
 
parameters cla cwtrend cchi cgy curvp constelab constepinf constebeta cmaw cmap calfa 
czcap  csadjcost ctou csigma chabb  cfc 
cindw cprobw cindp cprobp csigl clandaw 
crpi crdy cry crr 
crhoa crhob crhog crhoqs crhoms crhopinf crhow  
ctrend cg  ETA_C O_C ETA_I O_I ETA_STAR PSI1 PSI2 RHO_zeta RHO_zeta2 RHO_y_star RHO_pistar RHO_Rstar
constepistar constepiS pHss cstb constepiH constexi;


// fixed parameters
ctou=0.06/4;  //tasa de depreciación del capital -Tasa anual de depreciacion de 6% (de acuerdo al MMET) 
clandaw=1.1; //mark up de estado estacionario sobre salarios - 10% como en MMET
cg=0.111308403056381; // gasto exogeno per capita/pib per capita promedio 2005-2018
curvp=1; //kimball precios = 1 para que sea Dixit Stiglitz


// estimated parameters initialisation
calfa=1-0.7;  // participacion del capital en la produccion dsge banco
cgamma_init=1.00816735542229; // Crecimiento del PIB per capita de lp. Promedio 2005-2018
cbeta=0.99999999999;  //factor de descuento del consumo - consistente con R del banco, pie y gama , me da más que 1
csigma=1.0;  //elasticidad de sustitución del consumo dsge banco
cpie_init=1.01824098730178;  //inflacion de estado estacionario - promedio 2005-2018 del IPI trim
cfc=1.1;  //costos fijos en beneficios tales que el mark up de precios en ss sea 10%
cgy=0.51;  //parametro que mide cuanto afecta el shock a la productividad, via exportaciones, la evolucion de la productividad domestica - Igual al paper

csadjcost= 4.1*2;  //elasticidad de estado estacionario de la función de costos de ajuste del capital (elasticidad de la inversión a un cambio en q)
chabb=    0.7285;    //habitos - dsge banco
cprobw=    0.8203;  //calvo w - dsge banco
csigl=    1;  //Frisch elasticity - dsge banco
cprobp=   0.8982;  //calvo w -dsge banco
cindw=    0.5116; //indexacion de salarios - dsge banco
cindp=    0.4823; //indexacion de precios - dsge banco
czcap=    0.2696*2;  // funcion creciente en la elasticidad de la función de costos de ajuste del capital. Psi. The elasticity of the capital
//utilisation cost function has a mean of 0.2, and includes in its domain the value of 0.1 suggested by King
//and Rebelo (2000).
crpi=     1.5384;  //regla infla - DSGE banco
crr=      0.8324;  //regla tasa -DSGE banco
cry=      0.1230;  //regla prod -DSGE banco
crdy=     0.1230;   //regla crecimiento 
ETA_STAR=1;//banco
ETA_I=0.9;//banco es 1.5 pero se rompe calibracion del ratio consumo a producto. 
ETA_C=1.2;//banco
O_I=0.425;//banco
O_C=0.64246;//banco
PSI1=0.0101; //banco
PSI2=0;//banco
crhoa=    0.3642;  //persistencia de shock a la productividad del banco
crhob=    0.3832; //persistencia de shock a la prima de riesgo del banco del shock a las preferencias
crhog=   0.818213; //persistencia de shock al gasto exógeno, para usa daba 0.9957 mucho mas alto. Habra un tema con el detrended?
RHO_Rstar=0.980857;       //Antes 0.983763
RHO_y_star=0.866212;       //Antes 0.87635;
RHO_pistar=0.252705;      //Antes 0.247438; 
RHO_zeta=0.75;//0.7826;
RHO_zeta2=0.75;//0.7826;
crhoqs=   0.2834;  //persistencia precio del capital
crhoms=0; //persistencia del shock de politica monetaria
crhopinf=0; //persistencia del shock al markup de precios
crhow=0; //persistncia mark up de salarios
cmap = 0; //componente ma del shock a los marks ups de precios
cmaw  = 0; //componente ma del shock a los marks ups de salarios
cwtrend=0.0069675267425529;

cxi_init=1.00601511;
cpiS_init=1.00459364;
constebeta=0.000000001;
ctrend=(cgamma_init-1)*100;
constepiS=(cpiS_init-1)*100;  //Promedio 2005.Q2-2018.Q4
constepistar = (cpie_init/cpiS_init-1)*100;//  
constexi     =  (cxi_init-1)*100;     // Promedio 2005.Q2-2018.Q4 pensar #
constepinf=(cpie_init-1)*100;
constepiH=(cpie_init-1)*100;
constelab=0;
cchi = 0.1;
chntheta = 1;
cltheta = 0;
cla = 0;
pHss=1;
cstb=0;


// derived from steady state
clandap=cfc;//
cbetabar=cbeta*cgamma_init^(-csigma); //beta techo en estado estacionario está hallada
cr_init=cpie_init/(cbeta*cgamma_init^(-csigma)); //en estado estacionario, la tasa de interes nominal bruta es igual a la infla de ss sobre beta/gamma la halle y sale de euler
crstar_init=cpie_init/(cxi_init*cpiS_init*cbeta*cgamma_init^(-csigma));
crk=(cbeta^(-1))*(cgamma_init^csigma) - (1-ctou); //Sale de la eq del precio del capital
cw = (calfa^calfa*(1-calfa)^(1-calfa)/(clandap*crk^calfa))^(1/(1-calfa));//ok, sale de considerar que en estado estacionario, los costos marginales son 1/markup que es clandap
cikbar=(1-(1-ctou)/cgamma_init); //Sale de eq de mov el k. kbar en este caso es el cap total
cik=(1-(1-ctou)/cgamma_init)*cgamma_init; //Sale de la eq de mov del k pero en vez de en función del k total en función del k utilizado kut=kbar/gamma en ss
clk=((1-calfa)/calfa)*(crk/cw);//condiciones de primer orden de las firmas de  bienes intermedios
cky=cfc*(clk)^(calfa-1); //sale de definición de cfc, despejo y y hago k/ese y. 
ciy=cik*cky;//prop inversión en k por prop de k en y
pFss = ((1/O_C)-((1-O_C)/O_C)*pHss^(1-ETA_C))^(1/(1-ETA_C));
pIss = ((1-O_I)*(pHss)^(1-ETA_I)-O_I*(pFss)^(1-ETA_I))^(1/(1-ETA_I));
ccy=1-cg-cik*cky*(pIss/pHss)-cstb;
cxchy = (1-O_C)*(pHss)^(-ETA_C)*ccy;
cxihy = (1-O_I)*(pHss/pIss)^(-ETA_I)*ciy;
cxHstary = 1-cxchy-cxihy-cg;
cxcfy =O_C*(pFss)^(-ETA_C)*ccy;
cxify = O_I*(pFss/pIss)^(-ETA_I)*ciy;
cmy = cxcfy+cxify;
cRER=pFss/pHss;
crkky=crk*cky; //RENTAS DEL CAPITAL TOTAL? 
cwhlc=(1/clandaw)*(1-calfa)/calfa*crk*cky/ccy;
clab=(cwhlc*(1/(1-chabb/cgamma_init))*cgamma_init^((1-cchi)/cchi))^(1/(1+csigl));
cy=(1/cfc)*clab*(1/clk)^(calfa);
cwly=1-crk*cky;

conster=(cr_init-1)*100;
consterstar=(crstar_init-1)*100;// 

model(linear); 

//#usmodel_stst;
# cpie=1+constepinf/100;
# cpistar = 1+constepistar/100;
# cdep=1+constepiS/100;
# cgamma=1+ctrend/100 ;
# cbeta=1/(1+constebeta/100);
# cxi = 1+constexi/100;
# clandap=cfc;
# cbetabar=cbeta*cgamma^(-csigma);
# cr=cpie/(cbeta*cgamma^(-csigma));
# cRstar = cpie/(cxi*cdep*cbeta*cgamma^(-csigma));
# crk=(cbeta^(-1))*(cgamma^csigma) - (1-ctou);
# cw = (calfa^calfa*(1-calfa)^(1-calfa)/(clandap*crk^calfa))^(1/(1-calfa));
# cikbar=(1-(1-ctou)/cgamma);
# cik=(1-(1-ctou)/cgamma)*cgamma;
# clk=((1-calfa)/calfa)*(crk/cw);
# cky=cfc*(clk)^(calfa-1);
# ciy=cik*cky;
# pFss= ((1/O_C)-((1-O_C)/O_C)*pHss^(1-ETA_C))^(1/(1-ETA_C));
# pIss= ((1-O_I)*(pHss)^(1-ETA_I)-O_I*(pFss)^(1-ETA_I))^(1/(1-ETA_I));
# ccy=1-cg-cik*cky*(pIss/pHss)-cstb;
# cxchy = (1-O_C)*(pHss)^(-ETA_C)*ccy;
# cxihy = (1-O_I)*(pHss/pIss)^(-ETA_I)*ciy;
# cxHstary = 1-cxchy-cxihy-cg;
# cxcfy =O_C*(pFss)^(-ETA_C)*ccy;
# cxify = O_I*(pFss/pIss)^(-ETA_I)*ciy;
# cmy= cxcfy+cxify;
# cRER= pFss/pHss;
# crkky=crk*cky;
# cwhlc=(1/clandaw)*(1-calfa)/calfa*crk*cky/ccy;
# clab=(cwhlc*(1/(1-chabb/cgamma))*cgamma^((1-cchi)/cchi))^(1/(1+csigl));
# cy=(1/cfc)*clab*(1/clk)^(calfa);
# cwly=1-crk*cky;
# conster=(cr-1)*100;
# consterstar=(cRstar-1)*100;
//# consteu = (1-(1/clandaw)^(1/csigl));
# consteu = (clandaw-1)/csigl;

// flexible economy

	      0 =  calfa*rkf+(1-calfa)*(wf) - a ; //(1)
	      zcapf =  (1/(czcap/(1-czcap)))* rkf  ; //(2)
	      rkf =  wf+labf-kf ; //(3)
	      kf =  kpf(-1)+zcapf ; //(4)
	      invef = (1/(1+cbetabar*cgamma))* (  invef(-1) + cbetabar*cgamma*invef(1)+(1/(cgamma^2*csadjcost))*pkf ) +qs ; //(5)
          pkf = -rrf + b  +(crk/(crk+(1-ctou)))*rkf(1) +  ((1-ctou)/(crk+(1-ctou)))*pkf(1) ; //(6)
	      lamf = lamf(+1) + (rrf-b); //(7)
          lamf = lamf(+1) + (rstarf+xi) ; //(8)
          lamf = -csigma/(1-chabb/cgamma)*cf + csigma*(chabb/cgamma)/(1-chabb/cgamma)*cf(-1) ; //(9)
	      yf = cxchy*cf+cxihy*invef+g  +  1*crkky*zcapf+cxHstary*xHstarf ; //(10)
          yf = cfc*( calfa*kf+(1-calfa)*labf +a ); //(11)
 	      wf = csigl*labf  	- lamf + ls + xf; //(12)
	      kpf =  (1-cikbar)*kpf(-1)+(cikbar)*invef + (cikbar)*(cgamma^2*csadjcost)*qs ; //(13)
          xf = zf-1/(1-chabb/cgamma)*cf + (chabb/cgamma)/(1-chabb/cgamma)*cf(-1); //shock endogeno de preferencias //(14)
          zf = (1-cchi)*zf(-1) + cchi/(1-chabb/cgamma)*cf - cchi*(chabb/cgamma)/(1-chabb/cgamma)*cf(-1); //tendencia consumo //(15) 
          xHstarf = y_star; //(16)
// sticky price - wage economy

	      mc =  calfa*rk+(1-calfa)*(w) - 1*a ; //(17)
	      zcap =  (1/(czcap/(1-czcap)))* rk ; //(18)
	      rk =  w+lab-k ; //(19)
	      k =  kp(-1)+zcap ; //(20)
	      inve = (1/(1+cbetabar*cgamma))* (  inve(-1) + cbetabar*cgamma*inve(1)+(1/(cgamma^2*csadjcost))*pk ) +qs ; //(21)
          pk = -r+pinf(1)+b + (crk/(crk+(1-ctou)))*rk(1) +  ((1-ctou)/(crk+(1-ctou)))*pk(1) ; //(22)
	      lam = lam(+1) + (r-pinf(+1)-b) ; //(23)
          lam = -csigma/(1-chabb/cgamma)*c + csigma*(chabb/cgamma)/(1-chabb/cgamma)*c(-1); //(24)
	      yH = cxchy*c+cxihy*inve+g  +  1*crkky*zcap+cxHstary*xHstar ; //(25)
	      yH = cfc*( calfa*k+(1-calfa)*lab +a ); //(26)
	      pinf =  (1/(1+cbetabar*cgamma*cindp)) * ( cbetabar*cgamma*pinf(1) +cindp*pinf(-1) 
               +((1-cprobp)*(1-cbetabar*cgamma*cprobp)/cprobp)/((cfc-1)*curvp+1)*(mc+100*spinf)  )  + 0*spinf ; //(27)
	      w - w(-1) + pinf =
               +(cbetabar*cgamma)*(w(1)-w + pinf(1))
               +(cindw)*pinf(-1)
               -(cbetabar*cgamma*cindw)*pinf
               +(1-cprobw)*(1-cbetabar*cgamma*cprobw)/(cprobw*(1+(clandaw/(clandaw-1))*csigl))
                    *(-csigl*u + 100*sw) ; //(28)
	      r =  crpi*(1-crr)*pinf
               +cry*(1-crr)*(yH-yf)     
               +crdy*(yH-yf-yH(-1)+yf(-1))
               +crr*r(-1)
               +ms  ; //(29)
	      kp =  (1-cikbar)*kp(-1)+cikbar*inve + cikbar*cgamma^2*csadjcost*qs ; //(30)
          x = z - 1/(1-chabb/cgamma)*c + (chabb/cgamma)/(1-chabb/cgamma)*c(-1); //(31)
          z = (1-cchi)*z(-1) + cchi/(1-chabb/cgamma)*c - cchi*(chabb/cgamma)/(1-chabb/cgamma)*c(-1);  //(32)
// economia abierta
    //hogares
          lam = lam(+1) + (rstar+xi+piS(+1)-pinf(+1)) ; //(33)
    //bienes finales consumo
           c   = (1/((1-O_C)^(1/ETA_C)*cxchy^((ETA_C-1)/ETA_C)+O_C^(1/ETA_C)*cxcfy^((ETA_C-1)/ETA_C))^(ETA_C/(ETA_C-1)))*((xCH*(1-O_C)^(1/ETA_C)*cxchy^((ETA_C-1)/ETA_C)+xCF*O_C^(1/ETA_C)*cxcfy^((ETA_C-1)/ETA_C))); // //(34)
           xCF = -ETA_C*(pF)+c; //(35)
           xCH = -ETA_C*(pH)+c; //(36)
    //bienes finales inversion
           inve   = (1/((1-O_I)^(1/ETA_I)*cxihy^((ETA_I-1)/ETA_I)+O_I^(1/ETA_I)*cxify^((ETA_I-1)/ETA_I))^(ETA_I/(ETA_I-1)))*((xIH*(1-O_I)^(1/ETA_I)*cxihy^((ETA_I-1)/ETA_I)+xIF*O_I^(1/ETA_I)*cxify^((ETA_I-1)/ETA_I))); //  (37)
           xIF = -ETA_I*(pF-pI)+inve; //(38)
           xIH = -ETA_I*(pH-pI)+inve; //(39)
           
           
   //exportaciones y riesgo pais         
           xHstar =-ETA_STAR*pH+ETA_STAR*rer+y_star; // (40)
           xi = -PSI1*cRER*bstar-PSI2*(piS(+1)+piS)+zeta+zeta2; // (41)
   //agregacion y equilibrio de mercado 
           tb                    = pHss*cxHstary*cy*(pH+xHstar)-cRER*cmy*cy*(rer+m); // (42) pensar no tengo los coeficientes del ss 
           bstar                 = (cRstar*cxi/cpistar)*bstar(-1)+(1/cRER)*tb; // (43)
           m                     = (cxcfy/cmy)*xCF+(cxify/cmy)*xIF; // (44)
           yH = ccy*c+ciy*inve+g  +  1*crkky*zcap+cxHstary*xHstar-cmy*m; // (45)
   //por definir
           pH = pH(-1) + pinfh - pinf; // (46)
           //pF-pF(-1) = piF-pinf; // 
           //exp(pY)*exp(y)        = exp(c)+exp(pI)*exp(invest)+exp(pH)*exp(g)+tb; //
           


// exogenous disturbances

	      a = crhoa*a(-1)  + ea; //(47)
	      b = crhob*b(-1) + eb; //(48)
	      g = crhog*(g(-1)) + eg + cgy*ea; //(49)
	      qs = crhoqs*qs(-1) + eqs; //(50)
	      ms = crhoms*ms(-1) + em;  //(51)
	      spinf = crhopinf*spinf(-1) + epinfma - cmap*epinfma(-1); //(52)
	          epinfma=epinf; //(53)
	      sw = crhow*sw(-1) + ewma - cmaw*ewma(-1) ; //(54)
	          ewma=ew; //(55)
          zeta   = RHO_zeta*zeta(-1)+eps_zeta;             // (E.56) Country premium 
          zeta2  = RHO_zeta2*zeta2(-1)+eps_zeta2;             // (E.57) Country premium
          y_star = RHO_y_star*y_star(-1)+eps_y_star;     // (E.58) Rest of the world output
          pistar  = RHO_pistar*pistar(-1)+eps_pistar;     // (E.59) Rest of the world inflation
          rstar   = RHO_Rstar*rstar(-1)+eps_Rstar;         // (E.60) Rest of the world interest rate

// labor supply

u = labstar - lab; //(61)
ls = (0.999)*ls(-1) + els + cla*ea ; //(62)
w = csigl*labstar - lam + ls + x; //(63)

// measurement equations

dy=yH-yH(-1)+ctrend; //(64)
dy_star=y_star-y_star(-1)+ctrend; //(65)
dc=c-c(-1)+ctrend;  //(66)
dinve=inve-inve(-1)+ctrend; //(67)
dxHstar=xHstar-xHstar(-1)+ctrend; //(68)
dwAWE=w-w(-1)+cwtrend+eAWE; //(69)

pinfobs = 1*(pinf) + constepinf; //(70)
robs =    1*(r) + conster; //(71)
N = lab + constelab; //(72)
URlog = u + 100*consteu; //(73)
Rstar_obs    = 1*(rstar) + consterstar     ; //(74)
piS_obs      = piS + constepiS  ;   //(75)
pistar_obs   = pistar +constepistar     ; //(76)
piHobs = pinfh + constepiH; //(77)
xi_obs       = xi + constexi ;  % + me_xi_obs ; //(78)

winfAWE = w-w(-1)+pinf+eAWE; //(79)
winf = w-w(-1)+pinf; //(80)
rer= rer(-1)+ piS + pistar - pinf; //(81)
ygap = yH-yf; //(82)

auxi = eps_auxi;

end; 

steady_state_model;
dy=ctrend;
dy_star=ctrend;
dc=ctrend;
dinve=ctrend;
dxHstar=ctrend;
dw=0;
dwAWE=cwtrend;
URlog=100*(clandaw-1)/csigl;
pinfobs = constepinf;
piHobs = constepiH;
robs = (((1+constepinf/100)/((1/(1+constebeta/100))*(1+ctrend/100)^(-csigma)))-1)*100;
N = constelab;
Rstar_obs=  ((1+robs/100)*((1+constepiS/100)^(-1))*((1+constexi/100)^(-1))-1)*100;
rer_obs=((1/O_C)-((1-O_C)/O_C)*pHss^(1-ETA_C))^(1/(1-ETA_C));//
piS_obs=constepiS;
pistar_obs   = constepistar;
xi_obs       =  constexi ; 
end;



shocks;
var ea;
stderr 0.0467;
var eb;
stderr 0.0638;
var eg;
stderr 0.6090;
var eqs;
stderr 0.2140;
var em;
stderr 0.1954;
var epinf;
stderr 3.7549;
var ew;
stderr 0.2089;
var eps_zeta;
stderr 0.0162;
var eps_zeta2;
stderr 0.0601;
var eps_y_star;
stderr 0.009466;
var eps_pistar;
stderr 0.037808;
var eps_Rstar;
stderr 0.001021;
var eps_auxi; 
stderr .00001;
end;

steady;
check;


estimated_params;
// PARAM NAME, INITVAL, LB, UB, PRIOR_SHAPE, PRIOR_P1, PRIOR_P2, PRIOR_P3, PRIOR_P4, JSCALE
// PRIOR_SHAPE: BETA_PDF, GAMMA_PDF, NORMAL_PDF, INV_GAMMA_PDF
//stderr ea,0.48   ,0.01,10,UNIFORM_PDF,,,0,5;
stderr ea, inv_gamma_pdf, 0.1, inf;
//stderr eb,1.8  ,0.01,10,UNIFORM_PDF,,,0,5;
stderr eb, inv_gamma_pdf, 0.1, inf;
//stderr eg,0.49    ,0.01,10,UNIFORM_PDF,,,0,5;
stderr eg, inv_gamma_pdf, 0.1, inf;
//stderr eqs,0.42   ,0.01,10,UNIFORM_PDF,,,0,5;
stderr eqs, inv_gamma_pdf, 0.1, inf;
//stderr em,0.22   ,0.01,10,UNIFORM_PDF,,,0,5;
stderr em, inv_gamma_pdf, 0.1, inf;
//stderr epinf,0.06,0.01,10,UNIFORM_PDF,,,0,5;
stderr epinf, inv_gamma_pdf, 0.1, inf;
stderr ew, inv_gamma_pdf, 0.1, inf;
//stderr els,1.0 ,0.01,10,UNIFORM_PDF,,,0,5;
stderr els, inv_gamma_pdf, 0.1, inf;
//stderr eAWE,0.35  ,0.01,10,UNIFORM_PDF,,,0,5;
stderr eAWE, inv_gamma_pdf, 0.1, inf;
stderr eps_zeta2, inv_gamma_pdf, 0.1, inf;
stderr eps_zeta, inv_gamma_pdf, 0.1, inf;
//corr ea,els, 0.0,-0.99,0.99,UNIFORM_PDF,,,-1,1;
 crhoa, beta_pdf, 0.75/2, 0.1;
//crhoa,.975 ,.01,.9999,BETA_PDF,0.5,0.20;
crhob,.25  ,.01,.9999,BETA_PDF,0.5,0.20;
crhog,.975 ,.01,.9999,BETA_PDF,0.5,0.20;
crhoqs, beta_pdf, 0.70, 0.15;
//crhoqs,.7 ,.01,.9999,BETA_PDF,0.5,0.20;
crhoms,.06 ,.01,.9999,BETA_PDF,0.5,0.20;
crhopinf,.88 ,.01,.9999,BETA_PDF,0.5,0.2;
crhow,.975   ,.001,.9999,BETA_PDF,0.5,0.2;
cmap,.73     ,0.01,.9999,BETA_PDF,0.5,0.2;
cmaw,.53     ,0.01,.9999,BETA_PDF,0.5,0.2;
csadjcost,4 ,0.1,15,NORMAL_PDF,4,1.0;
// csigma,1.6    ,0.25,10,NORMAL_PDF,1.50,0.375;
chabb,0.55    ,0.001,0.99,BETA_PDF,0.7,0.1;
csigl,3.40     ,0.25,15,NORMAL_PDF,2,1.0;
cprobw, beta_pdf, 0.75, 0.1;
//cprobw,0.65   ,0.1,0.95,BETA_PDF,0.5,0.15;
//cprobp,0.66  ,0.1,0.95,BETA_PDF,0.5,0.15;
cprobp, beta_pdf, 0.75, 0.1;
cindp, beta_pdf, 0.5, 0.1;
cindw, beta_pdf, 0.5, 0.1;
//cindw,0.13    ,0.01,0.99,BETA_PDF,0.5,0.15;
//cindp,0.25    ,0.01,0.99,BETA_PDF,0.5,0.15;
czcap,0.75     ,0.01,1,BETA_PDF,0.5,0.15;
cfc,1.5      ,1.0,3,NORMAL_PDF,1.25,0.125;
crpi,1.8     ,1.0,3,NORMAL_PDF,1.5,0.25;
crr,0.86      ,0.01,0.975,BETA_PDF,0.75,0.10;
cry,0.13  ,0.001,0.5,NORMAL_PDF,0.125,0.05;
crdy,0.22 ,0.001,0.5,NORMAL_PDF,0.125,0.05;
//constepinf,1.01  ,0.01,2.0,GAMMA_PDF,0.625,0.1;//20;
//constebeta,0.0000000001  ,0.0000000001,2.0,GAMMA_PDF,0.25,0.1;//0.20;
constelab,0.0   ,-10.0,10.0,NORMAL_PDF,0.0,2.0;
//ctrend,0.32  ,0.01,0.8,NORMAL_PDF,0.4,0.10;
cwtrend,0.10  ,0.01,0.8,NORMAL_PDF,0.2,0.10;
cgy,0.55      ,0.01,2.0,NORMAL_PDF,0.5,0.25;
//calfa,0.18   ,0.01,1.0,NORMAL_PDF,0.3,0.05;
clandaw,1.25 ,1.01,10.0,NORMAL_PDF,1.5,0.250;
cchi,0.016 ,0.001,0.999,BETA_PDF,0.5,0.20;
// cla,0.8 ,0.001,0.999,BETA_PDF,0.5,0.20;

//ETA_STAR, , 0, , normal_pdf, 0.25, 0.2;  
ETA_C, , 0, , normal_pdf, 1.2, 0.3; 
ETA_I, , 0, , normal_pdf, 0.8, 0.15;
PSI1, inv_gamma_pdf, 0.005, inf; //
//RHO_zeta, beta_pdf, 0.75, 0.1;
//RHO_zeta2, beta_pdf, 0.75, 0.1;
end;


estimated_params_init(use_calibration);
end;


varobs dy dy_star dc dinve dxHstar N pinfobs robs URlog dwAWE Rstar_obs piS_obs pistar_obs piHobs xi_obs;
//estimation(optim=('MaxIter',200),datafile=SOE_data,mode_compute=5,first_obs=1,nobs=55,presample=4,lik_init=2,prefilter=0,mh_replic=0,mh_nblocks=2,mh_jscale=0.20,mh_drop=0.2,plot_priors=0);
estimation(optim=('MaxIter',200),datafile=SOE_data,mode_compute=5,mh_replic=0,plot_priors=1,mode_check);

//estimation(datafile=DATOS_DSGE_abr16, xls_sheet=Datos_transf_demean, plot_priors=1, mode_check, mode_compute=4, mh_replic=0);

//stoch_simul(irf=20) dy yH c inve lab labstar URlog pinfobs w robs ygap winf dwAWE;
//stoch_simul(irf=20) dy yH c inve lab pinfobs dwAWE w robs ;
stoch_simul(order=1, periods=0, irf=40, nograph);
//corrcoef([oo_.SmoothedShocks.ea oo_.SmoothedShocks.eb oo_.SmoothedShocks.eg oo_.SmoothedShocks.eqs oo_.SmoothedShocks.em oo_.SmoothedShocks.epinf oo_.SmoothedShocks.ew oo_.SmoothedShocks.els oo_.SmoothedShocks.eAWE ])

 