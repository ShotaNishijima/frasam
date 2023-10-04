// State space assessment model from Nielsen and Berg 2014, Fisheries Research.
//  --------------------------------------------------------------------------
// Copyright (c) 2014, Anders Nielsen <an@aqua.dtu.dk>,
// Casper Berg <cbe@aqua.dtu.dk>, and Kasper Kristensen <kkr@aqua.dtu.dk>.
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//   * Redistributions of source code must retain the above copyright
//     notice, this list of conditions and the following disclaimer.
//   * Redistributions in binary form must reproduce the above copyright
//     notice, this list of conditions and the following disclaimer in the
//     documentation and/or other materials provided with the distribution.
//   * Neither the name of the assessment tool SAM nor the
//     names of its contributors may be used to endorse or promote products
//     derived from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL ANDERS NIELSEN, CASPER BERG OR KASPER
// KRISTENSEN BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS;
// OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
// WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR
// OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
// ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//  --------------------------------------------------------------------------

#include <TMB.hpp>
#include <iostream>


/* Parameter transform */
template <class Type>
Type f(Type x){return Type(2)/(Type(1) + exp(-Type(2) * x)) - Type(1);}

template <class Type>
Type square(Type x){return x*x;}

// sqrt
template<class Type>
Type sqrt(Type x){
  return pow(x,Type(0.5));
}

template<class Type>
Type objective_function<Type>::operator() ()
{
  DATA_VECTOR(fleetTypes);
  DATA_VECTOR(sampleTimes);
  // DATA_VECTOR(years);
  DATA_INTEGER(nobs);
  // DATA_VECTOR(idx1);
  // DATA_VECTOR(idx2);
  DATA_ARRAY(obs);
  DATA_ARRAY(propMat);
  DATA_ARRAY(stockMeanWeight);
  DATA_ARRAY(catchMeanWeight);
  DATA_ARRAY(natMor);
  DATA_ARRAY(landFrac);
  DATA_ARRAY(disMeanWeight);
  DATA_ARRAY(landMeanWeight);
  DATA_ARRAY(propF);
  DATA_ARRAY(propM);
  DATA_INTEGER(minAge);
  DATA_INTEGER(maxAgePlusGroup);
  DATA_INTEGER(rhoMode);
  DATA_IARRAY(keyLogFsta);
  DATA_ARRAY(keyLogQ);
  DATA_ARRAY(keyLogB);
  DATA_ARRAY(keyVarF);
  DATA_ARRAY(keyVarLogN);
  DATA_ARRAY(keyVarObs);
  DATA_INTEGER(stockRecruitmentModelCode);
  DATA_SCALAR(scale);
  DATA_SCALAR(gamma);
  // DATA_VECTOR(fbarRange);
  DATA_INTEGER(sel_def); // selectivity difenition: devided by maxF(0), meanF(1), maxage(2)
  DATA_INTEGER(b_random); // if 1, b estimated by random effects

  PARAMETER_VECTOR(logQ);
  PARAMETER_VECTOR(logB);
  PARAMETER_VECTOR(logSdLogFsta);
  PARAMETER_VECTOR(logSdLogN);
  PARAMETER_VECTOR(logSdLogObs);
  PARAMETER(rec_loga);
  PARAMETER(rec_logb);
  PARAMETER(logit_rho);
  PARAMETER_ARRAY(U);
  PARAMETER(trans_phi1); //　#recruitment autocorrelation
  PARAMETER(logSD_b);

  PARAMETER_ARRAY(logwaa); // row: age, column: year (行列がstockMeanWeightと逆であることに注意！)
  PARAMETER_VECTOR(beta_w0);
  PARAMETER_VECTOR(alpha_w);
  PARAMETER_VECTOR(rho_w);
  PARAMETER(omicron); // log(SD) for process error in weight growth
  // PARAMETER(iota); // log(SD) for process error in maturity growth
  PARAMETER_VECTOR(logCV_w); // CV for observation error in weight (for age0 and older)

  DATA_IVECTOR(iy);
  DATA_INTEGER(nlogF);
  DATA_INTEGER(nlogN);
  DATA_ARRAY(propMat2);
  DATA_SCALAR(alpha);
  DATA_VECTOR(Fprocess_weight);
  DATA_SCALAR(F_RW_order); //0: first order, 1: second order
  DATA_SCALAR(lambda);
  DATA_SCALAR(lambda_Mesnil);
  DATA_IVECTOR(model_weight_maturity); // whether modeling weight and maturity
  DATA_SCALAR(scale_number);
  DATA_INTEGER(dist_wobs); //distribution for weight: 0: lognormal, 1: gamma
  DATA_ARRAY(weight_weight); // 体重データのモデリングの時に使用する重み（レトロの時に重みをゼロにするために使用）
  DATA_ARRAY(maturity_weight);  // 体重データのモデリングの時に使用する重み（レトロの時に重みをゼロにするために使用）

  array<Type> logF(nlogF,U.cols()); // logF (6 x 50 matrix)
  array<Type> logN(nlogN,U.cols()); // logN (7 x 50 matrix)
  array<Type> exp_logF(nlogF,U.cols()); // F (6 x 50 matrix)
  array<Type> exp_logN(nlogN,U.cols()); // N (7 x 50 matrix)

  array<Type> stockMeanWeight_true(stockMeanWeight.rows(),stockMeanWeight.cols()); // N (7 x 50 matrix)
  // stockMeanWeight_true.fill(1.0);
  // if(model_weight_maturity(0)<=0){
    for(int i=0;i<stockMeanWeight_true.rows();i++){
      for(int j=0;j<stockMeanWeight_true.cols();j++){
        stockMeanWeight_true(i,j)=stockMeanWeight(i,j);
      }
    }
  // }

  for(int i=0;i<nlogN;i++)
    for(int j=0;j<U.cols();j++){
      logN(i,j)=U(i,j);  // Uの1行からnlogNがlogN
      exp_logN(i,j)=exp(logN(i,j));
    }
  for(int i=0;i<nlogF;i++)
    for(int j=0;j<U.cols();j++){
      logF(i,j)=U(i+nlogN,j);  // UのnlogN+1からnlogN+nlogFがlogF
      exp_logF(i,j)=exp(logF(i,j));
    }

  int timeSteps=logF.dim[1]; // 年数
  int stateDimF=logF.dim[0]; // nlogF
  int stateDimN=logN.dim[0]; // nlogN
  //Type rho=f(logit_rho);
  vector<Type> sdLogFsta=exp(logSdLogFsta);  // Fのsd
  vector<Type> varLogN=exp(logSdLogN*Type(2.0)); //  Nの???散
  vector<Type> varLogObs=exp(logSdLogObs*Type(2.0)); // Obsの???散
  vector<Type> ssb(timeSteps); // SSB (年数??????
  vector<Type> ssn(timeSteps);
  vector<Type> logssb(timeSteps); //  logssb ???年数??????
  vector<Type> B_total(timeSteps);
  vector<Type> F_mean(timeSteps);
  vector<Type> Catch_biomass(timeSteps);
  vector<Type> Exploitation_rate(timeSteps);

  //First take care of F
  matrix<Type> fvar(stateDimF,stateDimF);  // Fの???散??????
  matrix<Type> fcor(stateDimF,stateDimF);  // Fの相関??????
  vector<Type> fsd(stateDimF);  // Fのsdのvector

  for(int i=0; i<stateDimF; ++i){
    for(int j=0; j<stateDimF; ++j){
      if(i!=j){if(rhoMode==0){fcor(i,j)=0.000001;
      }else{if(rhoMode==1){fcor(i,j)=0.999999;
      }else{
        if(rhoMode==2){fcor(i,j)=1/(1+exp(-logit_rho));
        }else{
          fcor(i,j)=pow(1/(1+exp(-logit_rho)),Type(abs(i-j)));
        }}}}else{
          fcor(i,j)=1.0;}  //  対???=1???非対???=rho
      }
    fsd(i)=sdLogFsta(CppAD::Integer(keyVarF(0,i)));  //
  }
  for(int i=0; i<stateDimF; ++i){
    for(int j=0; j<stateDimF; ++j){
      fvar(i,j)=fsd(i)*fsd(j)*fcor(i,j);  // var-covを定義
    }
  }
  using namespace density;  // 多変量正規分布を使う宣言
  MVNORM_t<Type> neg_log_densityF(fvar);  // var-cov matirx fvarを持つMVN
  Type ans=0;
  array<Type> logF_resid(stateDimF,timeSteps); //
  if (F_RW_order==0) {
    for(int i=1;i<timeSteps;i++){
      ans+=Fprocess_weight(i)*neg_log_densityF(logF.col(i)-logF.col(i-1)); // F-Process likelihood
      SIMULATE {
        logF.col(i) = logF.col(i-1) + neg_log_densityF.simulate();
      }
    }
  } else { // second order
    for(int i=2;i<timeSteps;i++){
      ans+=Fprocess_weight(i)*neg_log_densityF(logF.col(i)-(Type(1.0)+F_RW_order)*logF.col(i-1)+F_RW_order*logF.col(i-2)); // F-Process likelihood
      SIMULATE {
        logF.col(i) = (Type(1.0)+F_RW_order)*logF.col(i-1) - F_RW_order*logF.col(i-2) + neg_log_densityF.simulate();
      }
    }
  }

  for(int i=0;i<timeSteps;i++){ // calc ssb
    ssb(i)=0.0;
    ssn(i)=0.0;
    for(int j=0; j<stateDimN; ++j){
      ssb(i)+=exp(logN(j,i))*exp(-exp(logF((keyLogFsta(0,j)),i))*propF(i,j)-natMor(i,j)*propM(i,j))*propMat(i,j)*stockMeanWeight_true(i,j);  // ssbを
      ssn(i)+=exp(logN(j,i))*exp(-exp(logF((keyLogFsta(0,j)),i))*propF(i,j)-natMor(i,j)*propM(i,j))*propMat(i,j);
      // ssb(i)+=exp(logN(j,i))*propMat(i,j)*stockMeanWeight(i,j);  // ssbを??????
    }
    logssb(i)=log(ssb(i));  // log(ssb)
  }

  array<Type> saa(stateDimN,timeSteps); //selectivity at age
  vector<Type> Sel1(timeSteps);
  vector<Type> FY(stateDimN);
  for (int y=0;y<timeSteps;y++){
    for (int i=0;i<stateDimN;i++) {
      if(i<stateDimN-1){
        FY(i)=exp_logF(i,y);
      }else{
        FY(i)=alpha*exp_logF(i-1,y);
        }
    }
    if (sel_def==0) { //max
      Sel1(y)=max(FY);
    }else{
      if (sel_def==1) { //mean
        Sel1(y)=sum(FY)/stateDimN;
      }else{//maxage
        Sel1(y)=FY(stateDimN-1);
      }
    }
    for (int i=0;i<stateDimN;i++){
      saa(i,y)=FY(i)/Sel1(y);
    }
  }

  vector<Type> predN0(stateDimN);  // logNの予測値
  vector<Type> predN(stateDimN);  // logNの予測値
  vector<Type> recResid(timeSteps); //再生産関係からの残差

  int start_timeStep=1+minAge;
  if(stockRecruitmentModelCode==0){ // if RW
    start_timeStep=1;
    }
  //Now take care of N
  matrix<Type> nvar(stateDimN,stateDimN);  // logNのvcov
  for(int k=0; k<stateDimN; ++k){
    for(int j=0; j<stateDimN; ++j){
      if(k!=j){nvar(k,j)=0.0;}else{nvar(k,j)=varLogN(CppAD::Integer(keyVarLogN(0,k)));} // logNには相関なし varは加入とそれより上で異なる
    }
  }
  MVNORM_t<Type> neg_log_densityN(nvar);

  // //For initial step
  // matrix<Type> nvar0(stateDimN,stateDimN);  // logNのvcov
  // for(int k=0; k<stateDimN; ++k){
  //   for(int j=0; j<stateDimN; ++j){
  //     if(k!=j){nvar0(k,j)=0.0;}else{
  //       nvar0(k,j)=varLogN(CppAD::Integer(keyVarLogN(0,k)))/(1-pow(phi1,Type(2.0)));
  //       } // logNには相関なし varは加入とそれより上で異なる
  //   }
  // }
  // MVNORM_t<Type> neg_log_densityN0(nvar0);
  Type phi1 = (exp(trans_phi1)-Type(1.0))/(exp(trans_phi1)+Type(1.0));

  for(int i=start_timeStep;i<timeSteps;i++){
    if(stockRecruitmentModelCode==0){ // straight RW
      predN0(0)=logN(0,i-1);
    }else{
      if(stockRecruitmentModelCode==1){//ricker
        // predN(0)=rec_loga+log(ssb(i-1))-exp(rec_logb)*ssb(i-1);
        predN0(0)=rec_loga+log(ssb(i-minAge))-exp(rec_logb)*ssb(i-minAge);
      }else{
        if(stockRecruitmentModelCode==2){//BH
          // predN(0)=rec_loga+log(ssb(i-1))-log(1.0+exp(rec_logb)*ssb(i-1));
          predN0(0)=rec_loga+log(ssb(i-minAge))-log(Type(1.0)+exp(rec_logb)*ssb(i-minAge));
        }else{
          if(stockRecruitmentModelCode==3){ //HS
            predN0(0)=CppAD::CondExpLt(rec_logb,log(ssb(i-minAge)),rec_loga+rec_logb,rec_loga+log(ssb(i-minAge)));
            // vector<Type> rec_pred_HS(2);
            // rec_pred_HS(0)=rec_loga+rec_logb;
            // rec_pred_HS(1)=rec_loga+log(ssb(i-minAge));
            //
            // predN0=min(rec_pred_HS);
          } else {
            if(stockRecruitmentModelCode==4){ //Mesnil HS
              predN0(0)=ssb(i-minAge)+sqrt(square(exp(rec_logb))+square(gamma)/Type(4.0))-sqrt(square(ssb(i-minAge)-exp(rec_logb))+square(gamma)/Type(4.0));
              predN0(0)*=exp(rec_loga)/Type(2.0);
              predN0(0)=log(predN0(0));
            }else{
              if(stockRecruitmentModelCode==5){ //Constant R0
                predN0(0)=rec_loga;
              }else{
                if(stockRecruitmentModelCode==6){ //Proportional to SSB (no density-dependence)
                  predN0(0)=rec_loga+log(ssb(i-minAge));
                }else{
                  error("SR model code not recognized");
                }
              }
              }
          }
        }
      }
    }
    recResid(i)=logN(0,i)-predN0(0);
    if(i==0) {
      predN(0)=predN0(0);
    } else {
      predN(0)=predN0(0)+phi1*recResid(i-1);
    }

    for(int j=1; j<stateDimN; ++j){
      if (j<(stateDimN-1)) {
        predN(j)=logN(j-1,i-1)-exp(logF((keyLogFsta(0,j-1)),i-1))-natMor(i-1,j-1);  // population dynamics model
      }else{
        predN(j)=logN(j-1,i-1)-exp(logF((keyLogFsta(0,j-1)),i-1))-natMor(i-1,j-1);  // population dynamics model
      }
    }
    if(maxAgePlusGroup==1){
      predN(stateDimN-1)=log(exp(logN(stateDimN-2,i-1)-exp(logF((keyLogFsta(0,stateDimN-2)),i-1))-natMor(i-1,stateDimN-2))+
                             exp(logN(stateDimN-1,i-1)-alpha*exp(logF((keyLogFsta(0,stateDimN-1)),i-1))-natMor(i-1,stateDimN-1))); // plus group
    }
    ans+=neg_log_densityN(logN.col(i)-predN); // N-Process likelihood
    SIMULATE {
      logN.col(i) = predN + neg_log_densityN.simulate();
    }
  }


  // Now finally match to observations
  int f, ft, a, y, amax;
  int minYear=CppAD::Integer((obs(0,0)));
  Type predObs=0, zz, var;
  vector<Type> pred_log(nobs); //
  for(int i=0;i<nobs;i++){
    y=CppAD::Integer(obs(i,0))-minYear;   // 年のラベル
    f=CppAD::Integer(obs(i,1));    //  fleetのラベル
    ft=CppAD::Integer(fleetTypes(f-1));   // fleet typeが何にあたるか?0=caa, 2=survey biomass data, 3=survey SSB data, 4=survey recruitment data, 5=survey SSBm data
    a=CppAD::Integer(obs(i,2))-minAge;  // age
    amax=CppAD::Integer(obs(i,4))-minAge; //maxage
    if(a<(stateDimN-1)){
      zz=exp(logF((keyLogFsta(0,a)),y))+natMor(y,a);  // total mortality
    }else{
      zz=alpha*exp(logF((keyLogFsta(0,a)),y))+natMor(y,a);  // total mortality
    }

    if(ft==0){// residual fleet
      predObs=logN(a,y)-log(zz)+log(1-exp(-zz));
      if((keyLogFsta(f-1,a))>(-1)){
        if(a<(stateDimN-1)){
          predObs+=logF((keyLogFsta(0,a)),y);  // 漁獲方程式
        }else{
          predObs+=log(alpha)+logF((keyLogFsta(0,a)),y);  // 漁獲方程式
        }
      }
    }else{
      if(ft==1){//Not used (same as ft==4)
         predObs=logN(a,y)-zz*sampleTimes(f-1);
          if(CppAD::Integer(keyLogB(f-1,a))>(-1)){
            predObs*=exp(logB(CppAD::Integer(keyLogB(f-1,a))));
          }
          if(CppAD::Integer(keyLogQ(f-1,a))>(-1)){
            predObs+=logQ(CppAD::Integer(keyLogQ(f-1,a)));
          }
      }else{
        if(ft==2){// survey (biomass)
          predObs=0.0;
          for(int j=a; j<amax+1; ++j){
          // for(int j=0; j<stateDimN; ++j){
            // predObs=+exp(logN(j,y))*stockMeanWeight(iy(i),j);
            predObs+=exp(logN(j,y))*stockMeanWeight_true(iy(i),j); //
            }
          predObs=log(predObs)-log(scale);
          // predObs=logN(a,y)+log(stockMeanWeight(iy(i),a));
          if(CppAD::Integer(keyLogB(f-1,a))>(-1)){
            predObs*=exp(logB(CppAD::Integer(keyLogB(f-1,a))));
          }
          if(CppAD::Integer(keyLogQ(f-1,a))>(-1)){
            predObs+=logQ(CppAD::Integer(keyLogQ(f-1,a)));
          }
        }else{
          if(ft==3){// SSB survey
            predObs=0.0;
            for(int j=0; j<stateDimN; ++j){
              predObs+=exp(logN(j,y))*propMat2(iy(i),j)*stockMeanWeight_true(iy(i),j); //
            }
            predObs=log(predObs)-log(scale);
            if(CppAD::Integer(keyLogB(f-1,a))>(-1)){
              predObs*=exp(logB(CppAD::Integer(keyLogB(f-1,a))));
            }
            if(CppAD::Integer(keyLogQ(f-1,a))>(-1)){
              predObs+=logQ(CppAD::Integer(keyLogQ(f-1,a)));
            }
          }else{
            if(ft==4){// Number (e.g.,Recruitment) survey
              predObs=0.0;
              for(int j=a; j<amax+1; ++j){
                predObs=+exp(logN(j,y));
              }
              predObs=log(predObs);
              // predObs=logN(a,y)-zz*sampleTimes(f-1);
              if(CppAD::Integer(keyLogB(f-1,a))>(-1)){
                predObs*=exp(logB(CppAD::Integer(keyLogB(f-1,a))));
              }
              if(CppAD::Integer(keyLogQ(f-1,a))>(-1)){
                predObs+=logQ(CppAD::Integer(keyLogQ(f-1,a)));
              }
            }else{
              if (ft==5){ //Not used
                for(int j=0; j<stateDimN; ++j){
                  predObs+=exp(logN(a+j,y))*exp(-exp(logF((keyLogFsta(0,j)),iy(i)))-natMor(iy(i),j))*propMat2(iy(i),j)*stockMeanWeight_true(iy(i),j);
                }
                predObs=log(predObs);
                if(CppAD::Integer(keyLogB(f-1,a))>(-1)){
                  predObs*=exp(logB(CppAD::Integer(keyLogB(f-1,a))));
                }
                if(CppAD::Integer(keyLogQ(f-1,a))>(-1)){
                  predObs+=logQ(CppAD::Integer(keyLogQ(f-1,a)));
                }
              }else{
                if (ft==6){ // Biomass X selectivity
                  predObs=0.0;
                  for(int j=a; j<amax+1; ++j){
                    // predObs=+exp(logN(j,y))*stockMeanWeight(iy(i),j)*saa(j,y);
                    predObs+=exp(logN(j,y))*stockMeanWeight_true(iy(i),j)*saa(j,y); //
                  }
                  predObs=log(predObs)-log(scale);
                  // predObs=logN(a,y)+log(stockMeanWeight(iy(i),a));
                  if(CppAD::Integer(keyLogB(f-1,a))>(-1)){
                    predObs*=exp(logB(CppAD::Integer(keyLogB(f-1,a))));
                  }
                  if(CppAD::Integer(keyLogQ(f-1,a))>(-1)){
                    predObs+=logQ(CppAD::Integer(keyLogQ(f-1,a)));
                  }
                }else{
                  error("fleet type code not recognized");
                  }
                }
            }
          }
        }
      }
    }
    var=varLogObs(CppAD::Integer(keyVarObs(f-1,a)));
    ans+=-dnorm(log(obs(i,3)),predObs,sqrt(var),true);
    pred_log(i) = predObs;
    SIMULATE {
      obs(i,3) = exp( rnorm(predObs, sqrt(var)) ) ;
    }
  }

  if(b_random==1){
    ans+=-sum(dnorm(logB,0,exp(logSD_b),true));
    // for(int i=0;i<logB.size();i++){
    //   ans+=-dnorm(logB(i),0,exp(logSD_b),true);
    //   }
    }

  for(int i=0;i<timeSteps;i++){
    B_total(i)=0.0;
    F_mean(i)=0.0;
    Catch_biomass(i)=0.0;
    for(int j=0; j<stateDimN; j++){
      B_total(i)+=exp(logN(j,i))*stockMeanWeight_true(i,j);
      F_mean(i)+=exp(logF((keyLogFsta(0,j)),i));
      if (j<(stateDimN-1)) {
        zz=exp(logF((keyLogFsta(0,j)),i))+natMor(i,j);
        Catch_biomass(i)+=exp(logN(j,i))*stockMeanWeight_true(i,j)*exp(logF((keyLogFsta(0,j)),i))/zz;
      } else {
        zz=alpha*exp(logF((keyLogFsta(0,j)),i))+natMor(i,j);
        Catch_biomass(i)+=exp(logN(j,i))*stockMeanWeight_true(i,j)*alpha*exp(logF((keyLogFsta(0,j)),i))/zz;
      }
    }
    F_mean(i)/=stateDimN;
    Exploitation_rate(i)=Catch_biomass(i)/B_total(i);
  }

  ans = (Type(1.0)-lambda)*ans;
  for(int i=0; i<logB.size(); i++){
    ans += lambda*logB(i)*logB(i);
    }

  // Penalty for HS
  if(stockRecruitmentModelCode == 3){
    vector <Type> ssb_vector = ssb.rowwise().sum() ;
    Type maxSSB = max(ssb_vector) ;
    Type minSSB = min(ssb_vector) ;
    Type pen = 0 ;
    pen += CppAD::CondExpLt(exp(rec_logb), maxSSB, Type(0.0), square(rec_logb-log(maxSSB))) ;
    pen += CppAD::CondExpLt(minSSB, exp(rec_logb), Type(0.0), square(log(minSSB)-rec_logb)) ;
    pen *= lambda_Mesnil ;
    ans += pen ;
  }

  // Penalty for Mesnil
  if(stockRecruitmentModelCode == 4){
    vector <Type> ssb_vector = ssb.rowwise().sum() ;
    Type maxSSB = max(ssb_vector) ;
    Type minSSB = min(ssb_vector) ;
    Type pen = 0 ;
    pen += CppAD::CondExpLt(exp(rec_logb), maxSSB, Type(0.0), square(rec_logb-log(maxSSB))) ;
    pen += CppAD::CondExpLt(minSSB, exp(rec_logb), Type(0.0), square(log(minSSB)-rec_logb)) ;
    pen *= lambda_Mesnil ;
    ans += pen ;
  }

  // modeling somatic growth dynamics of weight
  Type sd_w=exp(omicron);
  array<Type> logwaa_pred(logwaa.rows(),logwaa.cols());
  array<Type> waa_true(logwaa.rows(),logwaa.cols());
  logwaa_pred.fill(0.0);
  waa_true.fill(-1);
  Type ans_w=0;
  vector<Type> wp(2);
  Type alpha_w_total, rho_w_total,beta_w0_total;
  vector<Type> N_sum(logwaa.cols());
  Type scale_par=1;
  Type shape=1;

  if(model_weight_maturity(0)==1) {
    for(int j=0;j<logwaa.cols();j++){
      N_sum(j)=Type(0.0);
      // ssb(j)=Type(0.0);
      for(int i=0;i<logwaa.rows();i++){
        waa_true(i,j)=exp(logwaa(i,j));
        N_sum(j)+=exp_logN(i,j)/scale_number;
        // ssb(j)+=exp_logN(i,j)*maa(i,j)*waa_true(i,j)/scale;
        stockMeanWeight_true(j,i)=waa_true(i,j);
        // observation likelihood
        if(i==0) { // age 0
          shape=pow(exp(logCV_w(0)),Type(-2.0));
          if (dist_wobs==0) { //lognormal
            shape=sqrt(log(Type(1.0)+pow(exp(logCV_w(0)),Type(2.0)))); //SD for lognormal distribution
          }
        } else {
          shape=pow(exp(logCV_w(1)),Type(-2.0));
          if (dist_wobs==0) { //lognormal
            shape=sqrt(log(Type(1.0)+pow(exp(logCV_w(1)),Type(2.0)))); //SD for lognormal distribution
          }
        }
        scale_par=waa_true(i,j)/shape;
        if (dist_wobs==1){ //gamma
          ans_w+=-weight_weight(i,j)*dgamma(stockMeanWeight(j,i),shape,scale_par,true); //using Gamma distribution
        }else{ // lognormal
          ans_w+= weight_weight(i,j)*(log(stockMeanWeight(j,i))-dnorm(log(stockMeanWeight(j,i)),logwaa(i,j)-Type(0.5)*shape*shape,shape,true)); //using lognormal distribution with bias correction
        }
      }
    }

    // process model for weight
    for(int j=1;j<logwaa.cols();j++){ //最初の年は除く（2年目から）
      for(int i=0;i<logwaa.rows();i++){
        alpha_w_total=alpha_w(0);
        if(alpha_w.size()>1){
          alpha_w_total+=alpha_w(1)*N_sum(j);
        }
        rho_w_total=rho_w(0);
        beta_w0_total=beta_w0(0);
        if(beta_w0.size()>1){
          beta_w0_total+=beta_w0(1)*ssb(j)/scale;
          // beta_w0_total+=beta_w0(1)*exp_logN(0,j)/scale_number;
          // beta_w0_total+=beta_w0(1)*ssn(j)/scale_number;
        }
        if(i==0){ // age 0
          logwaa_pred(i,j)=beta_w0_total;
        }else{
          if(i<logwaa.rows()-1){
            // from Brody-growth coefficient and the von-Bertalanffy weight model
            logwaa_pred(i,j)=alpha_w_total;
            logwaa_pred(i,j)+=rho_w_total*waa_true(i-1,j-1);
            logwaa_pred(i,j)=log(logwaa_pred(i,j));
          }else{
            if(maxAgePlusGroup==1) { //plus group
              wp(0)=alpha_w_total;
              wp(0)+=rho_w_total*waa_true(i-1,j-1);
              wp(1)=alpha_w_total;
              wp(1)+=rho_w_total*waa_true(i,j-1);
              logwaa_pred(i,j)=exp_logN(i-1,j-1)*wp(0)+exp_logN(i,j-1)*wp(1);
              logwaa_pred(i,j)=logwaa_pred(i,j)/(exp_logN(i-1,j-1)+exp_logN(i,j-1));
              logwaa_pred(i,j)=log(logwaa_pred(i,j));
            }else{
              logwaa_pred(i,j)=alpha_w_total;
              logwaa_pred(i,j)+=rho_w_total*waa_true(i-1,j-1);
              logwaa_pred(i,j)=log(logwaa_pred(i,j));
            }
          }
        }
        ans_w+=-dnorm(logwaa(i,j),logwaa_pred(i,j),sd_w,true);
      }
    }
    ans += ans_w;
  }




  SIMULATE {
    REPORT(logF);
    REPORT(logN);
    REPORT(obs);
  }
  ADREPORT(logN);
  ADREPORT(logF);
  ADREPORT(exp_logN);
  ADREPORT(exp_logF);
  ADREPORT(ssb);
  ADREPORT(B_total);
  ADREPORT(F_mean);
  ADREPORT(phi1);
  ADREPORT(Catch_biomass);
  ADREPORT(Exploitation_rate);
  ADREPORT(stockMeanWeight_true);

  REPORT(logF);
  REPORT(logN);
  REPORT(pred_log);
  REPORT(saa); //selectivity at age
  REPORT(Sel1);
  REPORT(stockMeanWeight_true);
  // REPORT(FY);

  return ans;
}
