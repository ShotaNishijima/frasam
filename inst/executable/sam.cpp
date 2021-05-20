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

  PARAMETER_VECTOR(logQ);
  PARAMETER_VECTOR(logB);
  PARAMETER_VECTOR(logSdLogFsta);
  PARAMETER_VECTOR(logSdLogN);
  PARAMETER_VECTOR(logSdLogObs);
  PARAMETER(rec_loga);
  PARAMETER(rec_logb);
  PARAMETER(logit_rho);
  PARAMETER_ARRAY(U);
  PARAMETER(trans_phi1); //　#recruitment autocorreltion

  DATA_IVECTOR(iy);
  DATA_INTEGER(nlogF);
  DATA_INTEGER(nlogN);
  DATA_ARRAY(propMat2);
  DATA_SCALAR(alpha);
  DATA_VECTOR(Fprocess_weight);
  DATA_SCALAR(F_RW_order); //0: first order, 1: second order

  array<Type> logF(nlogF,U.cols()); // logF (6 x 50 matrix)
  array<Type> logN(nlogN,U.cols()); // logN (7 x 50 matrix)
  array<Type> exp_logF(nlogF,U.cols()); // F (6 x 50 matrix)
  array<Type> exp_logN(nlogN,U.cols()); // N (7 x 50 matrix)

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
      if(i!=j){if(rhoMode==0){fcor(i,j)=0.01;
      }else{if(rhoMode==1){fcor(i,j)=0.99;
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
    for(int j=0; j<stateDimN; ++j){ 
      ssb(i)+=exp(logN(j,i))*exp(-exp(logF((keyLogFsta(0,j)),i))*propF(i,j)-natMor(i,j)*propM(i,j))*propMat(i,j)*stockMeanWeight(i,j);  // ssbを??????
      // ssb(i)+=exp(logN(j,i))*propMat(i,j)*stockMeanWeight(i,j);  // ssbを??????
    }
    logssb(i)=log(ssb(i));  // log(ssb)
  }

  // for(int i=0;i<timeSteps;i++){ // calc total biomass
  //   B_total(i)=0.0;
  //   for(int j=0; j<stateDimN; ++j){
  //     B_total(i)+=exp(logN(j,i))*exp(-exp(logF((keyLogFsta(0,j)),i))*propF(i,j)-natMor(i,j)*propM(i,j))*propMat(i,j)*stockMeanWeight(i,j);  // ssbを??????
  //     // ssb(i)+=exp(logN(j,i))*propMat(i,j)*stockMeanWeight(i,j);  // ssbを??????
  //   }
  //   logssb(i)=log(ssb(i));  // log(ssb)
  // }

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
              error("SR model code not recognized");
              }
          }
        }
      }
    }
    recResid(i)=logN(0,i)-predN0(0);
    if(i==start_timeStep) {
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
            predObs=+exp(logN(j,y))*stockMeanWeight(iy(i),j);
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
            for(int j=0; j<stateDimN; ++j){
              predObs+=exp(logN(a+j,y))*propMat2(iy(i),j)*stockMeanWeight(iy(i),j); //
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
              predObs=log(predObs)-log(scale);
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
                  predObs+=exp(logN(a+j,y))*exp(-exp(logF((keyLogFsta(0,j)),iy(i)))-natMor(iy(i),j))*propMat2(iy(i),j)*stockMeanWeight(iy(i),j);
                }
                predObs=log(predObs);
                if(CppAD::Integer(keyLogB(f-1,a))>(-1)){
                  predObs*=exp(logB(CppAD::Integer(keyLogB(f-1,a))));
                }
                if(CppAD::Integer(keyLogQ(f-1,a))>(-1)){
                  predObs+=logQ(CppAD::Integer(keyLogQ(f-1,a)));
                }
              }
            }
          }
        }
      }
    }
    var=varLogObs(CppAD::Integer(keyVarObs(f-1,a)));
    ans+=-dnorm(log(obs(i,3)),predObs,sqrt(var),true);
    SIMULATE {
      obs(i,3) = exp( rnorm(predObs, sqrt(var)) ) ;
    }
  }

  for(int i=0;i<timeSteps;i++){
    B_total(i)=0.0;
    F_mean(i)=0.0;
    Catch_biomass(i)=0.0;
    for(int j=0; j<stateDimN; ++j){
      B_total(i)+=exp(logN(j,i))*stockMeanWeight(i,j);
      F_mean(i)+=exp(logF((keyLogFsta(0,j)),i));
      if (j<(stateDimN-1)) {
        zz=exp(logF((keyLogFsta(0,j)),i))+natMor(i,j);
        Catch_biomass(i)+=exp(logN(j,i))*stockMeanWeight(i,j)*exp(logF((keyLogFsta(0,j)),i))/zz;
      } else {
        zz=alpha*exp(logF((keyLogFsta(0,j)),i))+natMor(i,j);
        Catch_biomass(i)+=exp(logN(j,i))*stockMeanWeight(i,j)*alpha*exp(logF((keyLogFsta(0,j)),i))/zz;
      }
    }
    F_mean(i)/=stateDimN;
    Exploitation_rate(i)=Catch_biomass(i)/B_total(i);
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

  return ans;
}
