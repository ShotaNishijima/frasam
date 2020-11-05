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


template<class Type>
Type objective_function<Type>::operator() ()
{
  // DATA_VECTOR(fleetTypes);
  // DATA_VECTOR(sampleTimes);
  // // DATA_VECTOR(years);
  // DATA_INTEGER(nobs);
  // // DATA_VECTOR(idx1);
  // // DATA_VECTOR(idx2);
  // DATA_ARRAY(obs);
  DATA_ARRAY(propMat);
  DATA_ARRAY(stockMeanWeight);
  DATA_ARRAY(catchMeanWeight);
  DATA_ARRAY(natMor);
  DATA_SCALAR(scale_num_to_mass);
  // DATA_ARRAY(landFrac);
  // DATA_ARRAY(disMeanWeight);
  // DATA_ARRAY(landMeanWeight);
  // DATA_ARRAY(propF);
  // DATA_ARRAY(propM);
  // DATA_INTEGER(minAge);
  // DATA_INTEGER(maxAgePlusGroup);
  // DATA_INTEGER(rhoMode);
  // DATA_IARRAY(keyLogFsta);
  // DATA_ARRAY(keyLogQ);
  // DATA_ARRAY(keyLogB);
  // DATA_ARRAY(keyVarF);
  // DATA_ARRAY(keyVarLogN);
  // DATA_ARRAY(keyVarObs);
  // DATA_INTEGER(stockRecruitmentModelCode);
  // // DATA_VECTOR(fbarRange);
  //
  // PARAMETER_VECTOR(logQ);
  // PARAMETER_VECTOR(logB);
  // PARAMETER_VECTOR(logSdLogFsta);
  // PARAMETER_VECTOR(logSdLogN);
  // PARAMETER_VECTOR(logSdLogObs);
  // PARAMETER(rec_loga);
  // PARAMETER(rec_logb);
  // PARAMETER(logit_rho);
  // PARAMETER_ARRAY(F);
  // PARAMETER_ARRAY(N);
  // PARAMETER(trans_phi1); //　#recruitment autocorreltion
  //
  // DATA_IVECTOR(iy);
  // DATA_INTEGER(nlogF);
  // DATA_INTEGER(nlogN);
  // DATA_ARRAY(propMat2);
  // DATA_SCALAR(alpha);
  // DATA_VECTOR(Fprocess_weight);
  // DATA_SCALAR(F_RW_order); //0: first order, 1: second order
  DATA_IARRAY(Fsigma_key); // ages x years x fleets
  DATA_IARRAY(phi_key); //years x fleets
  DATA_IARRAY(phi_mode); // 0: fixed at 0, 1: fixed at 1, 2: r=phi, 3: r=phi^|i-j|
  DATA_IVECTOR(Fprocess_order); // 1: 1st order RW, 2: 2nd order RW (length: fleet#)
  DATA_IARRAY(Fprocess_weight); //years x fleets
  DATA_IARRAY(Catch_occurrence); //ages x years x fleets
  DATA_INTEGER(recMode); //0: white noise, 1: proportional to SSB
  DATA_IVECTOR(reca_key); //length: years
  DATA_IVECTOR(recb_key); //length: years
  DATA_IARRAY(sigmaN_key); //ages x years
  DATA_IVECTOR(rho_key); //length: years
  DATA_INTEGER(rec_age); //length: years

  PARAMETER_ARRAY(logF); // logF (ages x years x fleets array)
  PARAMETER_ARRAY(logN);  // logN (ages x years matrix)
  PARAMETER_VECTOR(log_sigmaF);
  PARAMETER_VECTOR(trans_phi);
  PARAMETER_VECTOR(rec_loga);
  PARAMETER_VECTOR(rec_logb);
  PARAMETER_VECTOR(log_sigmaN);
  PARAMETER_VECTOR(trans_rho);

  int stateDimF, nyear, nfleet, stateDimN;
  stateDimF = logF.dim(0);
  nyear = logF.dim(1);
  nfleet = logF.dim(2);
  stateDimN = logN.dim(0);

  Type nll=0.0;

  // F processes treating possible multiple fleets

  vector<Type> phi(trans_phi.size());
  for(int i=0; i<phi.size(); i++){
    phi(i)=(exp(trans_phi(i))-Type(1.0))/(exp(trans_phi(i))+Type(1.0));
  }
  vector<Type> sigmaF = exp(log_sigmaF);

  Type nll_Fprocess=0.0;

  //First take care of F
  matrix<Type> fvar(stateDimF,stateDimF);  // F process variance
  matrix<Type> fcor(stateDimF,stateDimF);  // correlation matrix
  vector<Type> fsd(stateDimF);  // Fのsdのvector
  // matrix<Type> logF_mat(stateDimF,nyear);

  using namespace density;  // 多変量正規分布を使う宣言
  MVNORM_t<Type> neg_log_densityF(fvar);  // var-cov matrix fvarを持つMVN

  for(int k=0; k<nfleet; k++){ // per fleet
    for(int t=0; t<nyear; t++){ // per year
      for(int i=0; i<stateDimF; ++i){
        // logF_mat(i,t)=logF(i,t,k);
        fsd(i)=sigmaF(Fsigma_key(i,t,k));
        for(int j=0; j<stateDimF; ++j){
          if(i==j){
            fcor(i,j)=Type(1.0);
          }else{
            if(phi_mode(t,k)==0){
              fcor(i,j)=Type(0.0);
            }else{
              if(phi_mode(t,k)==1){
                fcor(i,j)=Type(1.0);
              }else{
                if(phi_mode(t,k)==2){
                  fcor(i,j)=phi(phi_key(t,k));
                }else{
                  fcor(i,j)=pow(phi(phi_key(t,k)),Type(abs(i-j)));
                }
              }
            }
          }
          fvar(i,j)=fsd(i)*fsd(j)*fcor(i,j);  // var-covを定義
        }
      }
      if(Fprocess_order(k)==1){ //1st RW
        if(t>0){
          nll_Fprocess+=Fprocess_weight(t,k)*neg_log_densityF(logF.col(k).col(t)-logF.col(k).col(t-1)); // F-Process likelihood
          SIMULATE{
            logF.col(k).col(t) = logF.col(k).col(t-1) + neg_log_densityF.simulate();
          }
        }
      }else{
        if(t>1){
          nll_Fprocess=Fprocess_weight(t,k)*neg_log_densityF(logF.col(k).col(t)-Type(2.0)*logF.col(k).col(t-1)+logF.col(k).col(t-2)); // F-Process likelihood
            SIMULATE {
              logF.col(k).col(t) = Type(2.0)*logF.col(k).col(t-1) - logF.col(k).col(t-2) + neg_log_densityF.simulate();
            }
        }
      }
    }
  }

  nll+=nll_Fprocess;

  array<Type> F(stateDimN,nyear,nfleet);
  array<Type> F_total(stateDimF,nyear); // logN (ages x years matrix)
  array<Type> logF_total(stateDimF,nyear); // logN (ages x years matrix)
  F_total.fill(0.0);

  for(int k=0; k<nfleet; k++){ // per fleet
    for(int t=0; t<nyear; t++){ //per age
      for(int i=0;i<stateDimF; i++){
        if(Catch_occurrence(i,t,k)==0){
          F(i,t,k)=Type(0.0);
        }else{
          F(i,t,k)=exp(logF(i,t,k));
        }
        F_total(i,t)+=F(i,t,k);
      }
    }
  }
  logF_total = log(F_total);

  //N process

  Type nll_Nprocess=0.0;

  array<Type> N(stateDimN,nyear);
  vector<Type> SSB(nyear);
  vector<Type> logSSB(nyear);
  vector<Type> B(nyear);
  vector<Type> logB(nyear);

  N=exp(logN);
  SSB.fill(0.0);
  B.fill(0.0);

  for(int i=0;i<nyear;i++){ // calc ssb
    for(int j=0; j<stateDimN; ++j){
      B(i)+=N(j,i)*stockMeanWeight(j,i);
      SSB(i)+=N(j,i)*stockMeanWeight(j,i)*propMat(j,i);
    }
    B(i)/=scale_num_to_mass;
    SSB(i)/=scale_num_to_mass;
    logB(i)=log(B(i));
    logSSB(i)=log(SSB(i));  // log(SSB)
  }

  // For recruitment

  vector<Type> rec_a = exp(rec_loga);
  vector<Type> rec_b = exp(rec_logb);
  vector<Type> rho(trans_rho.size());
  for(int i=0;i<trans_rho.size();i++){
    rho(i)=(exp(trans_rho(i))-Type(1.0))/(exp(trans_rho(i))+Type(1.0));
  }
  array<Type> sigmaN(stateDimN,nyear);
  sigmaN = exp(log_sigmaN);

  //Now take care of N
  matrix<Type> nvar(stateDimN,stateDimN);  // logNのvcov
  MVNORM_t<Type> neg_log_densityN(nvar);

  // vector<Type> R(nyear);
  vector<Type> predR(nyear);
  vector<Type> devianceR(nyear);
  predR.fill(0.0);
  devianceR.fill(0.0);
  int init_year;

  for(int t=0; t<nyear; t++){
    for(int i=0; i<stateDimN; i++) {
      if(i==0){ // recruitment
        switch(recMode){
          case 0: // white noise (constant recruitment)
          init_year=0;
          predR(t) += rec_a(reca_key(t));
          devianceR(t) += logN(i,t)-log(predR(t));
          break;
          
          case 1: //Random walk
          init_year=1;
          if(t>init_year-1){
            predR(t) += logN(i,t-1);
            devianceR(t) += logN(i,t)-log(predR(t));
          }
          break;
          
          case 2: //Proportional to SSB (density independent)
          init_year=rec_age;
            if(t>init_year-1){
              predR(t) += rec_a(reca_key(t))*SSB(t);
              devianceR(t) += logN(i,t)-log(predR(t));
            }
            break;
            
          case 3: //Beverton-Holt
          init_year=rec_age;
          if(t>init_year-1){
            predR(t) += rec_a(reca_key(t))*SSB(t)/(1+rec_b(recb_key(t))*SSB(t));
            devianceR(t) += logN(i,t)-log(predR(t));
          }
          break;
          
        case 4: //Ricker
          init_year=rec_age;
          if(t>init_year-1){
            predR(t) += rec_a(reca_key(t))*SSB(t)*exp(-rec_b(recb_key(t))*SSB(t));
            devianceR(t) += logN(i,t)-log(predR(t));
          }
          break;
          
          default:
          std::cout << "This recMode is not implemented. Input must be 0 to 6." << std::endl;
          exit(EXIT_FAILURE);
          break;
        }
        if(t==init_year){
          nll_Nprocess -= dnorm(logN(i,t),log(predR(t)),pow(sigmaN(sigmaN_key(i,t))/(Type(1.0)-pow(rho(rho_key(t)),Type(2.0))), Type(0.5)),true);
        }else{
          nll_Nprocess -= dnorm(logN(i,t),log(predR(t))+rho(rho_key(t))*devianceR(t-1),sigmaN(sigmaN_key(i,t)),true);
        }
      }
    }
  }


  // vector<Type> predN0(stateDimN);  // logNの予測値
  // vector<Type> predN(stateDimN);  // logNの予測値
  // vector<Type> recResid(timeSteps); //再生産関係からの残差
  //
  // int start_timeStep=1;
  // //Now take care of N
  // matrix<Type> nvar(stateDimN,stateDimN);  // logNのvcov
  // for(int k=0; k<stateDimN; ++k){
  //   for(int j=0; j<stateDimN; ++j){
  //     if(k!=j){nvar(k,j)=0.0;}else{nvar(k,j)=varLogN(CppAD::Integer(keyVarLogN(0,k)));} // logNには相関なし varは加入とそれより上で異なる
  //   }
  // }
  // MVNORM_t<Type> neg_log_densityN(nvar);
  //
  // // //For initial step
  // // matrix<Type> nvar0(stateDimN,stateDimN);  // logNのvcov
  // // for(int k=0; k<stateDimN; ++k){
  // //   for(int j=0; j<stateDimN; ++j){
  // //     if(k!=j){nvar0(k,j)=0.0;}else{
  // //       nvar0(k,j)=varLogN(CppAD::Integer(keyVarLogN(0,k)))/(1-pow(phi1,Type(2.0)));
  // //       } // logNには相関なし varは加入とそれより上で異なる
  // //   }
  // // }
  // // MVNORM_t<Type> neg_log_densityN0(nvar0);
  // Type phi1 = (exp(trans_phi1)-Type(1.0))/(exp(trans_phi1)+Type(1.0));
  //
  // for(int i=start_timeStep;i<timeSteps;i++){
  //   if(stockRecruitmentModelCode==0){ // straight RW
  //     predN0(0)=logN(0,i-1);
  //   }else{
  //     if(stockRecruitmentModelCode==1){//ricker
  //       // predN(0)=rec_loga+log(ssb(i-1))-exp(rec_logb)*ssb(i-1);
  //       predN0(0)=rec_loga+log(ssb(i))-exp(rec_logb)*ssb(i);
  //     }else{
  //       if(stockRecruitmentModelCode==2){//BH
  //         // predN(0)=rec_loga+log(ssb(i-1))-log(1.0+exp(rec_logb)*ssb(i-1));
  //         predN0(0)=rec_loga+log(ssb(i))-log(1.0+exp(rec_logb)*ssb(i));
  //       }else{
  //         if(stockRecruitmentModelCode==3){ //HS
  //           vector<Type> rec_pred_HS(2);
  //           rec_pred_HS(0)=rec_loga+rec_logb;
  //           rec_pred_HS(1)=rec_loga+log(ssb(i));
  //           predN0=min(rec_pred_HS);
  //         } else {
  //           error("SR model code not recognized");
  //         }
  //       }
  //     }
  //   }
  //   recResid(i)=logN(0,i)-predN0(0);
  //   if(i==start_timeStep) {
  //     predN(0)=predN0(0);
  //   } else {
  //     predN(0)=predN0(0)+phi1*recResid(i-1);
  //   }
  //
  //   for(int j=1; j<stateDimN; ++j){
  //     if (j<(stateDimN-1)) {
  //       predN(j)=logN(j-1,i-1)-exp(logF((keyLogFsta(0,j-1)),i-1))-natMor(i-1,j-1);  // population dynamics model
  //     }else{
  //       predN(j)=logN(j-1,i-1)-exp(logF((keyLogFsta(0,j-1)),i-1))-natMor(i-1,j-1);  // population dynamics model
  //     }
  //   }
  //   if(maxAgePlusGroup==1){
  //     predN(stateDimN-1)=log(exp(logN(stateDimN-2,i-1)-exp(logF((keyLogFsta(0,stateDimN-2)),i-1))-natMor(i-1,stateDimN-2))+
  //       exp(logN(stateDimN-1,i-1)-alpha*exp(logF((keyLogFsta(0,stateDimN-1)),i-1))-natMor(i-1,stateDimN-1))); // plus group
  //   }
  //   ans+=neg_log_densityN(logN.col(i)-predN); // N-Process likelihood
  //   SIMULATE {
  //     logN.col(i) = predN + neg_log_densityN.simulate();
  //   }
  // }


  // array<Type> exp_logF(nlogF,U.cols()); // F (6 x 50 matrix)
  // array<Type> exp_logN(nlogN,U.cols()); // N (7 x 50 matrix)
  //
  // for(int i=0;i<nlogN;i++)
  //   for(int j=0;j<U.cols();j++){
  //     logN(i,j)=U(i,j);  // Uの1行からnlogNがlogN
  //     exp_logN(i,j)=exp(logN(i,j));
  //   }
  //   for(int i=0;i<nlogF;i++)
  //     for(int j=0;j<U.cols();j++){
  //       logF(i,j)=U(i+nlogN,j);  // UのnlogN+1からnlogN+nlogFがlogF
  //       exp_logF(i,j)=exp(logF(i,j));
  //     }
  //
  //     int timeSteps=logF.dim[1]; // 年数
  // int stateDimF=logF.dim[0]; // nlogF
  // int stateDimN=logN.dim[0]; // nlogN
  // //Type rho=f(logit_rho);
  // vector<Type> sdLogFsta=exp(logSdLogFsta);  // Fのsd
  // vector<Type> varLogN=exp(logSdLogN*Type(2.0)); //  Nの???散
  // vector<Type> varLogObs=exp(logSdLogObs*Type(2.0)); // Obsの???散
  // vector<Type> ssb(timeSteps); // SSB (年数??????
  // vector<Type> logssb(timeSteps); //  logssb ???年数??????
  // vector<Type> B_total(timeSteps);
  // vector<Type> F_mean(timeSteps);
  //
  //
  //
  //
  // for(int i=0;i<timeSteps;i++){ // calc ssb
  //   ssb(i)=0.0;
  //   for(int j=0; j<stateDimN; ++j){
  //     ssb(i)+=exp(logN(j,i))*exp(-exp(logF((keyLogFsta(0,j)),i))*propF(i,j)-natMor(i,j)*propM(i,j))*propMat(i,j)*stockMeanWeight(i,j);  // ssbを??????
  //     // ssb(i)+=exp(logN(j,i))*propMat(i,j)*stockMeanWeight(i,j);  // ssbを??????
  //   }
  //   logssb(i)=log(ssb(i));  // log(ssb)
  // }
  //
  // // for(int i=0;i<timeSteps;i++){ // calc total biomass
  // //   B_total(i)=0.0;
  // //   for(int j=0; j<stateDimN; ++j){
  // //     B_total(i)+=exp(logN(j,i))*exp(-exp(logF((keyLogFsta(0,j)),i))*propF(i,j)-natMor(i,j)*propM(i,j))*propMat(i,j)*stockMeanWeight(i,j);  // ssbを??????
  // //     // ssb(i)+=exp(logN(j,i))*propMat(i,j)*stockMeanWeight(i,j);  // ssbを??????
  // //   }
  // //   logssb(i)=log(ssb(i));  // log(ssb)
  // // }
  //
  // vector<Type> predN0(stateDimN);  // logNの予測値
  // vector<Type> predN(stateDimN);  // logNの予測値
  // vector<Type> recResid(timeSteps); //再生産関係からの残差
  //
  // int start_timeStep=1;
  // //Now take care of N
  // matrix<Type> nvar(stateDimN,stateDimN);  // logNのvcov
  // for(int k=0; k<stateDimN; ++k){
  //   for(int j=0; j<stateDimN; ++j){
  //     if(k!=j){nvar(k,j)=0.0;}else{nvar(k,j)=varLogN(CppAD::Integer(keyVarLogN(0,k)));} // logNには相関なし varは加入とそれより上で異なる
  //   }
  // }
  // MVNORM_t<Type> neg_log_densityN(nvar);
  //
  // // //For initial step
  // // matrix<Type> nvar0(stateDimN,stateDimN);  // logNのvcov
  // // for(int k=0; k<stateDimN; ++k){
  // //   for(int j=0; j<stateDimN; ++j){
  // //     if(k!=j){nvar0(k,j)=0.0;}else{
  // //       nvar0(k,j)=varLogN(CppAD::Integer(keyVarLogN(0,k)))/(1-pow(phi1,Type(2.0)));
  // //       } // logNには相関なし varは加入とそれより上で異なる
  // //   }
  // // }
  // // MVNORM_t<Type> neg_log_densityN0(nvar0);
  // Type phi1 = (exp(trans_phi1)-Type(1.0))/(exp(trans_phi1)+Type(1.0));
  //
  // for(int i=start_timeStep;i<timeSteps;i++){
  //   if(stockRecruitmentModelCode==0){ // straight RW
  //     predN0(0)=logN(0,i-1);
  //   }else{
  //     if(stockRecruitmentModelCode==1){//ricker
  //       // predN(0)=rec_loga+log(ssb(i-1))-exp(rec_logb)*ssb(i-1);
  //       predN0(0)=rec_loga+log(ssb(i))-exp(rec_logb)*ssb(i);
  //     }else{
  //       if(stockRecruitmentModelCode==2){//BH
  //         // predN(0)=rec_loga+log(ssb(i-1))-log(1.0+exp(rec_logb)*ssb(i-1));
  //         predN0(0)=rec_loga+log(ssb(i))-log(1.0+exp(rec_logb)*ssb(i));
  //       }else{
  //         if(stockRecruitmentModelCode==3){ //HS
  //           vector<Type> rec_pred_HS(2);
  //           rec_pred_HS(0)=rec_loga+rec_logb;
  //           rec_pred_HS(1)=rec_loga+log(ssb(i));
  //           predN0=min(rec_pred_HS);
  //         } else {
  //           error("SR model code not recognized");
  //         }
  //       }
  //     }
  //   }
  //   recResid(i)=logN(0,i)-predN0(0);
  //   if(i==start_timeStep) {
  //     predN(0)=predN0(0);
  //   } else {
  //     predN(0)=predN0(0)+phi1*recResid(i-1);
  //   }
  //
  //   for(int j=1; j<stateDimN; ++j){
  //     if (j<(stateDimN-1)) {
  //       predN(j)=logN(j-1,i-1)-exp(logF((keyLogFsta(0,j-1)),i-1))-natMor(i-1,j-1);  // population dynamics model
  //     }else{
  //       predN(j)=logN(j-1,i-1)-exp(logF((keyLogFsta(0,j-1)),i-1))-natMor(i-1,j-1);  // population dynamics model
  //     }
  //   }
  //   if(maxAgePlusGroup==1){
  //     predN(stateDimN-1)=log(exp(logN(stateDimN-2,i-1)-exp(logF((keyLogFsta(0,stateDimN-2)),i-1))-natMor(i-1,stateDimN-2))+
  //       exp(logN(stateDimN-1,i-1)-alpha*exp(logF((keyLogFsta(0,stateDimN-1)),i-1))-natMor(i-1,stateDimN-1))); // plus group
  //   }
  //   ans+=neg_log_densityN(logN.col(i)-predN); // N-Process likelihood
  //   SIMULATE {
  //     logN.col(i) = predN + neg_log_densityN.simulate();
  //   }
  // }
  //
  //
  // // Now finally match to observations
  // int f, ft, a, y;
  // int minYear=CppAD::Integer((obs(0,0)));
  // Type predObs=0, zz, var;
  // for(int i=0;i<nobs;i++){
  //   y=CppAD::Integer(obs(i,0))-minYear;   // 年のラベル
  //   f=CppAD::Integer(obs(i,1));    //  fleetのラベル
  //   ft=CppAD::Integer(fleetTypes(f-1));   // fleet typeが何にあたるか???0=caa???2=survey biomass data, 3=survey SSB data, 4=survey recruitment data, 5=survey SSBm data???
  //   a=CppAD::Integer(obs(i,2))-minAge;  // age
  //   if(a<(stateDimN-1)){
  //     zz=exp(logF((keyLogFsta(0,a)),y))+natMor(y,a);  // total mortality
  //   }else{
  //     zz=alpha*exp(logF((keyLogFsta(0,a)),y))+natMor(y,a);  // total mortality
  //   }
  //
  //   if(ft==0){// residual fleet
  //     predObs=logN(a,y)-log(zz)+log(1-exp(-zz));
  //     if((keyLogFsta(f-1,a))>(-1)){
  //       if(a<(stateDimN-1)){
  //         predObs+=logF((keyLogFsta(0,a)),y);  // 漁獲方??????
  //       }else{
  //         predObs+=log(alpha)+logF((keyLogFsta(0,a)),y);  // 漁獲方??????
  //       }
  //     }
  //   }else{
  //     if(ft==1){// comm fleet
  //       predObs=logN(a,y)-zz*sampleTimes(f-1);
  //       if(CppAD::Integer(keyLogB(f-1,a))>(-1)){
  //         predObs*=exp(logB(CppAD::Integer(keyLogB(f-1,a))));
  //       }
  //       if(CppAD::Integer(keyLogQ(f-1,a))>(-1)){
  //         predObs+=logQ(CppAD::Integer(keyLogQ(f-1,a)));
  //       }
  //     }else{
  //       if(ft==2){// survey (biomass)
  //         predObs=logN(a,y)+log(stockMeanWeight(iy(i),a));
  //         if(CppAD::Integer(keyLogB(f-1,a))>(-1)){
  //           predObs*=exp(logB(CppAD::Integer(keyLogB(f-1,a))));
  //         }
  //         if(CppAD::Integer(keyLogQ(f-1,a))>(-1)){
  //           predObs+=logQ(CppAD::Integer(keyLogQ(f-1,a)));
  //         }
  //       }else{
  //         if(ft==3){// SSB survey
  //           for(int j=0; j<stateDimN; ++j){
  //             predObs+=exp(logN(a+j,y))*propMat2(iy(i),j)*stockMeanWeight(iy(i),j); //
  //           }
  //           predObs=log(predObs);
  //           if(CppAD::Integer(keyLogB(f-1,a))>(-1)){
  //             predObs*=exp(logB(CppAD::Integer(keyLogB(f-1,a))));
  //           }
  //           if(CppAD::Integer(keyLogQ(f-1,a))>(-1)){
  //             predObs+=logQ(CppAD::Integer(keyLogQ(f-1,a)));
  //           }
  //         }else{
  //           if(ft==4){// Recruitment survey
  //             predObs=logN(a,y)-zz*sampleTimes(f-1);
  //             if(CppAD::Integer(keyLogB(f-1,a))>(-1)){
  //               predObs*=exp(logB(CppAD::Integer(keyLogB(f-1,a))));
  //             }
  //             if(CppAD::Integer(keyLogQ(f-1,a))>(-1)){
  //               predObs+=logQ(CppAD::Integer(keyLogQ(f-1,a)));
  //             }
  //           }else{
  //             if (ft==5){
  //               for(int j=0; j<stateDimN; ++j){
  //                 predObs+=exp(logN(a+j,y))*exp(-exp(logF((keyLogFsta(0,j)),iy(i)))-natMor(iy(i),j))*propMat2(iy(i),j)*stockMeanWeight(iy(i),j);
  //               }
  //               predObs=log(predObs);
  //               if(CppAD::Integer(keyLogB(f-1,a))>(-1)){
  //                 predObs*=exp(logB(CppAD::Integer(keyLogB(f-1,a))));
  //               }
  //               if(CppAD::Integer(keyLogQ(f-1,a))>(-1)){
  //                 predObs+=logQ(CppAD::Integer(keyLogQ(f-1,a)));
  //               }
  //             }
  //           }
  //         }
  //       }
  //     }
  //   }
  //   var=varLogObs(CppAD::Integer(keyVarObs(f-1,a)));
  //   ans+=-dnorm(log(obs(i,3)),predObs,sqrt(var),true);
  //   SIMULATE {
  //     obs(i,3) = exp( rnorm(predObs, sqrt(var)) ) ;
  //   }
  // }
  //
  // for(int i=0;i<timeSteps;i++){
  //   B_total(i)=0.0;
  //   F_mean(i)=0.0;
  //   for(int j=0; j<stateDimN; ++j){
  //     B_total(i)+=exp(logN(j,i))*stockMeanWeight(i,j);
  //     F_mean(i)+=exp(logF((keyLogFsta(0,j)),i));
  //   }
  //   F_mean(i)/=stateDimN;
  // }
  //
  SIMULATE {
    REPORT(logF);
  //   REPORT(logN);
  //   REPORT(obs);
  }
  ADREPORT(logN);
  ADREPORT(logF);
  ADREPORT(N);
  ADREPORT(F);
  ADREPORT(nll);
  ADREPORT(nll_Fprocess);
  // ADREPORT(B_total);
  // ADREPORT(F_mean);
  // ADREPORT(phi1);

  return nll;
}
