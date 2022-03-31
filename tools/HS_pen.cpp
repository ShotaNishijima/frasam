
#include <TMB.hpp>
#include <iostream>

template <class Type>
Type square(Type x){return x*x;}

template<class Type>
Type objective_function<Type>::operator() ()
{
  DATA_VECTOR(SSB);
  DATA_VECTOR(R);

  PARAMETER(a);
  PARAMETER(b);

  // Type ny=SSB.size();
  Type ans=0;
  Type lambda=100000;
  vector<Type> predR(SSB.size());

  for(int i=0; i<SSB.size(); ++i){
    predR(i)=CppAD::CondExpLt(b,SSB(i),a*b,a*SSB(i));
    ans+=square(log(R(i))-log(predR(i)));
  }

  Type maxSSB=max(SSB);
  Type minSSB=min(SSB);
  Type pen=0;

  // pen+=CppAD::CondExpLt(b,maxSSB,Type(0.0),square(b-maxSSB));
  // pen+=CppAD::CondExpLt(minSSB,b,Type(0.0),square(minSSB-b));
  // pen+=CppAD::CondExpLt(b,maxSSB,Type(0.0),b-maxSSB);
  // pen+=CppAD::CondExpLt(minSSB,b,Type(0.0),minSSB-b);
  pen+=CppAD::CondExpLt(b,maxSSB,Type(0.0),square(log(b)-log(maxSSB)));
  pen+=CppAD::CondExpLt(minSSB,b,Type(0.0),square(log(minSSB)-log(b)));

  pen*=lambda;

  ans+=pen;

  return ans;
}
