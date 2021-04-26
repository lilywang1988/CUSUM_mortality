#include <Rcpp.h>
//#include <Rcpp/vector.h>
using namespace Rcpp;

double closest_t_cpp(NumericVector v,double tau0, double tau, int ndays);

// [[Rcpp::export]]
NumericVector O_E_CUSUM_hval_resample_cpp(int nloop,int n, NumericVector theta,
                                 NumericVector enrl_t, NumericVector chro_time,
                                 LogicalVector delta_vec,
                                 NumericVector xbeta,
                                 NumericVector lambda0 = 0.0004453,
                                 double tau0 = 1, 
                                 double tau = 1461, 
                                 double p = 0.95){
  int N = enrl_t.length();
  if(chro_time.length() != N||delta_vec.length() != N || xbeta.length() != N){
    stop("Different lengths of the following vectors: enrl_t, chro_time, delta_vec and xbeta");
  }
  int ndays = trunc(tau - tau0 + 1);
  double h_temp;
  NumericVector k_vec = (exp(theta) - 1)/(theta) - 1;
  IntegerVector sign_k_vec = sign(k_vec);
  NumericVector abs_k_vec = abs(k_vec);
  NumericVector h_sel(k_vec.length());
  int sN = floor(nloop * (1 - p));
  NumericMatrix h_quantiles(sN,k_vec.length());
  
  lambda0=rep_len(lambda0, ndays);
  NumericVector Lambda0=cumsum(lambda0);

   LogicalVector delta_vec_use(n);
   NumericVector chro_time_use(n);
   NumericVector enrl_t_use(n);
   NumericVector xbeta_use(n);
   IntegerVector select(n);

   for(int loop=1;loop<=nloop; loop++){
     //Rcout<<loop<<std::endl;
     select = sample(N,n,TRUE);
     
     xbeta_use = xbeta[select-1];
     enrl_t_use = enrl_t[select-1];
     chro_time_use = chro_time[select-1];
     delta_vec_use = delta_vec[select-1];
     
     NumericMatrix E_t_pre(n,ndays);
     NumericMatrix O_t_pre(n,ndays);
     NumericVector M_t_pre(ndays);
     NumericVector M_t_pre2(ndays);
     
     LogicalVector at_risk;
     at_risk = !((chro_time_use<tau0)|(enrl_t_use>tau));
     int true_in;
     double deduct;
     int index;
     for (int i=0;i<n;i++){
       if(at_risk[i]){
         //int true_fu=(start_day>enrl_t[i])?(time[i]-(start_day-enrl_t[i])):time[i];
         true_in = tau0>=enrl_t_use[i]?tau0:enrl_t_use[i];//enrl_t[i]
         //int true_in=enrl_t[i];
         for(int j=0;j<ndays;j++){
           if(j>=chro_time_use[i]-tau0 && delta_vec_use[i]){
             O_t_pre(i,j)=1;
           }
           index=int(true_in-2);
           deduct =true_in<=1?0:Lambda0[index];
           if(j>=chro_time_use[i]-tau0) {
             index=int(chro_time_use[i]-1);
             E_t_pre(i,j)=(Lambda0[index]-deduct)*exp(xbeta_use[i]);
           } else if(j>=true_in-tau0 && j<=chro_time_use[i]-tau0){
             index=int(j+tau0-1);
             E_t_pre(i,j)=(Lambda0[index]-deduct)*exp(xbeta_use[i]);
           }
         }
       }
     }
     NumericVector O_t(ndays);
     NumericVector E_t(ndays);
     NumericVector O_E_t(ndays);
     for(int j=0; j<ndays; j++){
       O_t[j]=sum(O_t_pre.column(j));
       E_t[j]=sum(E_t_pre.column(j));
     }
     O_E_t=O_t-E_t;
    //int sumit=sum(delta_vec);
    //Rcout<< E_t[ndays-1]<<std::endl;
    //Rcout<< O_t[ndays-1]<<std::endl;
    //print(delta_vec);
    //print(O_t);
    //print(O_E_t);
     for(int k=0;k<k_vec.length();k++){
       M_t_pre=sign_k_vec[k]*O_E_t-abs_k_vec[k]*E_t;
       M_t_pre.push_front(0);
       M_t_pre2 = NumericVector(cummin(M_t_pre));
       h_temp=max(M_t_pre-M_t_pre2);
       //Rcout<<h_temp<<std::endl;
       bool mainloop=false;
       for(int sn=0; sn<sN && mainloop==false; sn++){
         if(h_temp>h_quantiles(sn,k)){mainloop=true;
           for(int ssn=sN-1;ssn>sn;ssn--){h_quantiles(ssn,k)=h_quantiles(ssn-1,k);}
           h_quantiles(sn,k)=h_temp;
         }
       }
       // Rcout<<h_temp<<std::endl;
       // print(h_quantiles);
     }
     
  }
  //print(h_quantiles);
 return( h_quantiles.row(sN-1));
//  return(Lambda0);
}
//' Simulation via resmapling
//' 
// [[Rcpp::export]]
List O_E_CUSUM_sim_resample_cpp(int nloop,int n, NumericVector theta,
                       NumericVector rho_vec,
                       NumericVector h,
                       NumericVector enrl_t, NumericVector chro_time,
                       LogicalVector delta_vec,
                       NumericVector xbeta,
                       NumericVector lambda0= 0.0004453,
                       double tau0=1, double tau=1461){

  if(sum((rho_vec>1)|(rho_vec<0))!=0) {rho_vec=abs(rho_vec/max(rho_vec));
    warning("Rho_list has been corrected.");}//correct rho vavalues if invalid
  int N=enrl_t.length();
  if(chro_time.length()!=N||delta_vec.length()!=N||xbeta.length()!=N){
    stop("Different lengths of the following vectors: enrl_t, chro_time, delta_vec and xbeta");
  }
  int ndays = trunc(tau-tau0+1);
  IntegerVector time_list=seq_len(ndays);
  time_list.push_front(0);
  int hit_time;
  NumericVector k_vec=(exp(theta)-1)/(theta)-1;
  IntegerVector sign_k_vec = sign(k_vec);
  NumericVector abs_k_vec = abs(k_vec);
  
  NumericMatrix hit_table(rho_vec.length(),k_vec.length());
  NumericMatrix hit_utility(rho_vec.length(),k_vec.length());
  NumericMatrix hit_total(rho_vec.length(),k_vec.length());
  
 
  lambda0=rep_len(lambda0, ndays);
  NumericVector Lambda0=cumsum(lambda0);
  

  LogicalVector delta_vec_use(n);
  NumericVector chro_time_use(n);
  NumericVector enrl_t_use(n);
  NumericVector xbeta_use(n);
  IntegerVector select(n);
  
  for(int loop=1;loop<=nloop; loop++){
    //Rcout<<loop<<std::endl;
    select = sample(N,n,TRUE);
    xbeta_use = xbeta[select-1];
    enrl_t_use = enrl_t[select-1];
    chro_time_use = chro_time[select-1];
    delta_vec_use = delta_vec[select-1];
    
    NumericMatrix E_t_pre(n,ndays);
    NumericMatrix O_t_pre(n,ndays);
    NumericVector M_t_pre;
    NumericVector M_t_pre2;
    NumericVector M_t;
    NumericVector M_t2;
    NumericVector M_t_rho; 
    
    LogicalVector at_risk;
    at_risk = !((chro_time_use<tau0)|(enrl_t_use>tau));
    int true_in;
    double deduct;
    int index;
    for (int i=0;i<n;i++){
      if(at_risk[i]){
        //int true_fu=(start_day>enrl_t[i])?(time[i]-(start_day-enrl_t[i])):time[i];
        true_in = tau0>=enrl_t_use[i]?tau0:enrl_t_use[i];//enrl_t[i]
        //int true_in=enrl_t[i];
        for(int j=0;j<ndays;j++){
          if(j>=chro_time_use[i]-tau0 && delta_vec_use[i]){
            O_t_pre(i,j)=1;
          }
          index=int(true_in-2);
          deduct =true_in<=1?0:Lambda0[index];
          if(j>=chro_time_use[i]-tau0) {
            index=int(chro_time_use[i]-1);
            E_t_pre(i,j)=(Lambda0[index]-deduct)*exp(xbeta_use[i]);
          } else if(j>=true_in-tau0 && j<=chro_time_use[i]-tau0){
            index=int(j+tau0-1);
            E_t_pre(i,j)=(Lambda0[index]-deduct)*exp(xbeta_use[i]);
          }
        }
      }
    }
    NumericVector O_t(ndays);
    NumericVector E_t(ndays);
    NumericVector O_E_t(ndays);
    for(int j=0; j<ndays; j++){
      O_t[j]=sum(O_t_pre.column(j));
      E_t[j]=sum(E_t_pre.column(j));
    }
    O_E_t=O_t-E_t;
    //int sumit=sum(delta_vec);
    // Rcout<< O_E_t[ndays-1]<<std::endl;
    //print(delta_vec);
    //print(O_t);
    //print(O_E_t);
    for(int k=0;k<k_vec.length();k++){
      M_t_pre=sign_k_vec(k)*O_E_t-abs_k_vec(k)*E_t;
      //Rcout<<"M_t_pre.length()="<<M_t_pre.length()<<std::endl;
      M_t_pre.push_front(0);
      //Rcout<<"M_t_pre.length()="<<M_t_pre.length()<<std::endl;
      M_t_pre2 = NumericVector(cummin(M_t_pre));
      M_t = M_t_pre2-M_t_pre+h(k);
      //print(M_t);
      for(int r=0;r<rho_vec.length();r++){
        M_t2=(1-rho_vec(r))*h(k)-M_t_pre;
        //print(M_t2);
        M_t_rho=pmin(M_t,M_t2);
        //Rcout<<"min(M_t_rho)="<<min(M_t_rho)<<' '<<"M_t_rho(0)="<<M_t_rho(0)<<" "<<"M_t(0)"<<M_t(0)<<" "<<"M_t_pre(0)"<<M_t_pre(0)<<" "<<"k="<<k<<"rho="<<rho_vec(r)<<std::endl;
        if(min(M_t_rho)<0){
          // if(k==0) Rcout<<"yes"<< r<<std::endl;
          hit_total(r,k) += 1;
          hit_time = min(as<NumericVector>(time_list[M_t_rho<0]));
          //Rcout<<hit_time<<std::endl;
        } else {hit_time = ndays;}
        hit_table(r,k) += hit_time;
        //Rcout<<"yes"<< hit_table(r,k)<<std::endl;
        // print(time_list[M_t_rho.row(k)<0]);
        //Rcout<<"hittime="<<hit_time<<"k="<<k<<"r="<<r<<std::endl;
        
        if(loop == nloop){
          hit_table(r,k) = hit_table(r,k)/nloop;
          hit_utility(r,k)=hit_total(r,k)/nloop;
        }
      }
      
    }
    
  }
  
 return(List::create(Named("hit_table")=hit_table,_["hit_utility"]=hit_utility));
//return(List::create(_["hello"]=0));
}
//' Simulation via resampling and record the second signal
//' 
// [[Rcpp::export]]
List O_E_CUSUM_sim_resample2_cpp(int nloop,int n, NumericVector theta,
                                NumericVector rho_vec,
                                NumericVector h,
                                NumericVector enrl_t, NumericVector chro_time,
                                LogicalVector delta_vec,
                                NumericVector xbeta,
                                NumericVector lambda0 = 0.0004453,
                                double tau0 = 1, double tau = 1461){
  
  if(sum((rho_vec>1)|(rho_vec<0))!=0) {rho_vec=abs(rho_vec/max(rho_vec));
    warning("Rho_list has been corrected.");}//correct rho vavalues if invalid
  int N=enrl_t.length();
  if(chro_time.length()!=N||delta_vec.length()!=N||xbeta.length()!=N){
    stop("Different lengths of the following vectors: enrl_t, chro_time, delta_vec and xbeta");
  }
  int ndays = trunc(tau-tau0+1);
  IntegerVector time_list=seq_len(ndays);
  time_list.push_front(0);
  int hit_time;
  int hit_time2;
  NumericVector k_vec=(exp(theta)-1)/(theta)-1;
  IntegerVector sign_k_vec = sign(k_vec);
  NumericVector abs_k_vec = abs(k_vec);
  
  NumericMatrix hit_table(1,k_vec.length());
  NumericMatrix hit_utility(1,k_vec.length());
  NumericMatrix hit_total(1,k_vec.length());
  
  NumericMatrix hit_table2(rho_vec.length(),k_vec.length());
  NumericMatrix hit_utility2(rho_vec.length(),k_vec.length());
  NumericMatrix hit_total2(rho_vec.length(),k_vec.length());
  
  
  lambda0=rep_len(lambda0, ndays);
  NumericVector Lambda0=cumsum(lambda0);
  
  
  LogicalVector delta_vec_use(n);
  NumericVector chro_time_use(n);
  NumericVector enrl_t_use(n);
  NumericVector xbeta_use(n);
  IntegerVector select(n);
  
  for(int loop=1;loop<=nloop; loop++){
    //Rcout<<loop<<std::endl;
    select = sample(N,n,TRUE);
    xbeta_use = xbeta[select-1];
    enrl_t_use = enrl_t[select-1];
    chro_time_use = chro_time[select-1];
    delta_vec_use = delta_vec[select-1];
    
    NumericMatrix E_t_pre(n,ndays);
    NumericMatrix O_t_pre(n,ndays);
    NumericVector M_t_pre;
    NumericVector M_t_pre2;
    NumericVector M_t;
    NumericVector M_t2;
    NumericVector M_t_rho; 
    NumericVector M_t_tmp;
    NumericVector M_t2_tmp;
    NumericVector time_list_tmp;
    NumericVector M_t_pre2_tmp;
    
    LogicalVector at_risk;
    at_risk = !((chro_time_use<tau0)|(enrl_t_use>tau));
    int true_in;
    double deduct;
    int index;
    for (int i=0;i<n;i++){
      if(at_risk[i]){
        //int true_fu=(start_day>enrl_t[i])?(time[i]-(start_day-enrl_t[i])):time[i];
        true_in = tau0>=enrl_t_use[i]?tau0:enrl_t_use[i];//enrl_t[i]
        //int true_in=enrl_t[i];
        for(int j=0;j<ndays;j++){
          if(j>=chro_time_use[i]-tau0 && delta_vec_use[i]){
            O_t_pre(i,j)=1;
          }
          index=int(true_in-2);
          deduct =true_in<=1?0:Lambda0[index];
          if(j>=chro_time_use[i]-tau0) {
            index=int(chro_time_use[i]-1);
            E_t_pre(i,j)=(Lambda0[index]-deduct)*exp(xbeta_use[i]);
          } else if(j>=true_in-tau0 && j<=chro_time_use[i]-tau0){
            index=int(j+tau0-1);
            E_t_pre(i,j)=(Lambda0[index]-deduct)*exp(xbeta_use[i]);
          }
        }
      }
    }
    NumericVector O_t(ndays);
    NumericVector E_t(ndays);
    NumericVector O_E_t(ndays);
    for(int j=0; j<ndays; j++){
      O_t[j]=sum(O_t_pre.column(j));
      E_t[j]=sum(E_t_pre.column(j));
    }
    O_E_t=O_t-E_t;
    //int sumit=sum(delta_vec);
    // Rcout<< O_E_t[ndays-1]<<std::endl;
    //print(delta_vec);
    //print(O_t);
    //print(O_E_t);
    for(int k=0;k<k_vec.length();k++){
      M_t_pre=sign_k_vec(k)*O_E_t-abs_k_vec(k)*E_t;
      //Rcout<<"M_t_pre.length()="<<M_t_pre.length()<<std::endl;
      M_t_pre.push_front(0);
      //Rcout<<"M_t_pre.length()="<<M_t_pre.length()<<std::endl;
      M_t_pre2 = NumericVector(cummin(M_t_pre));
      M_t = M_t_pre2-M_t_pre + h(k);
      //print(M_t);
      if(min(M_t)<0){
        hit_total(0,k) += 1;
        hit_time = min(as<NumericVector>(time_list[M_t < 0]));
        
        if(hit_time<ndays-1){
          //Rcout<<"hit_time"<<std::endl;
          time_list_tmp = seq(hit_time+1, ndays);
          //Rcout <<"ndays"<<ndays;
          
          for(int r=0; r < rho_vec.length(); r++){
            //Rcout << M_t_pre2_tmp.length()<<std::endl;
            M_t_pre2_tmp = NumericVector(cummin(M_t_pre[seq(hit_time, ndays-1)]));
            M_t_tmp = M_t_pre2_tmp - M_t_pre[seq(hit_time, ndays-1)] + h(k);
            M_t2_tmp = (1-rho_vec(r))*h(k)-M_t_pre[seq(hit_time, ndays-1)]+M_t_pre[hit_time-1];
            M_t_rho = pmin(M_t_tmp,M_t2_tmp);
            if(min(M_t_rho)<0){
              hit_total2(r,k) += 1;
              // if(k==0) Rcout<<"yes"<< r<<std::endl;
              hit_time2=min(as<NumericVector>(time_list_tmp[M_t_rho<0]));
              //Rcout<<"hit_time2"<<std::endl;
            }else{ hit_time2 = ndays;}
            
            hit_table2(r,k) += hit_time2;
            //Rcout<<"yes"<< hit_table(r,k)<<std::endl;
            // print(time_list[M_t_rho.row(k)<0]);
            //Rcout<<"hittime="<<hit_time<<"k="<<k<<"r="<<r<<std::endl;
            
          }
        }else{
          for(int r=0; r < rho_vec.length(); r++){
            hit_table2(r,k) += ndays;
          }
        }
        
      }else{ 
        hit_time = ndays;
        for(int r=0;r<rho_vec.length();r++){
          hit_table2(r,k) += ndays;
        }
      }
      
      hit_table(0,k) += hit_time;
      //Rcout<<"yes"<< hit_table(1,k)<<std::endl;
      
    }
    if(loop==nloop){
      hit_table = hit_table / nloop;
      hit_utility = hit_total / nloop;
      hit_table2 = hit_table2 / nloop;
      hit_utility2 = hit_total2 / nloop;
    }
}
  
  
  
  return(List::create(Named("hit_table")=hit_table,_["hit_utility"]=hit_utility,
                      Named("hit_table2")=hit_table2,_["hit_utility2"]=hit_utility2));
  //return(List::create(_["hello"]=0));
}
  
