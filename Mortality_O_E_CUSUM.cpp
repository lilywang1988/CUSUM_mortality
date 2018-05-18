#include <Rcpp.h>
#include <R.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector O_E_CUSUM_hval(int nloop,int yr_size,NumericVector theta1, NumericVector theta0,double mu,double tau,double yr_er,double p,double yr_int=1,double start_yr=0){
  if(yr_int>tau) stop("Error: your yr_int>tau, this may be invalid!");
  int size=trunc(yr_size*tau);
  int start_day=start_yr*365+1;
  NumericVector c1=(exp(theta1)-exp(theta0))/(theta1-theta0)-1;
  IntegerVector sign_c1=sign(c1);
  NumericVector abs_c1=abs(c1);
  NumericVector h_sel(c1.length());
  double h_temp;
  int sN=floor(nloop*p);
  NumericMatrix h_quantiles(sN,c1.length());
  double gamma_subject=-log(1-yr_er)/365;
  int ndays=trunc(tau*365);
  IntegerVector time_list=seq_len(ndays-start_day+1);
  time_list.push_front(0);
  for(int loop=1;loop<=nloop; loop++){
    Rcout<<loop<<std::endl;
    srand (loop);
    NumericVector enrl_gp=rexp(size*2, yr_size);
    NumericVector enrl_t=cumsum(enrl_gp);
    enrl_t=floor(as<NumericVector>(enrl_t[enrl_t<tau])*365);
    int N3=enrl_t.length();
    //Rcout<<N3<<std::endl;
    NumericVector pre_time=trunc(rexp(N3,(gamma_subject)*exp(mu)));
    LogicalVector delta_list=((pre_time<=365*yr_int) & ((enrl_t+pre_time)<(tau*365)));
    NumericVector time=trunc(pmin(pre_time,pmin(365*yr_int,tau*365-enrl_t)));
    NumericVector cho_time=enrl_t+time;
    NumericMatrix E_t_pre(N3,ndays-start_day+1);
    NumericMatrix O_t_pre(N3,ndays-start_day+1);
    NumericVector M_t_pre(ndays-start_day+1);
    NumericVector M_t_pre2(ndays-start_day+1);
    for (int i=0;i<N3;i++){
      if(cho_time[i]>=start_day){
        //int true_fu=(start_day>enrl_t[i])?(time[i]-(start_day-enrl_t[i])):time[i];
        int true_in=(start_day>enrl_t[i])?start_day:enrl_t[i];//enrl_t[i]
        //int true_in=enrl_t[i];
        for(int j=0;j<ndays-start_day+1;j++){
          if(j>=cho_time[i]-start_day && delta_list[i]){
            O_t_pre(i,j)=1;
          }
          if(j>=cho_time[i]-start_day) {
            E_t_pre(i,j)=gamma_subject*(cho_time[i]-true_in);// start_day<=cho_time && enrl_time< cho_Time
          } else if(j>=true_in-start_day && j<=cho_time[i]-start_day){
            E_t_pre(i,j)=gamma_subject*(j-(true_in-start_day));
          }
        }
      
      }
     }
      NumericVector O_t(ndays-start_day+1);
      NumericVector E_t (ndays-start_day+1);
      NumericVector O_E_t (ndays-start_day+1);
      for(int j=0; j<ndays-start_day+1; j++){
        O_t[j]=sum(O_t_pre.column(j));
        E_t[j]=sum(E_t_pre.column(j));
      }
      O_E_t=O_t-E_t;
      //print(O_t);
      //print(E_t);
    for(int k=0;k<c1.length();k++){
        M_t_pre=sign_c1[k]*O_E_t-abs_c1[k]*E_t;
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
 // return(0);
}
// [[Rcpp::export]]
// O_E_CUSUM_rho is for performance check on different rho values: ARL and percent of hits
List O_E_CUSUM_rho(int nloop,NumericVector h, NumericVector rho_list,int yr_size,NumericVector theta1, NumericVector theta0,double mu,double tau,double yr_er,double p,double yr_int=1,double start_yr=0,double tauL=4.5){
  if(yr_int>tau) stop("Error: your yr_int>tau, this may be invalid.");
  if(sum((rho_list>1)|(rho_list<0))!=0) {rho_list=abs(rho_list/max(rho_list));
                                warning("Rho_list has been corrected.");}//correct rho vavalues if invalid
  int size=trunc(yr_size*tau);
  int start_day=start_yr*365+1;
  int hit_time;
  
  NumericVector c1=(exp(theta1)-exp(theta0))/(theta1-theta0)-1;
  IntegerVector sign_c1=sign(c1);
  NumericVector abs_c1=abs(c1);
  NumericVector h_sel(c1.length());
  int sN=floor(nloop*p);
  NumericMatrix h_quantiles(sN,c1.length());
  double gamma_subject=-log(1-yr_er)/365;
  int ndays=trunc(tau*365);
  int ndaysL=trunc(tauL*365);
  IntegerVector time_list=seq_len(ndays-start_day+1);
  time_list.push_front(0);
  NumericVector M_t_pre;
  NumericVector M_t_pre2;
  NumericVector M_t;
  NumericVector M_t2;
  NumericVector M_t_rho; 
  NumericMatrix hit_table(rho_list.length(),c1.length());
  NumericMatrix hit_utility(rho_list.length(),c1.length());
  NumericMatrix hit_total(rho_list.length(),c1.length());
  NumericMatrix hit_count(rho_list.length(),c1.length()); 
  for(int loop=1;loop<=nloop; loop++){
    Rcout<<loop<<std::endl;
    srand (loop);
    NumericVector enrl_gp=rexp(size*2, yr_size);
    NumericVector enrl_t=cumsum(enrl_gp);
    enrl_t=floor(as<NumericVector>(enrl_t[enrl_t<tau])*365);
    int N3=enrl_t.length();
    //Rcout<<N3<<std::endl;
    NumericVector pre_time=trunc(rexp(N3,(gamma_subject)*exp(mu)));
    LogicalVector delta_list=((pre_time<=365*yr_int) & ((enrl_t+pre_time)<(tau*365)));
    NumericVector time=trunc(pmin(pre_time,pmin(365*yr_int,tau*365-enrl_t)));
    NumericVector cho_time=enrl_t+time;
    NumericMatrix E_t_pre(N3,ndays-start_day+1);
    NumericMatrix O_t_pre(N3,ndays-start_day+1);
    for (int i=0;i<N3;i++){
      if(cho_time[i]>=start_day){
        //int true_fu=(start_day>enrl_t[i])?(time[i]-(start_day-enrl_t[i])):time[i];
        int true_in=(start_day>enrl_t[i])?start_day:enrl_t[i];//enrl_t[i]
        //int true_in=enrl_t[i];
        for(int j=0;j<ndays-start_day+1;j++){
          if(j>=cho_time[i]-start_day && delta_list[i]){
            O_t_pre(i,j)=1;
          }
          if(j>=cho_time[i]-start_day) {
            E_t_pre(i,j)=gamma_subject*(cho_time[i]-true_in);// start_day<=cho_time && enrl_time< cho_Time
          } else if(j>=true_in-start_day && j<=cho_time[i]-start_day){
            E_t_pre(i,j)=gamma_subject*(j-(true_in-start_day));
          }
        }
        
      }
    }
    NumericVector O_t(ndays-start_day+1);
    NumericVector E_t (ndays-start_day+1);
    NumericVector O_E_t (ndays-start_day+1);
    for(int j=0; j<ndays-start_day+1; j++){
      O_t[j]=sum(O_t_pre.column(j));
      E_t[j]=sum(E_t_pre.column(j));
    }
    O_E_t=O_t-E_t;
    //print(O_t);
   // print(E_t);
    for(int k=0;k<c1.length();k++){
      M_t_pre=sign_c1(k)*O_E_t-abs_c1(k)*E_t;
      //Rcout<<"M_t_pre.length()="<<M_t_pre.length()<<std::endl;
      M_t_pre.push_front(0);
      //Rcout<<"M_t_pre.length()="<<M_t_pre.length()<<std::endl;
      M_t_pre2 = NumericVector(cummin(M_t_pre));
      M_t = M_t_pre2-M_t_pre+h(k);
      for(int r=0;r<rho_list.length();r++){

          M_t2=(1-rho_list(r))*h(k)-M_t_pre;
          M_t_rho=pmin(M_t,M_t2);
          //Rcout<<"min(M_t_rho)="<<min(M_t_rho)<<' '<<"M_t_rho(0)="<<M_t_rho(0)<<" "<<"M_t(0)"<<M_t(0)<<" "<<"M_t_pre(0)"<<M_t_pre(0)<<" "<<"k="<<k<<"rho="<<rho_list(r)<<std::endl;
         if(min(M_t_rho)<0){
          // if(k==0) Rcout<<"yes"<< r<<std::endl;
          hit_time=min(as<NumericVector>(time_list[M_t_rho<0]));
           //Rcout<<hit_time<<std::endl;
          hit_table(r,k)+=hit_time;
          hit_total(r,k)+=(hit_time<=ndays-start_day+1)?1.0:0.0;
          hit_count(r,k)+=(hit_time<=ndaysL-start_day+1)?1.0:0.0;
          };
       //print(time_list[M_t_rho.row(k)<0]);
       // Rcout<<"hittime="<<hit_time<<"k="<<k<<"r="<<r<<std::endl;
        
        if(loop==nloop){
          hit_table(r,k)/=hit_total(r,k);
          hit_utility(r,k)=hit_total(r,k)/nloop;
          hit_count(r,k)/=nloop;
          //print(hit_table);
        }
      }

     }
      // Rcout<<h_temp<<std::endl;
      // print(h_quantiles);
    }
  //print(h_quantiles);
  return(List::create(Named("hit_table")=hit_table, _["hit_count"]=hit_count, _["hit_utility"]=hit_utility));
}
  

// [[Rcpp::export]]
// The quantile function to keep the first interesting part of the data
NumericVector Quantile(NumericVector newvalues, NumericVector oldvalues,int sN){
  if(sN!=oldvalues.length()) stop("Error: invalid number of oldvalues for input of the Quantile()");
  for(NumericVector::iterator value=newvalues.begin(); value!=newvalues.end(); value++){
    bool mainloop=false;
    for(NumericVector::iterator ovalue=oldvalues.begin(); ovalue!=newvalues.end() && mainloop==false; ovalue++){
      if(*value>*ovalue){oldvalues.insert(ovalue,*value);oldvalues.erase(oldvalues.end()-1);mainloop=true;}
    }
  }
  return(oldvalues);
}




/*** R
if(F){

options(max.print=1000000)

theta1=c(log(2),-log(2))
theta0=c(log(1),-log(1))
mu=0
start_yr=1
yr_er=0.1
tau=4.5
p=0.08
O_E_CUSUM_hval(1000,100,theta1,theta0,mu,tau,yr_er,p,1,start_yr)
#List O_E_CUSUM_rho(int nloop,NumericVector h, 
#NumericVector rho_list,int yr_size,NumericVector theta1, 
#NumericVector theta0,double mu,double tau,double yr_er,double p,double yr_int=1,
#double start_yr=0,double tauL=4.5)
h=c(6.463784009,5.602556288)
#h=c(6.463784009)
#mu=log(2)
O_E_CUSUM_rho(1000,h,seq(0,1,0.1),100,c(log(2),-log(2)),c(0,0),mu,10,0.1,0.08,1,1,4.5)
}
*/
