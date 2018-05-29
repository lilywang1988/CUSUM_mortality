#include <Rcpp.h>
//#include <Rcpp/vector.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector O_E_CUSUM_hval(int nloop,int yr_size,NumericVector theta1, NumericVector theta0,double mu,double tau,double yr_er,double p,double yr_int=1,double start_yr=0){
  if(yr_int>tau) stop("Error: your yr_int>tau, this may be invalid!");
  int size=trunc(yr_size*tau);
  int start_day=trunc(start_yr*365)+1;
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
    //srand (loop);
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
List O_E_CUSUM_rho(int nloop,NumericVector h, NumericVector rho_list,int yr_size,NumericVector theta1, NumericVector theta0,double mu,double tau,double yr_er,double yr_int=1,double start_yr=0,double tauL=4.5){
  if(yr_int>tau) warning("Warning: your yr_int>tau");
  if(sum((rho_list>1)|(rho_list<0))!=0) {rho_list=abs(rho_list/max(rho_list));
                                warning("Rho_list has been corrected.");}//correct rho vavalues if invalid
  int size=trunc(yr_size*tau);
  int start_day=trunc(start_yr*365)+1;
  int hit_time;
  
  NumericVector c1=(exp(theta1)-exp(theta0))/(theta1-theta0)-1;
  IntegerVector sign_c1=sign(c1);
  NumericVector abs_c1=abs(c1);
  //NumericVector h_sel(c1.length());
  //int sN=floor(nloop*p);
 // NumericMatrix h_quantiles(sN,c1.length());
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
    //srand (loop);
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
// O_E_CUSUM_rho_t is for is for some spefic rho_t/rho simulations. 
List O_E_CUSUM_rho_t(int nloop,NumericVector h, NumericVector rho_t,int yr_size,NumericVector theta1, NumericVector theta0,double mu,double tau,double yr_er,double yr_int=1,double start_yr=0,double tauL=4.5){
    if(yr_int>tau) warning("Warning: your yr_int>tau");
    if(sum((rho_t>1)|(rho_t<0))) {stop("Error: invalid rho_t.");}
    int size=trunc(yr_size*tau);
    int start_day=trunc(start_yr*365)+1;
    int hit_time;
    
    NumericVector c1=(exp(theta1)-exp(theta0))/(theta1-theta0)-1;
    IntegerVector sign_c1=sign(c1);
    NumericVector abs_c1=abs(c1);
    //NumericVector h_sel(c1.length());
    //int sN=floor(nloop*p);
    //NumericMatrix h_quantiles(sN,c1.length());
    double gamma_subject=-log(1-yr_er)/365;
    int ndays=trunc(tau*365);
    int ndaysL=trunc(tauL*365);
    IntegerVector time_list=seq_len(ndays-start_day+1);
    int rho_length=rho_t.length();
    if(rho_length!=ndays-start_day+2){ //warning("Warning: improper length of rho_t");
      if(rho_length<ndays-start_day+2)
        {double rho_tail=rho_t(rho_length-1);
        for(int l=0;l<ndays-start_day+2-rho_length;l++){rho_t.push_back(rho_tail);}}
        else{rho_t.erase(ndays-start_day+2,rho_length-1);}}
    //print(rho_t);
    //Rcout<<"rho_t.length()="<<rho_t.length()<<"ndays-start_day+2="<<ndays-start_day+2<<std::endl;
    time_list.push_front(0);
    NumericVector M_t_pre;
    NumericVector M_t_pre2;
    NumericVector M_t;
    NumericVector M_t2;
    NumericVector M_t_rho; 
    NumericVector hit_table(c1.length());
    NumericVector hit_utility(c1.length());
    NumericVector hit_total(c1.length());
    NumericVector hit_count(c1.length()); 
    for(int loop=1;loop<=nloop; loop++){
      Rcout<<loop<<std::endl;
      //srand (loop);
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
        M_t2=-M_t_pre+(1-rho_t)*h(k);
        M_t_rho=pmin(M_t,M_t2);
        //Rcout<<"min(M_t_rho)="<<min(M_t_rho)<<' '<<"M_t_rho(0)="<<M_t_rho(0)<<" "<<"M_t(0)"<<M_t(0)<<" "<<"M_t_pre(0)"<<M_t_pre(0)<<" "<<"k="<<k<<std::endl;
        if(min(M_t_rho)<0){
          //Rcout<<"yes"<<std::endl;
          hit_time=min(as<NumericVector>(time_list[M_t_rho<0]));
          //Rcout<<hit_time<<std::endl;
          hit_table(k)+=hit_time;
          hit_total(k)+=(hit_time<=ndays-start_day+1)?1.0:0.0;
          hit_count(k)+=(hit_time<=ndaysL-start_day+1)?1.0:0.0;
          //print(time_list[M_t_rho.row(k)<0]);
          // Rcout<<"hittime="<<hit_time<<"k="<<k<<"r="<<r<<std::endl;
        }
        if(loop==nloop){
          hit_table(k)/=hit_total(k);
          hit_utility(k)=hit_total(k)/nloop;
          hit_count(k)/=nloop;
          //print(hit_table);
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


// [[Rcpp::export]]
// O_E_CUSUM_rho_vs_rho_t is for to compare the performance of rho and rho_t or time_dependent
List O_E_CUSUM_rho_vs_rho_t(int nloop,NumericVector h, NumericVector rho_list,IntegerMatrix rho_t_matrix,int yr_size,NumericVector theta1, NumericVector theta0,double mu,double tau,double yr_er,double yr_int=1,double start_yr=0,double tauL=4.5){
  if(yr_int>tau) warning("Warning: your yr_int>tau");
  if(sum((rho_list>1)|(rho_list<0))!=0) {rho_list=abs(rho_list/max(rho_list));
    warning("Rho_list has been corrected.");}//correct rho vavalues if invalid
  if(rho_list.length()!=rho_t_matrix.nrow()) stop("Error: rho_list.length()!=rho_t_matrix.nrow(). ");
  int size=trunc(yr_size*tau);
  int start_day=trunc(start_yr*365)+1;
  int hit_time;
  
  NumericVector c1=(exp(theta1)-exp(theta0))/(theta1-theta0)-1;
  IntegerVector sign_c1=sign(c1);
  NumericVector abs_c1=abs(c1);
  
  if(c1.length()!=rho_t_matrix.ncol()) stop("c1.length()!=rho_t_matrix.ncol()");
  
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
  NumericVector M_t_rho_t;
  NumericMatrix hit_table(rho_list.length(),c1.length());
  NumericMatrix hit_utility(rho_list.length(),c1.length());
  NumericMatrix hit_total(rho_list.length(),c1.length());
  NumericMatrix hit_count(rho_list.length(),c1.length()); 
  NumericMatrix hit_t_table(rho_list.length(),c1.length());
  NumericMatrix hit_t_utility(rho_list.length(),c1.length());
  NumericMatrix hit_t_total(rho_list.length(),c1.length());
  NumericMatrix hit_t_count(rho_list.length(),c1.length()); 
  NumericVector ARL_rho;
  NumericVector ARL_rho_t;
  IntegerVector ARL_rho_index;
  IntegerVector ARL_rho_t_index;
  for(int loop=1;loop<=nloop; loop++){
    Rcout<<loop<<std::endl;
    //srand (loop);
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
    O_E_t.push_front(0);
    //print(O_t);
    // print(E_t);
    for(int k=0;k<c1.length();k++){
      M_t_pre=sign_c1(k)*O_E_t-abs_c1(k)*E_t;
      M_t_pre2 = NumericVector(cummin(M_t_pre));
      M_t = M_t_pre2-M_t_pre+h(k);
      for(int r=0;r<rho_list.length();r++){
        
        M_t2=(1-rho_list(r))*h(k)-M_t_pre;
        M_t_rho=pmin(M_t,M_t2);
        
        M_t_rho_t= clone(M_t_rho);
        int start=((rho_t_matrix(r,k)+1)<(ndays-start_day+1))?(rho_t_matrix(r,k)+1):(ndays-start_day+1);
        for(int l=start;l<ndays-start_day+2;l++){
          M_t_rho_t(l)=M_t(l);
        }
        
        
        
        //Rcout<<"min(M_t_rho)="<<min(M_t_rho)<<' '<<"M_t_rho(0)="<<M_t_rho(0)<<" "<<"M_t(0)"<<M_t(0)<<" "<<"M_t_pre(0)"<<M_t_pre(0)<<" "<<"k="<<k<<"rho="<<rho_list(r)<<std::endl;
        if(min(M_t_rho)<0){

          // if(k==0) Rcout<<"yes"<< r<<std::endl;
          hit_time=min(as<NumericVector>(time_list[M_t_rho<0]));
          ARL_rho.push_back(hit_time);
          ARL_rho_index.push_back(r);
          //Rcout<<hit_time<<std::endl;
          hit_table(r,k)+=hit_time;
          hit_total(r,k)+=(hit_time<=ndays-start_day+1)?1.0:0.0;
          hit_count(r,k)+=(hit_time<=ndaysL-start_day+1)?1.0:0.0;
        }
        if(min(M_t_rho_t)<0){
          // if(k==0) Rcout<<"yes"<< r<<std::endl;
          hit_time=min(as<NumericVector>(time_list[M_t_rho_t<0]));
          ARL_rho_t.push_back(hit_time);
          ARL_rho_t_index.push_back(r);
          //Rcout<<hit_time<<std::endl;
          hit_t_table(r,k)+=hit_time;
          hit_t_total(r,k)+=(hit_time<=ndays-start_day+1)?1.0:0.0;
          hit_t_count(r,k)+=(hit_time<=ndaysL-start_day+1)?1.0:0.0;
        }
        
        
        //print(time_list[M_t_rho.row(k)<0]);
        // Rcout<<"hittime="<<hit_time<<"k="<<k<<"r="<<r<<std::endl;
        
        if(loop==nloop){
          hit_table(r,k)/=hit_total(r,k);
          hit_utility(r,k)=hit_total(r,k)/nloop;
          hit_count(r,k)/=nloop;
          hit_t_table(r,k)/=hit_t_total(r,k);
          hit_t_utility(r,k)=hit_t_total(r,k)/nloop;
          hit_t_count(r,k)/=nloop;          
          
          //print(hit_table);
        }
      }
      
    }
    // Rcout<<h_temp<<std::endl;
    // print(h_quantiles);
  }
  //print(h_quantiles);
  return(List::create(_["hit_table"]=hit_table, _["hit_count"]=hit_count, _["hit_utility"]=hit_utility,
                      _["hit_t_table"]=hit_t_table, _["hit_t_count"]=hit_t_count, _["hit_t_utility"]=hit_t_utility,
                        _["ARL_rho"]=ARL_rho,_["ARL_rho_index"]=ARL_rho_index,_["ARL_rho_t"]=ARL_rho_t,_["ARL_rho_t_index"]=ARL_rho_t_index));
}

/*** R
#CUSUM_data_gen:generate data for O_E_calc as example
CUSUM_data_gen=function(mu,yr_size,tau,yr_er,yr_int=1,start_yr=1,seed=1,beta=0,change_yr=F,change_rate=0){
  set.seed(seed)
  size=trunc(yr_size*tau)
  gamma_subject=-log(1-yr_er)/365
  ndays=trunc(tau*365)
  start_day=trunc(start_yr*365)+1
  time_list=0:(ndays-start_day+1)
  Lambda0=(0:ndays)*gamma_subject
  enrl_gp=rexp(size*2, yr_size)
  enrl_t=cumsum(enrl_gp)
  enrl_t=floor((enrl_t[enrl_t<tau])*365)
  N3=length(enrl_t)
  pre_time=trunc(rexp(N3,(gamma_subject)*exp(mu)))
  if(change_yr!=F){
    ind_change <- which((pre_time+enrl_t)>(change_yr*365))
    pre_time[ind_change]<- ceiling((pre_time[ind_change]+enrl_t[ind_change]-pmax(change_yr*365,enrl_t)[ind_change])*exp(mu)/exp(change_rate)+pmax(change_yr*365,enrl_t)[ind_change]-enrl_t[ind_change])
  }
  delta_list=((pre_time<=365*yr_int) & ((enrl_t+pre_time)<(tau*365)))
  time=trunc(pmin(pre_time,pmin(365*yr_int,tau*365-enrl_t)));
  cho_time=enrl_t+time;
  x=rnorm(N3)
  xbeta=x*beta
  return(list(time_list=time_list,delta_list=delta_list,cho_time=cho_time,enrl_t=enrl_t,xbeta=xbeta,Lambda0=Lambda0,total_size=N3))
  
}

#O_E_calc: main function for the O-E CUSUM curve mornitoting; Lambda0 must start from time 0, end at least at the end of the followup day (unit:day). 
O_E_CUSUM_calc=function(h, rho_t,restart,delta, enrl_t, cho_time, xbeta,theta1,theta0,Lambda0,tau,yr_int=1,start_yr=0,Hstart=F){
  if(yr_int>tau) warning("Warning: your yr_int>tau")
  if(sum((rho_t>1)|(rho_t<0))) {stop("Error: invalid rho_t.")}
  if(length(delta)!=length(enrl_t)||length(enrl_t)!=length(cho_time)||length(cho_time)!=length(xbeta)){
    stop("Error: invalid input data.")
  } else {N3=length(delta)}

  start_day=trunc(start_yr*365)+1
  c1=(exp(theta1)-exp(theta0))/(theta1-theta0)-1
  sign_c1=sign(c1)
  abs_c1=abs(c1)
  c1_length=length(c1)
  if(sum(Hstart)>0) restart=Hstart|restart  #Hstart has higher priority to restart all
  rho_t=matrix(rho_t)
  if(ncol(rho_t)==1){
    rho_t=matrix(replicate(c1_length,rho_t),ncol=c1_length)
  }
  ndays=trunc(tau*365)
  time_list=0:(ndays-start_day+1)
  rho_nrow=nrow(rho_t)
  if(rho_nrow!=ndays-start_day+2){ #warning("Warning: improper length of rho_t")
    if(rho_nrow<ndays-start_day+2)
    { 
      rho_t=rbind(rho_t,t(replicate(ndays-start_day+2-rho_nrow,rho_t[rho_nrow,])))
    }
    else{
      rho_t=rho_t[1:(ndays-start_day+2),]
    }
  }
  #print(rho_t)
  M_t_pre=NULL
  O_t_pre=E_t_pre=matrix(0,nrow=N3,ncol=ndays-start_day+1);
  for (i in 1:N3){
    if(cho_time[i]>=start_day){
      true_in=max(start_day,enrl_t[i]) 
      for(j in 1:(ndays-start_day+1)){
        if(j>=cho_time[i]-start_day && delta[i]){
          O_t_pre[i,j]=1;
        }
        if(j>=cho_time[i]-start_day) {
          E_t_pre[i,j]=Lambda0[cho_time[i]-true_in+1]*exp(xbeta[i])
        } else if(j>=true_in-start_day && j<=cho_time[i]-start_day){
          E_t_pre[i,j]=Lambda0[j-(true_in-start_day)+1]*exp(xbeta[i])
        }
      }
      
    }
  }

  O_t=c(0,colSums(O_t_pre))
  E_t=c(0,colSums(E_t_pre))
  O_E_t=O_t-E_t
  cross=NULL
  start_times=vector("list",c1_length)
  M_restart=matrix(0,nrow=ndays-start_day+2,ncol=c1_length)
  for(k in 1:c1_length){
    M_t_pre=sign_c1[k]*O_E_t-abs_c1[k]*E_t
    M_min=0
    cross[k]=0
    start_t=1
    for(t in 1:(ndays-start_day+2)){
      M_min=min(M_min,M_t_pre[t]);
      M_pt1=M_min-M_t_pre[t]+h[k];
      
      if((restart[k]&&cross[k]>0)||Hstart)  {
        M_pt2=-(M_t_pre[t]-M_t_pre[start_t+1])+(1-rho_t[max(t-start_t+1,1),k])*h[k]
        M_restart[t,k]=min(M_pt1,M_pt2)
      }else{
        M_restart[t,k]=M_pt1;
      }
      
      
      if(M_restart[t,k]<0){

      if(restart[k]){
        cross[k]=cross[k]+1
        start_times[[k]][cross[k]]=t-1
        start_t=t
      } else if(!restart[k]&&cross[k]==0){
        cross[k]=1
        start_times[[k]][1]=t-1
      }
        
        if(restart[k]==T &&t<ndays-start_day+2){
          M_min=M_t_pre[t+1]
        }
        
      }
    }
  }
  
  return(list(M_restart=M_restart,O_E_t=O_E_t,time_list=time_list,signal_times=start_times,signal.n=cross,theta1=theta1,theta0=theta0,h=h))
}

std_CUSUM_calc<-function(h, rho_t,restart,delta, enrl_t, cho_time, xbeta,theta1,theta0,Lambda0,tau,yr_int=1,start_yr=0,Hstart=F){
  if(yr_int>tau) warning("Warning: your yr_int>tau")
  if(sum((rho_t>1)|(rho_t<0))) {stop("Error: invalid rho_t.")}
  if(length(delta)!=length(enrl_t)||length(enrl_t)!=length(cho_time)||length(cho_time)!=length(xbeta)){
    stop("Error: invalid input data.")
  } else {N3=length(delta)}
  
  L=h*abs(theta1-theta0)
  start_day=trunc(start_yr*365)+1
  c1=(exp(theta1)-exp(theta0))/(theta1-theta0)-1
  sign_c1=sign(c1)
  abs_c1=abs(c1)
  c1_length=length(c1)
  if(sum(Hstart)>0) restart=Hstart|restart  #Hstart has higher priority to restart all
  rho_t=matrix(rho_t)
  if(ncol(rho_t)==1){
     rho_t=matrix(replicate(c1_length,rho_t),ncol=c1_length)
  }
  ndays=trunc(tau*365)
  time_list=0:(ndays-start_day+1)
  rho_nrow=nrow(rho_t)
    if(rho_nrow!=ndays-start_day+2){ #warning("Warning: improper length of rho_t")
      if(rho_nrow<ndays-start_day+2)
      { 
        rho_t=rbind(rho_t,t(replicate(ndays-start_day+2-rho_nrow,rho_t[rho_nrow,])))
      }
      else{
        rho_t=rho_t[1:(ndays-start_day+2),]
      }
    }
#print(rho_t)
    O_t_pre=E_t_pre=matrix(0,nrow=N3,ncol=ndays-start_day+1)
    for (i in 1:N3){
      if(cho_time[i]>=start_day){
        true_in=max(start_day,enrl_t[i]) 
        for(j in 1:(ndays-start_day+1)){
          if(j>=cho_time[i]-start_day && delta[i]){
             O_t_pre[i,j]=1;
          }
          if(j>=cho_time[i]-start_day) {
             E_t_pre[i,j]=Lambda0[cho_time[i]-true_in+1]*exp(xbeta[i])
          } else if(j>=true_in-start_day && j<=cho_time[i]-start_day){
            E_t_pre[i,j]=Lambda0[j-(true_in-start_day)+1]*exp(xbeta[i])
          }
        }
          
      }
    }
      
    O_t=c(0,colSums(O_t_pre))
    E_t=c(0,colSums(E_t_pre))
    cross=NULL
    start_times=vector("list",c1_length)
    S_restart=matrix(0,nrow=ndays-start_day+2,ncol=c1_length)
    for(k in 1:c1_length){
      R_t=(theta1[k]-theta0[k])*O_t-(exp(theta1[k])-exp(theta0[k]))*E_t
      R_min=0
      cross[k]=0
      start_t=1
      for(t in 1:(ndays-start_day+2)){
        R_min=min(R_min,R_t[t])
        R_pt1=R_t[t]-R_min
        if((restart[k]&&cross[k]>0)||Hstart)  {
          R_pt2=R_t[t]-R_t[start_t+1]+L[k]*rho_t[t-start_t+1,k]
          S_restart[t,k]=max(R_pt1,R_pt2)
        }else{
          S_restart[t,k]=R_pt1
        }
            
        if(S_restart[t,k]>L[k]){
              
          if(restart[k]){
            cross[k]=cross[k]+1
            start_times[[k]][cross[k]]=t-1
            start_t=t
                
          } else if(!restart[k]&&cross[k]==0){
            cross[k]=1
            start_times[[k]][1]=t-1
          }
              
          if(restart[k]==T && t<ndays-start_day+2){
             R_min=R_t[t+1]
          }
        }
      }
    }
        
      return(list(S_restart=S_restart,time_list=time_list,signal_times=start_times,signal.n=cross,theta1=theta1,theta0=theta0,L=L))
        
}

#plot the CUSUM curves in R ggplot2
CUSUM_plot=function(result,O_E=T,adjust=T){ 
  require(ggplot2)
  require(reshape)
  theta1=result$theta1
  theta0=result$theta0
  c1_length=length(theta1)
  name_temp1=round(exp(theta1),digits=2)
  name_temp0=round(exp(theta0),digits=2)
  sign0=sign1=rep("",c1_length)
  for(d in 1:c1_length){
    if(name_temp1[d]<1)  {name_temp1[d]=round(1/name_temp1[d],2); sign1[d]="-"}
    if(name_temp0[d]<1)  {name_temp0[d]=round(1/name_temp0[d],2); sign0[d]="-"}
  }
  legend_names<-paste0(sign0,"log",name_temp0," vs ", sign1,"log",name_temp1)
  legend_title<-"H0 vs H1: "
  if(O_E) {
     names=c("O_E_t",paste0("test",1:c1_length))
     cross_t<--1000
     indicator_x=names[1]
     for(j in 1:c1_length) {
       cross_t<-c(cross_t,result$signal_times[[j]])
       indicator_x<-c(indicator_x,rep(names[j+1],result$signal.n[j]))
      }
     cross_data<-data.frame(cross_t,indicator_x=factor(indicator_x))
     if(adjust==T)M_t_adj=result$M_restart%*%diag(theta1-theta0) else M_t_adj=(result$M_restart)%*%diag(sign(theta1))
     M_t_preplot=(M_t_adj+as.matrix(replicate(c1_length,result$O_E_t)))
     colnames(M_t_preplot)<-paste0("test",1:c1_length)                               
     M_t_data=data.frame(time_list=result$time_list,O_E_t=result$O_E_t,M_t_preplot)
     M_t_plot=melt( M_t_data,id="time_list")
     M_t_plot$variable=as.factor(M_t_plot$variable)
     result_plt<-ggplot(data=M_t_plot,aes(x=time_list,y=value,color=variable))+geom_line(size=0.3)+  scale_linetype_manual(values=c("solid",rep("dashed",c1_length)))+
       xlab("time")+ylab("O-E")+ 
       geom_vline(data=cross_data,aes(xintercept=cross_t,color=indicator_x),linetype="dashed",size=0.5)+xlim(c(0,length(result$time_list)*1.005))+theme_light()+theme(plot.title = element_text(hjust = 0.5,size=15))+
       scale_color_manual(name  =legend_title,values=1:(c1_length+1),labels=c("O-E",  legend_names))+ggtitle("O-E CUSUM")
     result_plt
       #+theme_light()
  }else{
    names=paste0("test",1:c1_length)
    cross_t=indicator_x=NULL
    for(j in 1:c1_length) {
      cross_t<-c(cross_t,result$signal_times[[j]])
      indicator_x<-c(indicator_x,rep(names[j],result$signal.n[j]))
    }
    cross_data<-data.frame(cross_t,variable=factor(indicator_x))
    S_preplot=result$S_restart
    colnames(S_preplot)<-paste0("test",1:c1_length)                               
    S_data=data.frame(time_list=result$time_list,S_preplot)
    S_plot=melt( S_data,id="time_list")
    S_plot$variable=as.factor(S_plot$variable)
    L_bar<-data.frame(variable=factor(names),result$L)
    result_plt<-ggplot(data=S_plot,aes(x=time_list,y=value))+geom_line(size=0.3)+facet_wrap( ~ variable, ncol=2,scales="free",labeller = )+
      ggtitle("Standard CUSUM")+theme_bw()+theme(plot.title = element_text(hjust = 0.5,size=15))+
      geom_hline(data = L_bar, aes(yintercept = result.L,color=variable))+geom_vline(aes(xintercept=cross_t,color=variable),data=cross_data,linetype="dashed",size=0.5)+
      scale_color_manual(name  = legend_title,values=2:(c1_length+1),labels=  legend_names)
    result_plt
    
  }
}


if(F){

options(max.print=1000000)

  theta1=c(log(2),-log(2))
  theta0=c(log(1),-log(1))
  mu=log(1)
  start_yr=1
  yr_er=0.1
  tau=4.5
  h=c(6.463784009,5.602556288)
  #data=CUSUM_data_gen(mu,100,tau,yr_er)
  data=CUSUM_data_gen(mu,100,tau,yr_er,change_yr=2,change_rate=log(0.5))
#either scalor, two-length vector or long tables for input of rho_t
  result1=O_E_CUSUM_calc(h,0.5,c(T,T),data$delta_list, data$enrl_t, data$cho_time, data$xbeta,theta1,theta0,data$Lambda0,tau,yr_int=1,start_yr=1)
  result2=std_CUSUM_calc(h,0.5,c(T,T),data$delta_list, data$enrl_t, data$cho_time, data$xbeta,theta1,theta0,data$Lambda0,tau,yr_int=1,start_yr=1)
  CUSUM_plot(result2,F)
theta1=c(log(2),-log(2))
theta0=c(log(1),-log(1))
mu=log(2)
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
O_E_CUSUM_rho(1000,h,seq(0,1,0.1),100,c(log(2),-log(2)),c(0,0),mu,10,0.1,1,1,4.5)
O_E_CUSUM_rho(100,h,0.5,100,c(log(2),-log(2)),c(0,0),mu,10,0.1,1,1,4.5)
  ARL=2
  rho1_days=trunc(ARL*365)
  rho2_days=trunc(tau*365)
mu=log(2)
  res1=O_E_CUSUM_rho_t(1000,h,0.5,100,c(log(2),-log(2)),c(0,0),mu,10,0.1,1,1,4.5)
  res2=O_E_CUSUM_rho_t(1000,h,c(rep(0.5,230),rep(0,tau*365)),100,c(log(2),-log(2)),c(0,0),mu,10,0.1,1,1,4.5)
  res3=O_E_CUSUM_rho_t(1000,h,0,100,c(log(2),-log(2)),c(0,0),mu,10,0.1,1,1,4.5)
mu=log(1)
  res4=O_E_CUSUM_rho_t(1000,h,0.5,100,c(log(2),-log(2)),c(0,0),mu,10,0.1,1,1,4.5)
  res5=O_E_CUSUM_rho_t(1000,h,c(rep(0.5,230),rep(0,tau*365)),100,c(log(2),-log(2)),c(0,0),mu,10,0.1,1,1,4.5)
  res6=O_E_CUSUM_rho_t(1000,h,0,100,c(log(2),-log(2)),c(0,0),mu,10,0.1,1,1,4.5)

  theta1=c(log(2),-log(2))
  theta0=c(log(1),-log(1))
  mu=log(2)
  start_yr=1
  yr_er=0.1
  tau=4.5
  p=0.08
  h=c(6.463784009,5.602556288)
  res_rho_vs_rho_t<-O_E_CUSUM_rho_vs_rho_t(1000,h[1],seq(0,0.9,0.1),as.matrix(rep(224.658/2,10)),100,log(2),c(0),mu,10,0.1,1,1,4.5)
  res_rho_vs_rho_t$hit_table
  res_rho_vs_rho_t$hit_t_table
  
  res_rho_vs_rho_t$hit_count
  res_rho_vs_rho_t$hit_t_count
}


*/
