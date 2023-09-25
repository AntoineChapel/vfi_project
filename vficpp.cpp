#include <iostream>
#include <armadillo>
#include <cmath>
#include <chrono>



int main(){
  
  auto started = std::chrono::high_resolution_clock::now();
  //If I add a comment, what happens ?
  const int A = 10;
  const float alpha = 0.5;
  const float beta = 0.9;
  const float tol = 1e-6;
  const int maxiter = 1000;

  float kmin = 1;
  float kmax = 25;
  int prec = 1000;

  arma::mat kgrid(1, prec, arma::fill::randu);
  arma::mat gk(1, prec, arma::fill::randu);
  arma::mat Vk0(1, prec, arma::fill::ones);
  arma::mat Vk(1, prec, arma::fill::randu);
  arma::mat Vkprim(1, prec, arma::fill::randu);
  arma::mat value_array(prec, prec, arma::fill::randu);

  int n_iter = 0;
  float k, kprim, c, kstar;
  float norm = 1e5;
  
  kgrid = arma::linspace(kmin, kmax, prec);
  gk = arma::linspace(kmin, kmax, prec);


  Vk = Vk0;

  while (n_iter < maxiter && norm > tol){
    for (int iprim = 0; iprim < prec; iprim++){
      kprim = kgrid(iprim);
      for (int i = 0; i < prec; i++){
        k = kgrid(i);
        c = A*pow(k, alpha) - kprim;
        if (c > 0){
          value_array(i, iprim) = log(c) + beta*Vk(iprim);
        }
        else {
          value_array(i, iprim) = -1e6;
        }
      }
    }
 

    for (int rownum = 0; rownum < prec; rownum++){
      gk(rownum) = kgrid(value_array.row(rownum).index_max());
      Vkprim(rownum) = value_array.row(rownum).max();
    }

    norm = abs(Vkprim - Vk).max();
    Vk = Vkprim;
    n_iter++;
    std::cout << "Iteration " << n_iter << " norm: " << norm << std::endl;
  }
  
  kstar = kgrid(abs(gk - kgrid).index_min());
  std::cout << "The steady state value of capital is: " << kstar << std::endl;

  auto done = std::chrono::high_resolution_clock::now();
  std::cout << "VFI in C++ executed in " << std::chrono::duration_cast<std::chrono::milliseconds>(done-started).count() << " milliseconds" << std::endl;
  return 0;
}
