#include <iostream>
#include <cmath>
#include <vector>
#include <chrono>
#include <functional>


std::vector<double> linspace(double start, double end, int n);
void print_matrix(std::vector<std::vector<double>> matrix);
void print_vector(std::vector<double> vec);
std::vector<double> vector2_ope(std::vector<double> a, std::vector<double> b, std::string ope);
std::vector<double> vector1_ope(std::vector<double> a, std::function<double(double)> func);
double max(std::vector<double> vec);
int argmax(std::vector<double> vec);
int argmin(std::vector<double> vec);
double abs_f(double x);


int main(){
  auto started = std::chrono::high_resolution_clock::now();
  const int A = 10;
  const double alpha = 0.5;
  const double beta = 0.9;
  
  double tol = 1e-6;
  int maxiter = 1000;

  double kmin = 1.0d;
  double kmax = 25.0d;
  int prec = 500;
  static std::vector<double> kgrid(prec);
  std::vector<double> gk(prec, 1);
  std::vector<double> Vk(prec, 1);
  std::vector<double> Vkprim(prec);
  std::vector<double> zeros_init(prec, 0);
  std::vector<std::vector<double>> value_array(prec, zeros_init);
  std::vector<double> value_array_rowiprim(prec);

  int n_iter = 0;
  double k, kprim, c, kstar;
  double norm = 1e5;

  kgrid = linspace(kmin, kmax, prec);

  while (n_iter < maxiter && norm > tol){
    for (int iprim = 0; iprim < prec; iprim++){
      kprim = kgrid[iprim];
      for (int i = 0; i < prec; i++){
        k = kgrid[i];
        c = A*pow(k, alpha) - kprim;
        if (c > 0){
          value_array[i][iprim] = (std::log(c) + beta*Vk[iprim]);
        }
        else {
          value_array[i][iprim] = -1e6;
        }
      }
    }

    for (int rownum = 0; rownum < prec; rownum++){
      int argmax_k = argmax(value_array[rownum]); 
      gk[rownum] = kgrid[argmax_k];
      Vkprim[rownum] = max(value_array[rownum]);
    }
    
    std::vector<double> update, absupdate;
    update = vector2_ope(Vkprim, Vk, "-");
    absupdate = vector1_ope(update, &abs_f);
    norm = max(absupdate);
    Vk = Vkprim;
    n_iter++;
    std::cout << "Iteration " << n_iter << " norm: " << norm << std::endl;
  }
  
  std::vector<double> dist_gk_kgrid, abs_dist_gk_kgrid;
  dist_gk_kgrid = vector2_ope(gk, kgrid, "-");
  
  abs_dist_gk_kgrid = vector1_ope(dist_gk_kgrid, &abs_f);
  int ind_kopt = argmin(abs_dist_gk_kgrid);
  kstar = kgrid[ind_kopt]; 
  std::cout << "The steady state value of capital is: " << kstar << std::endl;
  auto done = std::chrono::high_resolution_clock::now();
  std::cout << "VFI in C++ executed in " << std::chrono::duration_cast<std::chrono::milliseconds>(done-started).count() << " milliseconds" << std::endl;
return 0;

}

std::vector<double> linspace(double start, double end, int n){
  std::vector<double> vec(n);
  for (int i = 0; i < n; i++){
    vec[i] = start + i*(end - start)/(n-1);
  }
  return vec;
}


std::vector<double> vector1_ope(std::vector<double> a, std::function<double(double)> func){
  int size_a = size(a);
  std::vector<double> c(size_a);
  for (int i = 0; i < size_a; i++){
    c[i] = func(a[i]);
  }
  return c;
}

std::vector<double> vector2_ope(std::vector<double> a, std::vector<double> b, std::string ope){
  int size_a = size(a);
  int size_b = size(b); 
  std::vector<double> c(size_a);
  if (size_a != size_b){
    std::cerr << "WARNING: a and b are not of the same size.\n";
  }
  else {
    for (int i = 0; i < size_a; i++){
      if (ope=="+"){
        c[i] = a[i] + b[i];
      }
      else if (ope=="-"){
        c[i] = a[i] - b[i];
      }
    }
  }
  return c;
}


double abs_f(double x){
  return std::abs(x);
}

double max(std::vector<double> vec){
  int size_vec = size(vec);
  double maxval = vec[0];
  for (int i=1; i < size_vec; i++){
    if (vec[i] > maxval){
      maxval = vec[i];
    }
  }
  return maxval;
}

int argmax(std::vector<double> vec){
  int size_vec = size(vec);
  int maxindex = 0;
  double maxval = vec[maxindex];
  for (int i=1; i < size_vec; i++){
    if (vec[i] > maxval){
      maxval = vec[i];
      maxindex = i;
    }
  }
  return maxindex;
}

int argmin(std::vector<double> vec){
  int size_vec = size(vec);
  int minindex = 0;
  double minval = vec[minindex];
  for (int i=1; i < size_vec; i++){
    if (vec[i] < minval){
      minval = vec[i];
      minindex = i;
    }
  }
  return minindex;
}


/* linear algebra printing functions
void print_matrix(std::vector<std::vector<double>> matrix){
  for (std::vector<double> row : matrix){
    for (double val : row){
      std::cout << val << " ";
    }
    std::cout << " " << std::endl;
  }
}

void print_vector(std::vector<double> vec){
  for (double val : vec){
    std::cout << val << " ";
  }
  std::cout << " " << std::endl;
}

*/
