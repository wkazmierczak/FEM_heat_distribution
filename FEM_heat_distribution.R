n <- 7

# base functions starting from 0, because it was easier for me to track indexes
e <- function(i, x){
  interval_width = 2/n
  if (x > interval_width*(i-1) && x <= interval_width*(i)){
    return (n/2*x - i + 1)
  }
  else if (x > interval_width*(i) && x < interval_width*(i+1)){
    return (-n/2*x + i + 1)
  }
  else {
    return (0)
  }
}

e_derivative <- function(i, x){
  interval_width = 2/n
  if (x > interval_width*(i-1) && x <= interval_width*(i)){
    return (n/2)
  }
  else if (x > interval_width*(i) && x < interval_width*(i+1)){
    return (-n/2)
  }
  else {
    return (0)
  }
}

# draw plot for base functions
x_values <- seq(0, 2, length.out = 1000)

plot(x_values, sapply(x_values, function(x) e(0, x)), type = "l", col = 100, lty = 1,
     xlab = "x", ylab = "e(i, x)", main = "Plot of the function e(i, x)")

for (i in 1:n) {
  y_values <- sapply(x_values, function(x) e(i, x))
  lines(x_values, y_values, col = i, lty = i)
}

legend("topright", legend = paste("i =", 0:n), col = 0:n, lty = 0:n, cex = 0.8)


# product of derivatives of the base functions
integrand_k_u_v <- function(i, j){
  return (function(x){
      if (x>=0 && x<= 1){
        return(e_derivative(i, x)*e_derivative(j, x))
      }
      else if(x > 1 && x<= 2){
        return(2*e_derivative(i, x)*e_derivative(j, x))
      }
      else{
        return (0)
      }
    }
  )
}

gaussian_legendre_2_points_interval <- function(f, a, b) {
  standard_points <- c(-sqrt(3)/3, sqrt(3)/3)
  weights <- c(1, 1)
  
  points_suited_for_interval <- 0.5 * (b - a) * standard_points + 0.5 * (a + b)
  
  result <- 0.5 * (b - a) * sum(weights[1] * f(points_suited_for_interval[1]), weights[2] * f(points_suited_for_interval[2]))
  
  return(result)
}

# count B(e_i, e_j)
B <- function(i, j){
  interval_width = 2/n
  
  a = interval_width * max(0, i-1, j-1)
  b = interval_width * min(i+1, j+1)
  
  return(gaussian_legendre_2_points_interval(integrand_k_u_v(i,j), a, b) - e(i, 0)*e(j, 0))
}

# count L(e_j)
L <- function(j){
  return((-17)*e(j, 0))
}

matrix_equation_res <- function(){
    M <- matrix(0, n, n)
    
    for (i in 1:(n)){
      for (j in 1:(n)){
        M[i,j] <- B(j-1, i-1) # B-indexes because of indexes of base functions
      }
    }
    
    C <- sapply(0:(n-1), L)
    
    return(solve(M, C))
}

lin_combination_of_base_fun <- function(x, vector){
  return(sum(vector * sapply(0:(n-1), e, x = x))) # multiply i-elem of vector e(i-1, x), because of indexing
  }


FEM_heat_distribution <- function(){
  w_vector = matrix_equation_res()
  
  approximated_function <- function(x){
    return (lin_combination_of_base_fun(x, w_vector))
  }

  # draw plot
  x_values <- seq(0, 2, length.out = 1000)
  y_values <- sapply(x_values, approximated_function)
  
  plot(x_values, y_values, type = "l", col = "blue", lwd = 2, xlab = "x", ylab = "y", main = "Approximated heat distribution function")
  
  return (approximated_function)
}


FEM_heat_distribution()
