a <- read.table("http://www.stat.ucla.edu/~nchristo/statistics100C/jura.txt", header=TRUE)
y <- a$Pb
x1 <- a$Cd
x2 <- a$Cu
x3 <- a$Zn
x4 <- a$Cr

one_vector <-rep(1, 359)
x <- cbind(one_vector, x1,x2,x3,x4)
x_p <- t(x)
x_x_p <- x_p %*% x
x_x_p_inverse <- solve(t(x) %*% x)
xp_y <- x_p %*% y
H <- x%*%x_x_p_inverse%*%x_p
beta_h <- x_x_p_inverse %*% xp_y

#part(i)
x1 <- cbind(one_vector, x1)
x2 <- cbind(x2, x3,x4)
x1_p <- t(x1)
x1_x1_p <- x1_p %*% x1
x1_x1_p_inverse <- solve(t(x1) %*% x1)
H1 <- x1%*%x1_x1_p_inverse%*%x1_p

x2_star <- x2 - H1%*%x2
y_star <- y - H1%*%y

beta_21 <- solve(t(x2_star)%*%(x2_star))%*%t(x2_star)%*%y_star
beta_12 <- solve(t(x1)%*%x1)%*%t(x1)%*%y-solve(t(x1)%*%x1)%*%t(x1)%*%x2%*%beta_21
beta_partial <- rbind(beta_12,beta_21)
library("matrixcalc")
is.idempotent.matrix(H-H1, tol = 1e-08)

#part(j)
c <- t(rbind(c(1,2),c(2,-1),c(-3,1),c(5,3),c(-1,0)))
g <- t(rbind(c(5,10)))

beta_lagrange <- beta_h - x_x_p_inverse %*% t(c) %*% solve(c %*% x_x_p_inverse %*% t(c)) %*% (c %*% beta_h - g)

c1 <- t(rbind(c(1,2),c(2,-1)))
c2 <- t(rbind(c(-3,1),c(5,3),c(-1,0)))

yr <- y - x1 %*% solve(c1)%*%g
x2r <- x2 - x1%*%solve(c1)%*%c2
beta_2_c <- solve(t(x2r)%*%(x2r))%*%t(x2r)%*%yr
beta_1_c <- solve(c1)%*%(g-c2%*%beta_2_c)
beta_cononical <- rbind(beta_1_c, beta_2_c)
      