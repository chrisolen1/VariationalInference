library(extraDistr)
library(ggplot2)

sigma = 3
K = 5
n = 1000

set.seed(247)
mu = rnorm(K, mean=0, sd=sigma)
cs = rcat(n, rep(1/K, K))
x = rnorm(n, mean=mu[cs], sd=1)

df = data.frame(x=x, mu=as.factor(cs))
ggplot(df, aes(x=x, color=mu, fill=mu)) + geom_histogram(alpha=0.5)

mk = rnorm(K)
sk2 = rgamma(K,5)
phis = rdirichlet(n, c(1,1,1,1,1))

ELBO = function(mk,  sk2, phis) {
  t = sk2+mk^2
  a = -(1/(2*sigma^2))*sum(t)
  b = 2*sum(sweep(sweep(phis, MARGIN=2, mk, "*"), MARGIN=1, x, "*"))-0.5*sum(sweep(phis, MARGIN=2, t, "*"))
  c = -0.5*sum(log(2*pi*sk2))
  d = sum(phis*log(phis))
  return(a+b+c+d)
}

iter = 30
elbos = rep(NA, iter+1)
elbos[1] = ELBO(mk,sk2,phis)
for (i in 1:iter) {
  phis.new = matrix(nrow=n, ncol=K)
  for (j in 1:n) {
    phis.new[j,] = exp(x[j]*mk-0.5*(sk2+mk^2))
    phis.new = phis.new/rowSums(phis.new)
  }
  phis = phis.new
  
  mk.new = rep(NA, K)
  sk2.new = rep(NA, K)
  for (k in 1:K) {
    sk2.new[k] = 1/(1/sigma^2+sum(phis[,k]))
    mk.new[k] = sk2.new[k]*sum(phis[,k]*x)
  }
  sk2 = sk2.new
  mk = mk.new
  
  elbos[i+1] = ELBO(mk,sk2,phis)
  cat("Iteration: ", i, "ELBO-diff: ", abs(elbos[i+1]-elbos[i]), "\n")
  if (abs(elbos[i+1])<0.1) break
}

ggplot(df, aes(x=x, color=mu, fill=mu)) +
  geom_histogram(alpha=0.5) +
  geom_vline(data=data.frame(x=mk), 
             aes(xintercept=x, 
                 color=as.factor(c(2,4,1,3,5))),
             linetype="dashed", size=1)


mu
mk

sk2

