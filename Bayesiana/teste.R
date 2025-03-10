set.seed(123)
n <- 30
x <- rpois(n, 5)
xhat <- mean(x)

score <- function(lambda) {-n+sum(x)/lambda}
curve(score, from=2, to =8, xlab='Lambda', ylab='Score')
abline(a=1, b=-0.2, col='red')
legend('topright', lty=1, col=c('black', 'red'), c('Score function', 'Tangent at root'))

n/xhat ## expected information

n <- 1000
x <- rpois(n, 5)
xhat <- mean(x)

score <- function(lambda) {-n+sum(x)/lambda}
curve(score, from=2, to =8, xlab='Lambda', ylab='Score')
abline(a=1, b=-0.2, col='red')
legend('topright', lty=1, col=c('black', 'red'), c('Score function', 'Tangent at root'))

n/xhat ## expected information
