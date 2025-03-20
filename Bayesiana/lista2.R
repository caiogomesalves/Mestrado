# limpar todas as variáveis
rm(list = ls(all.names = TRUE))


library(pscl)

# Comparação de estimadores (Verossimilhança Bernoulli)

n <- 30
theta <- seq(0.0001,0.9999,length.out=1000)
EQM.MV <- theta*(1-theta)/n

## diferentes valores para a e b, parametros da priori Beta(a,b)
m.ab <- rbind(cbind(1,1),cbind(5,10),cbind(10,5),cbind(10,10))

pdf("compare-n=30.pdf")
par(mfrow=c(2,2))

for(i in 1:4)
{
  EQM.EAP <- (n*theta*(1-theta) + (m.ab[i,1]-theta*(m.ab[i,1] + m.ab[i,2]))^2)/((n+m.ab[i,1]+m.ab[i,2])^2)

  EQM.MAP <- ((n * theta * (1 - theta)) + (m.ab[i, 1] - 1 - theta * (m.ab[i, 1] + m.ab[i, 2] - 2))^(2))/((n + m.ab[i, 1] + m.ab[i, 2] - 2)^(2))

  #EQM EMV
  plot(theta,EQM.MV,type="l",lwd=3,xlab=expression(theta),ylab="EQM",cex=1.3,
       cex.lab=1.3,cex.axis=1.3,
       main=paste("a = ",m.ab[i,1],", b = ",m.ab[i,2],", n = ",n),col=1)

  lines(theta,EQM.EAP,type="l",lwd=3,xlab=expression(theta),ylab="EQM",cex=1.3,
        cex.lab=1.3,cex.axis=1.3,lty=2,col=2)

  lines(theta,EQM.MAP,type="l",lwd=3,xlab=expression(theta),ylab="EQM",cex=1.3,
        cex.lab=1.3,cex.axis=1.3,lty=2,col=3)

  legend("center", legend = rbind("MV","EAP", "MAP"),bty="n",lty=c(1,2, 2),
         cex=1.2,lwd=2,col=c(1,2, 3))
}
dev.off()

n <- 100
EQM.MV <- theta*(1-theta)/n

pdf("compara-n100.pdf")
par(mfrow=c(2,2))

for(i in 1:4)
{
  EQM.EAP <- (n*theta*(1-theta) + (m.ab[i,1]-theta*(m.ab[i,1] + m.ab[i,2]))^2)/((n+m.ab[i,1]+m.ab[i,2])^2)

  #EQM EMV
  plot(theta,EQM.MV,type="l",lwd=3,xlab=expression(theta),ylab="EQM",cex=1.3,
       cex.lab=1.3,cex.axis=1.3,
       main=paste("a = ",m.ab[i,1],", b = ",m.ab[i,2],", n = ",n),col=1)

  lines(theta,EQM.EAP,type="l",lwd=3,xlab=expression(theta),ylab="EQM",cex=1.3,
        cex.lab=1.3,cex.axis=1.3,lty=2,col=2)

  legend("center", legend = rbind("MV","EAP"),bty="n",lty=c(1,2),
         cex=1.2,lwd=2,col=c(1,2))
}
dev.off()

n <- 2000
EQM.MV <- theta*(1-theta)/n

pdf("compara-n2000.pdf")
par(mfrow=c(2,2))
for(i in 1:4)
{
  EQM.EAP <- (n*theta*(1-theta) + (m.ab[i,1]-theta*(m.ab[i,1] + m.ab[i,2]))^2)/((n+m.ab[i,1]+m.ab[i,2])^2)

  #EQM EMV
  plot(theta,EQM.MV,type="l",lwd=3,xlab=expression(theta),ylab="EQM",cex=1.3,
       cex.lab=1.3,cex.axis=1.3,
       main=paste("a = ",m.ab[i,1],", b = ",m.ab[i,2],", n = ",n),col=1)

  lines(theta,EQM.EAP,type="l",lwd=3,xlab=expression(theta),ylab="EQM",cex=1.3,
        cex.lab=1.3,cex.axis=1.3,lty=2,col=2)

  legend("center", legend = rbind("MV","EAP"),bty="n",lty=c(1,2),
         cex=1.2,lwd=2,col=c(1,2))
}
dev.off()


# Superfícies de EQM: n x theta

# EMV
theta<- seq(0.0001,0.9999,length.out=200)
n <- seq(30,50,length.out=200)

eqm.mv<-function(theta,n){
  theta*(1-theta)/n
}

## grafico de superficie do EMV e EAP

z <- outer(theta,n,eqm.mv)
pdf("Surface_EQM_MV.pdf")
par(mfrow=c(1,1))

persp(theta,n,z,theta=65,phi=20,expand=0.5,ticktype = "detailed",shade=0.75,
      col=5,zlab="EQM",main="Estimador de máxima verossimilhança")

dev.off()

# EAP

theta<- seq(0.0001,0.9999,length.out=200)
n <- seq(30,50,length.out=200)

a <- 10
b <- 10

eqm.EAP <- function(theta,n){
  (n*theta*(1-theta) + (a-theta*(a + b))^2)/((n+a+b)^2)
}

z <- outer(theta,n,eqm.EAP)

pdf("Surface_EQM_EAP.pdf")
par(mfrow=c(1,1))

persp(theta,n,z,theta=65,phi=20,expand=0.5,ticktype = "detailed",shade=0.75,
      col=5,zlab="EQM",main=paste("Esperança a posteriori, a = ",a,", b = ",b))
dev.off()

pdf("priori-beta-10-10.pdf")
s.beta <- rbeta(1000,10,10)
plot(density(s.beta),main="priori Beta(10,10)",xlab="theta",ylab="densidade")
dev.off()

#----Meu método----

func_grafico_eqm <- function(n, a, b) {
    m <- cbind(a, b)
    if (nrow(m) == 1) {
        theta <- seq(0, 1, length.out = 1000)
        eqm_mle <- (theta * (1 - theta))/n
        eqm_eap <- (n * theta * (1 - theta) + (a - theta * (a + b))^2)/((n + a + b)^2)
        eqm_map <- ((n * theta * (1 - theta)) + (a - 1 - theta * (a + b - 2))^(2))/((n + a + b - 2)^(2))
        plot(theta, eqm_mle, type = "l", xlab = expression(theta), ylab = "EQM")
        lines(theta, eqm_eap, type = "l", lty = 2, col = "red", lwd = 2)
        lines(theta, eqm_map, type = "l", lty = 3, col = "blue", lwd = 2)
    }
    else {
        par(mfrow = c(round(nrow(m)/2), 2))
        theta <- seq(0, 1, length.out = 1000)
        eqm_mle <- (theta * (1 - theta))/n
        for (i in 1:nrow(m)) {
            eqm_eap <- (n * theta * (1 - theta) + (m[i, 1] - theta * (m[i, 1] + m[i, 2]))^2)/((n + m[i, 1] + m[i, 2])^2)
            eqm_map <- ((n * theta * (1 - theta)) + (m[i, 1] - 1 - theta * (m[i, 1] + m[i, 2] - 2))^(2))/((n + m[i, 1] + m[i, 2] - 2)^(2))
            plot(theta, eqm_mle, type = "l", xlab = expression(theta), ylab = "EQM", main = bquote(a == .(m[i, 1]) ~~~ b == .(m[i, 2])))
            lines(theta, eqm_eap, type = "l", lty = 2, col = "red", lwd = 2)
            lines(theta, eqm_map, type = "l", lty = 3, col = "blue", lwd = 2)
        }
    }
}

func_grafico_eqm(100, rep(c(10, 2), 3), rep(4:6, each = 2))

par(mfrow = c(1, 1))
