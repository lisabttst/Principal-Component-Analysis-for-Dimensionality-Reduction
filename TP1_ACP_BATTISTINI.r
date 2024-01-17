# TP1
# Jules PAOLO
# Coralie GIBOUDEAU
# Lisa BATTISTINI
# Antoine LEROY

library(psych)
library(Matrix)
library(matlib)

## Création des matrices X et D
X <- t(matrix(c(0, 1, 0, 1, 1, 0, 0, 1, 1, 0, 1, 0), 2, byrow = TRUE))
colnames(X) = c("X", "Y")
D <- 1/6*diag(6)


## Fonction qui réalise l'ACP d'une matrice quelconque de taille n.p
ACP <- function(X, norm = FALSE) {
  
  n = dim(X)[1]
  p = dim(X)[2]
  
  # Convertit les valeurs en nombres
  X <- as.numeric(unlist(X))
  X <- matrix(X, nrow=n, ncol=p)
  
  # Centre X par défaut, réduit si norm == TRUE
  for (i in 1:p) {
    X[,i] <- X[,i] - mean(X[,i])
    if(norm == TRUE) {
      X[,i] <- X[,i] / (((n-1)/n) * var(X[,i]))
    }
  }
  
  # Matrice de covariance A
  A <- ((n-1)/n)*cov(X)
  
  # Diagonalise A
  lambda <- eigen(A)$values
  u <- eigen(A)$vectors
  
  # Détermine les composantes principales
  comp <- X %*% u
  
  return(list("Norm"=X, "Cov"= A, "Val_p"=lambda, "Vect_p"=u, "Comp"=comp))
}


## Stock les variables de l'ACP dans une liste L
L <- ACP(X, norm = FALSE)


## Vérification des valeurs propres
lambda <- L$Val_p
cat("Tr(A) :", tr(L$Cov), "_ et somme des vp :", lambda[1] + lambda[2])
cat("det(A) :", det(L$Cov), "_ et produit des vp :", lambda[1] * lambda[2])


## Vérifications des vecteurs propres par calcul et par diagonalisation
u <- L$Vect_p
u1 <- c(L$Cov[2,2]-lambda[1], -L$Cov[1,2])
u2 <- c(L$Cov[2,2]-lambda[2], -L$Cov[1,2])

u1 <- -(1/norm(u1, type="2")) * u1
u2 <- -(1/norm(u2, type="2")) * u2

print("Calcul vecteurs propres")
print(u1)
print(u2)
print("Vecteurs propres données par la fonction ACP")
print(u)


## Vérifications égalité inertie statistique
is_u1 <- 1/6*sum((L$Norm %*% u[,1])**2)
is_u2 <- 1/6*sum((L$Norm %*% u[,2])**2)
is_u1
is_u2

t(u1) %*% L$Cov %*% u1
lambda[1] * t(u1) %*% u1
lambda[1] * norm(u1, type="2")**2
lambda[1]


## Vérification de l'inertie totale
it <- (1/6*sum(L$Norm**2))
cat("Inertie totale :", it, "_ et somme des inerties statistiques :", is_u1 + is_u2)
cat("Somme des VP : ", lambda[2] + lambda[1])


## Trace les deux axes principaux C1 et C2
Z <- L$Norm
plot(Z[,1], Z[,2], xlim=c(-1.5,1.5), ylim=c(-1.5,1.5), xlab="X", ylab="Y")
arrows(0,0,u1[1],u1[2], length=0.1, col="red")
arrows(0,0,u2[1],u2[2], length=0.1, col="blue")
legend(-1.3, -0.6, legend=c("Axe 1", "Axe 2"), col=c("red", "blue"), lty=1:1, cex=0.8)


## Taux d'inertie expliquée par chacun des deux axes
T1 <- lambda[1] / (lambda[1]+lambda[2])
T2 <- lambda[2] / (lambda[1]+lambda[2])
cat("Taux d'inertie expliquée par la 1ère composante :", T1*100, "% _ et par la 2ième :", T2*100, "%")


## Plot les projections des 6 individus sur les composantes C1 et C2
plot(L$Comp, col="red", panel.first=grid(), xlab="C1", ylab="C2")
L$Comp


## Calcule la qualité de projection des individus
qual <- rep(0,dim(L$Comp)[2])
for (i in 1:dim(L$Comp)[1]) {
    qual[i] <-  sum(L$Comp[i,][2]**2) / sum(L$Comp[i,]**2)
}
print(qual)


## Contribution à une composante
contribution <- L$Comp
colnames(contribution) = c("Contribution sur c1", "Contribution sur c2")
for (i in 1:dim(L$Comp)[1]) {
  for (j in 1:dim(L$Comp)[2])
  contribution[i,j] <- (1/6)*L$Comp[i,j]**2 / lambda[j]
}
contribution
