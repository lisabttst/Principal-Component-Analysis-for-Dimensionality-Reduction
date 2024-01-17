# TP1
library(psych)
library(Matrix)
library(matlib)

X <- t(matrix(c(0, 1, 0, 1, 1, 0, 0, 1, 1, 0, 1, 0), 2, byrow = TRUE))
D <- 1/6*diag(6)
colnames(X) = c("X", "Y")
X
D

# Poids des individus 1/6, donc D = 1/6 * I6

# Moyenne des variables : 
m_X_Y <- t(rep(1/6, 6) %*% X)
m_X_Y

# Centrer les données : 
Z <- (X - 0.5)
Z
# Matrice de variance-covariance : 
# S2(X) = S2(Y) = 0.25 et cov(X,Y) = 0.083
A <- (5/6)*cov(Z)
A

# Valeurs propres :
lambda <- eigen(A)$values
lambda
cat("Tr(A) :", tr(A), "_ et somme des vp :", lambda[1] + lambda[2])
cat("det(A) :", det(A), "_ et produit des vp :", lambda[1] * lambda[2])

# Vecteurs propres
u <- eigen(A)$vectors
u1 <- c(A[2,2]-lambda[1], -A[1,2])
u2 <- c(A[2,2]-lambda[2], -A[1,2])

u1 <- -(1/norm(u1, type="2")) * u1
u2 <- -(1/norm(u2, type="2")) * u2

u1
u2
u

# Vérifier égalité inertie statistique
is_u1 <- 1/6*sum((Z %*% u[,1])**2)
is_u2 <- 1/6*sum((Z %*% u[,2])**2)
is_u1
is_u2

t(u1) %*% A %*% u1
lambda[1] * t(u1) %*% u1
lambda[1] * norm(u1, type="2")**2
lambda[1]

# Inertie totale
it <- (1/6*sum(Z**2))
cat("Inertie totale :", it, "_ et somme des inerties statistiques :", is_u1 + is_u2)
cat("Somme des VP : ", lambda[2] + lambda[1])

# Trace les deux axes principaux C1 et C2
plot(X[,1], X[,2], xlim=c(-1.5,1.5), ylim=c(-1.5,1.5), xlab="X", ylab="Y")
arrows(0,0,u1[1],u1[2], length=0.1)
arrows(0,0,u2[1],u2[2], length=0.1)

# 3.7 Taux d'inertie expliquée
T1 <- lambda[1] / (lambda[1]+lambda[2])
T2 <- lambda[2] / (lambda[1]+lambda[2])
cat("Taux d'inertie expliquée par la 1ère composante :", T1*100, "% _ et par la 2ième :", T2*100, "%")

# Calcule les projections des 6 individus sur u1 et u2
coord_proj <- Z %*% u
coord_proj
plot(coord_proj, col="red", panel.first=grid())

# Qualité de projection des individus
qual <- 



