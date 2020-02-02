
#---------------------------------------AFD-------------------------------------------------
library(ade4)
library(FactoMineR)

My_data<-read.csv2("chienloup.csv")

#les données 
data_table<-as.matrix(My_data[2:7])

#Centralisation des variables 
data_table<-data_table-t(replicate(n,g))


n<-nrow(data_table)
p<-ncol(data_table)
n1<-30
n2<-12
nbr_classes<-2

Chien<-data_table[1:n1,]
Loup<-data_table[(n1+1): 42,]




#-----------Question 1: calcul et comparaison desvariances-----------

#les centres d'inertie
g<-colMeans(data_table)
g1<-colMeans(Chien)
g2<-colMeans(Loup)

# La matrice de variance-covariance V 
V<-matrix(0,nrow=p,ncol=p)
for ( i in 1:n)
{
  V<-V+(data_table[i,]-g)%*%t(data_table[i,]-g)
  
}
V<-V/n

#La matrice inter-classes B (Between)
nk<-c(n1,n2)
B<-(g1-g)%*%t(g1-g)*n1+(g2-g)%*%t(g2-g)*n2
B<-B/n
print(B)

# Les matrices intra-classes 
W1<-cov(Chien)*(n1-1)/n1  # on utilise cet astuce parce que  
W2<-cov(Loup)*(n2-1)/n2   #la covariance sous R est non biasée.

W<-n1*W1+n2*W2
W<-W/n
print(W)

#compariason des variances
print(V>B)
print(V>W)

#Les variables initiales partagent le groupe.

print(V/(B+W))

#---------Question 2: étude de V-1 et B -----------------------------------------

# La matrice B est une matrice symétrique  positive ==> Diagonalisable 
#et de dimension p*p
# La matrice V est une matrice symétrique  positive ==> Diagonalisable
#et de dimension p*p

V_inverse<-solve(V) # la matrice V_(-1) est une matrice symétrique
                    #(car V est symetrique ) et de dimension p*p


d<-V_inverse%*%B    # La matrice n'est pas forcément diagonalisable :
                    # On peut pas la diagonaliser directement

#----------------Question 3: Représentation de chaque groupe d'individu--------

#Diagonalisation de la matrice (V-1B)
#Passage par l'intermidiare C

C<-matrix(0,nrow=p,ncol=nbr_classes)
C[,1] <-(g1-g)*sqrt(n1/n)
C[,2]<-(g2-g)*sqrt(n2/n)

eigen_matrix<-eigen(t(C)%*%V_inverse%*%C)


lambda<-eigen_matrix$values
barplot(lambda)
vectors<-(eigen_matrix$vectors)

#Centralisation des variables 
data_table<-data_table-t(replicate(n,g))

#Les axes factoriels
F1<-data_table%*%vectors[,1]
F2<-data_table%*%vectors[,2]
F1C<-F1[1:n1]
F2C<-F2[1:n1]
 plot(F1C,F2C,pch=15,col="red",xlim=c(90,200),ylim = c(-70,0))
F1L<-F1[n1+1:n]
F2L<-F2[n1+1:n]
points(F1L,F2L,pch=18,col="blue")
abline(v=133,lty=2)
title("Discirmination des classes Chien et Loup")
 
#Les fonctions des packages ade4 et factoMineR
FAMD(My_data,ncp=2,graph=T)
dudi.mix(My_data)

#-----------------Question 6: rapport de corrélation----------------------------
mu<-integer(p)
for (j in 1:p)
{
Z<-data_table%*%vectors[,j]
SCT<-sum((Z-mean(Z))*(Z-mean(Z)))
SCE<-n1*(mean(Z[1:n1])-mean(Z)) + n2*(mean(Z[(n1+1):n])-mean(Z))
 mu[j]<-SCE/SCT
}

barplot(mu)



#-----------Question 7: comparaison des valeurs et vecteurs propres. -------------------



W_inverse<-solve(W)

eigen(W_inverse%*%B)

eigen(V_inverse%*%B)
------------



#------------Question 9: Les équations des axes factoriels  ---------------------

# multiplication des tables de données par les directions principles

#------------Question 10: Indicateurs de la qualité de projection---------------

#Qualité de projection d'un nuage
  #axe 1:
QN1<-lambda[1]/(sum(lambda))
  #axe 2:
QN2<-lambda[2]/(sum(lambda))

#contribution absolue du centre g_q à l'axe a-s:
   #axe 1:
CA1<-n1*(a1%*%V_inverse%*%g1)*(a1%*%V_inverse%*%g1)/n
  #axe 2:
CA2<-n2*(a2%*%V_inverse%*%g2)*(a2%*%V_inverse%*%g2)/n

#Contribution relative du centre de gravité g_q à l'axe a_s
  #axe1:
CR1<-n1*(a1%*%V_inverse%*%g1)*(a1%*%V_inverse%*%g1)/(n*lambda[1])
  #axe 2:
CA2<-n2*(a2%*%V_inverse%*%g2)*(a2%*%V_inverse%*%g2)/(n*lambda[2])

----




