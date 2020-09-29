#1. 

r1 <- c(2, -3, -7, 5, 2, -2)
r2 <- c(1, -2, -4, 3, 1, -2)
r3 <- c(2, 0, -4, 2, 1, 3)
r4 <- c(1, -5, -7, 6, 2, -7)
rbind(r1, r2, r3, r4)
A <- matrix(c(2,1,2,1,-3,-2,0,-5,-7,-4,-4,-7,5,3,2,6,2,1,1,2),4); w <- c(-2,-2,3,-7)
A; w

r1 <- r1 - r2; rbind(r1, r2)
r2 <- r2 - r1; rbind(r1, r2)
r3 <- r3 - 2*r1; rbind(r1, r2, r3)
r4 <- r4 - r1; rbind(r1, r2, r3, r4)
r2 <- -r2; rbind(r1, r2)
r3 <- r3 - 2*r2; rbind(r1, r2, r3)
r4 <- r4 + 4*r2; rbind(r1, r2, r3, r4)
r1<- r1 + r2; rbind(r1, r2)
r3<- -r3; rbind(r1, r2, r3)
r4<- r4-r3; rbind(r1, r2, r3, r4)
r1<-r1-r3
rbind(r1,r2,r3,r4)

#we found many solutions
a5<-1
a1<- 1-(-2)*a3-a4; a2<-2-a3+a4;

a3 <- 2; a4 <-3;  v<- c(a1,a2,a3,a4,a5); A%*% v
a3 <- 4; a4 <-3;  v<- c(a1,a2,a3,a4,a5); A%*% v
a3 <- 4; a4 <-1;  v<- c(a1,a2,a3,a4,a5); A%*% v


#2. 

#Let's work in Z_5
"%+5%" <- function(x,y) (x+y) %%5  #addition
"%-5%" <- function(x,y) (x-y) %%5  #subtraction
"%*5%" <- function(x,y) (x*y) %%5  #multiplication
"%/5%" <- function(x,y) (x*y*y*y) %%5  #division

f1<-c(3,0,4,0,2,2)
f2<-c(1,1,3,3,2,1)
f3<-c(0,2,1,1,4,2)
f4<-c(1,0,2,0,3,4)
F<-rbind(f1,f2,f3,f4);F

temp_f <- f4
f4 <- f1
f1 <-temp_f
f2 <- f2%-5%f1
rbind(f1,f2,f3,f4)

f4<-f4%-5%(3 %*5% f1)
f3<-f3%-5%(2 %*5% f2)
rbind(f1,f2,f3,f4)

f4<-f4%/5%3
f3<-f3%-5%(4 %*5% f4)
rbind(f1,f2,f3,f4)

temp_f1 <-f3
f3<-f4
f4<-temp_f1
f4<-f4%/5%2
rbind(f1,f2,f3,f4)

f1<-f1%-5%(2%*5%f3)
f3<-f3%-5%f4
rbind(f1,f2,f3,f4)

f2<-f2%-5%f3
rbind(f1,f2,f3,f4)

f1<-f1%-5%f4
f2<-f2%-5%(4%*5%f4)
rbind(f1,f2,f3,f4)

#basis of image are the pivot columns of the original matrix
F
F[,1]; F[,2]; F[,3];F[,5]

#basis of kernel 
#To find a basis for the kernel, use the two nonpivotal columns.
#Vectors of the form (x1 y1 z1 1 w1 0) and (x2 y2 z2 0 w2 1) must be independent,
#because no linear combination can have 0 as its 2nd and 4th component.
#The top row of the row reduced matrix applied to (x1 y1 z1 1 w1 0) says that x1=0
#The second row of the row reduced matrix applied to (x1 y1 z1 1 w1 0) says that y1+3=0; so y1 = -3
#the third row says z1+1*0 = 0; so z1 = 0 
#fourth = 1
#the fifth row says w1+4(0)=0
rbind(f1,f2,f3,f4)
k1 <- c(0,-3,0,1,0,0)  
#first basis vector for the kernel
F%*%k1%*5%1     #yes, it is in the kernel

#The top row of the row reduced matrix applied to (x1 y1 z1 0 w1 1) says that x1=0
#The second row of the row reduced matrix applied to (x1 y1 z1 0 w1 1) says that y1=0
#the third row says z1+1= 0; so z1 = -1
#fourth = 0
#the fifth row says w1+4(1)=0; w1=-4
rbind(f1,f2,f3,f4)
k2 <- c(0,0,-1,0,-4,1)  
#second basis vector for the kernel
F%*%k2%*5%1     #yes, it is in the kernel


#3. 
v1<-c(1,1,1)
v2<-c(1,0,1)
v3<-c(1,0,0)
rref(cbind(v1,v2,v3,c(0,0,0)))

#Step 1: make the first unit vector by normalization
v_1 <- v1/Norm(v1); v_1; Norm(v_1)
#Step 2: convert w2 to a vector that is orthogonal to v1.
x <- v2 - (v2%.%v_1)*v_1; x%.% v_1
#Then convert x to a unit vector
v_2 <- x/Norm(x); round(v_2)
#Step 3: convert w3 to a vector that is orthogonal to both v1 and v2.
x_ <- v3 - (v3%.%v_1)*v_1 - (v3%.%v_2)*v_2; x_%.% v_1; x_%.% v_2
#Then convert x to a unit vector
v_3 <- x_/Norm(x_); round(v_3)

C <- cbind(v_1,v_2, v_3) ;C  #basis vectors are the columns
round(t(C)%*%C, digits = 6)   #the identity matrix

det(C)
