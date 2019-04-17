rm(list=ls())
source("functions.R")

n.cue = 1
n.alt = 3
n.node = n.cue + n.alt

omega = matrix(0,n.node,n.node)
omega[1,] = omega[,1] = c(0,1,1,1)

# Debreu Example ----------------------------------------------------------

omega.debreu = omega
colnames(omega.debreu) = rownames(omega.debreu) = c("R","B_F","B_K","D_C")
omega.debreu[2,3] = omega.debreu[3,2] = -10

mu.debreu = c(0,3,3,3.40548)
beta.debreu = 1
vareps.debreu = 0

res.debreu.FK = ecb(omega.debreu[-c(4),-c(4)],mu.debreu[-c(4)],beta.debreu,vareps.debreu, plot=T)
res.debreu.FC = ecb(omega.debreu[-c(3),-c(3)],mu.debreu[-c(3)],beta.debreu,vareps.debreu, plot=T)
res.debreu.KC = ecb(omega.debreu[-c(2),-c(2)],mu.debreu[-c(2)],beta.debreu,vareps.debreu, plot=T)
res.debreu.FKC = ecb(omega.debreu,mu.debreu,beta.debreu,vareps.debreu, plot=T)

# Choice probabilties for choosing B_K or D_C given that we already have B_F
res.debreu.KC.F = ecb(omega.debreu,mu.debreu,beta.debreu,vareps.debreu, plot=T, n.cue=2)


# Pens Example A ----------------------------------------------------------

omega.pens.a = omega
colnames(omega.pens.a) = rownames(omega.pens.a) = c("R","$","P_+","P_-")
omega.pens.a[2,3] = omega.pens.a[3,2] = -.59
omega.pens.a[2,4] = omega.pens.a[4,2] = -2.75

mu.pens.a = c(0,4.25,3.675,.825)
beta.pens.a = 1
vareps.pens.a = 0

res.pens.a.mpp = ecb(omega.pens.a[-c(4),-c(4)],mu.pens.a[-c(4)],beta.pens.a,vareps.pens.a, plot=T)
res.pens.a.mpm = ecb(omega.pens.a[-c(3),-c(3)],mu.pens.a[-c(3)],beta.pens.a,vareps.pens.a, plot=T)
res.pens.a.ppm = ecb(omega.pens.a[-c(2),-c(2)],mu.pens.a[-c(2)],beta.pens.a,vareps.pens.a, plot=T)
res.pens.a.mppm = ecb(omega.pens.a,mu.pens.a,beta.pens.a,vareps.pens.a, plot=T)

# Choice probabilties for choosing $ or P_- given that we already have P_+
omega.pens.a.pp = omega.pens.a[c(1,3,2,4),c(1,3,2,4)]
mu.pens.a.pp = mu.pens.a[c(1,3,2,4)]
res.pens.a.mpm.pp = ecb(omega.pens.a.pp,mu.pens.a.pp,beta.pens.a,vareps.pens.a, plot=T, n.cue=2)

# Pens Example B ----------------------------------------------------------

omega.pens.b = omega
colnames(omega.pens.b) = rownames(omega.pens.b) = c("R","$","P_+","P_-")
omega.pens.b[3,4] = omega.pens.b[4,3] = 2.75

mu.pens.b = c(0, 1.750, 1.175, -1.675)
beta.pens.b = 1
vareps.pens.b = 0

res.pens.b.mpp = ecb(omega.pens.b[-c(4),-c(4)],mu.pens.b[-c(4)],beta.pens.b,vareps.pens.b, plot=T)
res.pens.b.mpm = ecb(omega.pens.b[-c(3),-c(3)],mu.pens.b[-c(3)],beta.pens.b,vareps.pens.b, plot=T)
res.pens.b.ppm = ecb(omega.pens.b[-c(2),-c(2)],mu.pens.b[-c(2)],beta.pens.b,vareps.pens.b, plot=T)
res.pens.b.mppm = ecb(omega.pens.b,mu.pens.b,beta.pens.b,vareps.pens.b, plot=T)

# Choice probabilties for choosing $ or P_- given that we already have P_+
omega.pens.b.pp = omega.pens.b[c(1,3,2,4),c(1,3,2,4)]
mu.pens.b.pp = mu.pens.b[c(1,3,2,4)]
res.pens.b.mpm.pp = ecb(omega.pens.b.pp,mu.pens.b.pp,beta.pens.b,vareps.pens.b, plot=T, n.cue=2)


# Camera Example ----------------------------------------------------------

omega.camera = omega
colnames(omega.camera) = rownames(omega.camera) = c("C","C_L","C_M","C_H")
omega.camera[2,4] = omega.camera[4,2] = -15
omega.camera[2,3] = omega.camera[3,2] = omega.camera[4,3] = omega.camera[3,4] = 0


mu.camera = c(0,7.5,7.5,7.5)
beta.camera = 1
vareps.camera = 0

res.camera.LM = ecb(omega.camera[-c(4),-c(4)],mu.camera[-c(4)],beta.camera,vareps.camera, plot=T)
res.camera.LH = ecb(omega.camera[-c(3),-c(3)],mu.camera[-c(3)],beta.camera,vareps.camera, plot=T)
res.camera.MH = ecb(omega.camera[-c(2),-c(2)],mu.camera[-c(2)],beta.camera,vareps.camera, plot=T)
res.camera.LMH = ecb(omega.camera,mu.camera,beta.camera,vareps.camera, plot=T)





