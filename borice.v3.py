# -*- coding: utf-8 -*-

# v3: seperate t (outcrossing rate) for each population
# v3: Output posterior density for Dxy for each pairwise contrast of populations
# v3: Allow runs on no data to investigate priors
# v3: MCMC replication ID entered as command line argument

# full sib half sib structure.  delta[j] is now the "paternal index": 0 = selfed, >=1 is outcrossed, number = sire ID
# Essentially same as dt2 but with different structure for likelihood

# UpdateDelta() allows full sibs
# Include updates to AFL (assume exponential prior)
# Run a bunch of cycles for allele freqs first, do deltaupdate only at thinning freq

import copy
from math import pow, log, exp, gamma
from random import random
from random import randint
from numpy.random import beta
import scipy.stats 
import sys

rep="run1" # sys.argv[1]

MinCdxy = 4 # min called individuals at a SNP in population to count for Dxy

Bigbad=-999999

inx=open("Control.v3.txt", "r")
for line_idx, line in enumerate(inx):
    cols = line.split('\t') 
    u1=cols[0].split("=")
    if u1[0]=="InputDataType":
        InputDataType=int(u1[1])
        if InputDataType != 1:
            print("Using wrong program version for this input data type ")
            break
    elif u1[0]=="FILEPREFIX":
        FILEPREFIX=u1[1]        # main input file is [FILEPREFIX].genotypes.txt 
    elif u1[0]=="RUNTYPE":
        RUNTYPE=u1[1]            # 
    elif u1[0]=="NoSubpops":
        NoSubpops=int(u1[1])        # number of subpopulations (details in popfile)
    elif u1[0]=="NoFamilies":
        NoFamilies=int(u1[1])        # number of families
    elif u1[0]=="popfile":
        popfile=u1[1]            # file containing location of each maternal family
    elif u1[0]=="FullOutput":
        FullOutput=int(u1[1])        # set to 1 to output posterior probs for all genotypes
    elif u1[0]=="IHPriorModel":
        IHPriorModel=int(u1[1])        # 0=IH priors determined by t; 1= IH prior with uniform F
    elif u1[0]=="ChainLength":
        ChainLength=int(u1[1])        # run duration
    elif u1[0]=="afstep":
        afstep=float(u1[1])        # allele frequency step window size
    elif u1[0]=="tstep":
        tstep=float(u1[1])        # t step window size
    elif u1[0]=="BetaGstep":
        BetaGstep=float(u1[1])        # step size for BetaGu (relevant only for InputDataType=2)
    elif u1[0]=="AFboundary":
        AFboundary=float(u1[1])        # minimum for Minor Allele Frequency
    elif u1[0]=="thinningfreq":
        thinningfreq=int(u1[1])        # records taken from chain every thinningfreq steps
    elif u1[0]=="burnin":
        burnin=int(u1[1])        # length of burn-in
    elif u1[0]=="meanFst":
        meanFst=float(u1[1])        # Prior for degree of subpopulation differentiation
    elif u1[0]=="vfacFst":
        vfacFst=float(u1[1])        # variance factor for Fst Prior
    elif u1[0]=="alfstep":
        alfstep=float(u1[1])
    elif u1[0]=="Tuning_Steps":
        TuningSteps=int(u1[1])

    
inx.close()

FpopPrior = [meanFst,vfacFst*meanFst*(1-meanFst)] # mean, variance of beta prior
alphaFpop =FpopPrior[0]*(FpopPrior[0]*(1-FpopPrior[0])/FpopPrior[1] - 1)
betaFpop =(1-FpopPrior[0])*(FpopPrior[0]*(1-FpopPrior[0])/FpopPrior[1] - 1)
print ("F pop prior ",alphaFpop, betaFpop)


# Mendelism [mom genotype][dad genotype][off]

OProb=[ [[1.0,0,0], [0.5,0.5,0], [0,1.0,0]], [[0.5,0.5,0], [0.25,0.5,0.25], [0,0.5,0.5]], [[0,1.0,0], [0.0,0.5,0.5], [0,0,1.0]] ] # UPDATE: NOW A 3x3 ARRAY (DIPLOID MALE GENOTYPE) for outcrossed progeny

SProb=[ [1.0,0,0], [0.25,0.5,0.25], [0,0,1.0] ] # selfed progeny
Fmom = [0,0.5,0.75,0.875,0.9375,0.96875,0.984375,1.0] # F given IH


# Genotyping probilities, "RA" values are updated for each snp 
GProb={"RR":{},"RA":{},"AA":{}}
GProb["RR"]["R"]=1.0 # [true genotype][observed data]
GProb["RR"]["H"]=0.0
GProb["RR"]["A"]=0.0
GProb["RR"]["N"]=1.0
GProb["RA"]["R"]=0.5
GProb["RA"]["H"]=0.5
GProb["RA"]["A"]=0.5
GProb["RA"]["N"]=1.0
GProb["AA"]["R"]=0.0
GProb["AA"]["H"]=0.0
GProb["AA"]["A"]=1.0
GProb["AA"]["N"]=1.0

########## Likelihood calculators


def fam1(Data,famID,snpID,delt,use1): # family j probability for one SNP

    pop=SubPop[famID]
    q = Q[pop][snpID]
    f = Fmom[IH[famID]] 
    PriorMom={}
    PriorMom[0]= q*q*(1-f) + f*q
    PriorMom[1]= 2*q*(1-q)*(1-f)
    PriorMom[2]= 1.0-PriorMom[0]-PriorMom[1]
    PriorDad={}
    PriorDad[0]= q*q*(1-Fadultmean) + Fadultmean*q
    PriorDad[1]= 2*q*(1-q)*(1-Fadultmean)
    PriorDad[2]= 1.0-PriorDad[0]-PriorDad[1]

    dg=0
    Famsize=len(Data)
    

    for k in range(Famsize): # check for dataless fams
        if Data[k][0] < 1.0 or Data[k][1] < 1.0 or Data[k][2] < 1.0:
            dg = 1
    if dg==0:
        if use1==0:
            return 1.0 # no data for family
        else:
            return -9,-9,-9


    else: # calculate mom prob

        if dg==1 and Famsize==1: # only a mom that has data (she produced no genotyped progeny)
            if use1==0:
                return PriorMom[0]*Data[0][0] + PriorMom[1]*Data[0][1] + PriorMom[2]*Data[0][1] 
            else:
                return PriorMom[0]*Data[0][0],PriorMom[1]*Data[0][1],PriorMom[2]*Data[0][1] 

        elif dg==1 and Famsize>1: 

            Sires = max(delt) # UPDATE

            FamLikelihood=0.0
            mpp={}
            for MG in range(3): # condition on maternal genotype

                if MG==0:
                     Fk=Data[0][0]
                elif MG==1:
                     Fk=Data[0][1]
                elif MG==2:
                     Fk=Data[0][2]


                ProductX=[1.0 for xval in range(Sires+1)]
                for xval in range(Sires+1):

                    Ls=[1.0,1.0,1.0]
                    for offspring in range(1,Famsize): # now over sires

                        if xval==0 and delt[offspring-1]==0:  # selfed
                            ProductX[xval] *= ( SProb[MG][0]*Data[offspring][0] + SProb[MG][1]*Data[offspring][1] + SProb[MG][2]*Data[offspring][2] )
                        
                        elif delt[offspring-1]==xval: # outbred with Sire xval
                            Ls[0]*= OProb[MG][0][0]*Data[offspring][0] + OProb[MG][0][1]*Data[offspring][1] + OProb[MG][0][2]*Data[offspring][2] # dad is RR
                            Ls[1]*= OProb[MG][1][0]*Data[offspring][0] + OProb[MG][1][1]*Data[offspring][1] + OProb[MG][1][2]*Data[offspring][2] # dad is RA
                            Ls[2]*= OProb[MG][2][0]*Data[offspring][0] + OProb[MG][2][1]*Data[offspring][1] + OProb[MG][2][2]*Data[offspring][2] # dad is AA
                    if xval>0:
                        ProductX[xval]*= (PriorDad[0]*Ls[0] + PriorDad[1]*Ls[1] + PriorDad[2]*Ls[2])
                    Fk*= ProductX[xval]
                FamLikelihood+=PriorMom[MG]*Fk
                mpp[MG]=PriorMom[MG]*Fk

            if use1==0:
                return FamLikelihood
            else:
                return mpp[0],mpp[1],mpp[2]



def allfams1(snpID): # whole dataset probability for one SNP:: 

    ln_l = 0.0
    for j in range(1,NoFamilies+1):
        famL=fam1(FullData[snpID][j],j,snpID,delta[j],0)
        if famL>0.0:
            ln_l+=log(famL)
        else:
            print("allfams1, Impossible fam ",j,snpID,FullData[snpID][j])
            print("delta ",delta[j])
            return float('-inf')
    return ln_l


def allsnps1(famID): # fall snps for one family
    ln_l = 0.0
    for snpID in Q[0]:
        famL=fam1(FullData[snpID][famID],famID,snpID,delta[famID],0)
        if famL>0.0:
            ln_l+=log(famL)
        else:
            return float('-inf')
    return ln_l


def FLL():
    BigLL=0.0 # LL of full data given all latents
    for k in range(NoSnps):
        BigLL+=allfams1(SNPList[k])
    return BigLL


########## UPDATES TO parameters


def UpdateALPHA(): # 

    take=[0,0]
    ALF[0]
    nALF = ALF[0] + (random()-0.5)*alfstep
    if nALF > 0.0:
        LL0 =-ALF[0] # exponential prior with mean = 1
        LL1 =-nALF
        for j in range(1,NoFamilies+1):
            dads={}
            for offspring in range(1,len(FullData[NamedSNP][j])):
                if delta[j][offspring-1]>0:
                    dads[delta[j][offspring-1]]=1
            d=len(dads)
            if d>0:
                LL0 +=( log(gamma(ALF[0])) + float(d)*log(ALF[0]) - log(gamma(ALF[0]+float(d))) )
                LL1 +=( log(gamma(nALF)) + float(d)*log(nALF) - log(gamma(nALF+float(d))) )
    #    print ALF[0],LL0
    #    print nALF,LL1
        if LL1>=LL0: 
            ALF[0] = nALF
            take[1]+=1
        else:
            if random()< exp(LL1-LL0):
                ALF[0] = nALF
                take[1]+=1
            else:
                take[0]+=1
    #print "y ",ALF[0]
    return take



def Updatet(): # currently uniform prior

    take=[0,0]
    for spop in range(1,NoSubpops+1):

        tp= t[spop] + (random()-0.5)*tstep
        if tp<0:
            tp = -tp
        elif tp>1.0:
            tp = tp - 1.0

        cx=[0,0]
        for j in range(1,NoFamilies+1):
            if SubPop[j]==spop:
                for offspring in range(1,len(FullData[NamedSNP][j])):
                    if delta[j][offspring-1]>0:
                        cx[ 0 ]+=1
                    else:
                        cx[ 1 ]+=1

        LL0 =float(cx[0])*log(t[spop]) + float(cx[1])*log(1.0-t[spop])
        LL1 =float(cx[0])*log(tp) + float(cx[1])*log(1.0-tp)

        if LL1>=LL0: 
            t[spop]=tp
            take[1]+=1
        else:
            if random()< exp(LL1-LL0):
                t[spop]=tp
                take[1]+=1
            else:
                take[0]+=1
    return take


def UpdateQA(): # currently uniform prior; this is "ancestral Q"
    take=[0,0]
    for snpID in Q[0]:
        q=Q[0][snpID]
        qp=Q[0][snpID]+ (random()-0.5)*afstep

        if qp > AFboundary and qp < 1.0-AFboundary:
            LL0=0.0
            LL1=0.0
            for spop in range(1,NoSubpops+1):
                LL0+=log( scipy.stats.beta.pdf(Q[spop][snpID],q*(1.0-Fpop[spop])/Fpop[spop],(1-q)*(1.0-Fpop[spop])/Fpop[spop]) )
                LL1+=log( scipy.stats.beta.pdf(Q[spop][snpID],qp*(1.0-Fpop[spop])/Fpop[spop],(1-qp)*(1.0-Fpop[spop])/Fpop[spop]) )

            if LL1>=LL0: 
                Q[0][snpID]=qp
                take[1]+=1
            else:
                if random()< exp(LL1-LL0):
                    Q[0][snpID]=qp
                    take[1]+=1
                else:
                    take[0]+=1

    return take


def UpdateQsubpops(): # conditional on "ancestral Q"
    take=[0,0]
    for snpID in Q[0]:
        for spop in range(1,NoSubpops+1):
            q=Q[spop][snpID]
#             qp= beta( Q[0][snpID]*(1.0-Fpop[spop])/Fpop[spop] , (1-Q[0][snpID])*(1.0-Fpop[spop])/Fpop[spop] ) # sample new value from the prior
            qp=Q[spop][snpID]+ (random()-0.5)*afstep # use prior ratio below
            if qp > AFboundary and qp < 1.0-AFboundary:
                PriorRatio= scipy.stats.beta.pdf(qp,Q[0][snpID]*(1.0-Fpop[spop])/Fpop[spop],(1-Q[0][snpID])*(1.0-Fpop[spop])/Fpop[spop])
                PriorRatio= PriorRatio/scipy.stats.beta.pdf(q,Q[0][snpID]*(1.0-Fpop[spop])/Fpop[spop],(1-Q[0][snpID])*(1.0-Fpop[spop])/Fpop[spop])
                LL0=0.0
                LL1=0.0

                for j in range(1,NoFamilies+1):
                    if SubPop[j]==spop:
                        LL0+=log( fam1(FullData[snpID][j],j,snpID,delta[j],0) )
                        Q[spop][snpID]=qp # temp change
                        LL1+=log( fam1(FullData[snpID][j],j,snpID,delta[j],0) )
                        Q[spop][snpID]=q # revert 
                if LL1-LL0+log(PriorRatio)>0:
                    Q[spop][snpID]=qp
                    take[1]+=1
                elif random()< PriorRatio*exp(LL1-LL0):
                    Q[spop][snpID]=qp
                    take[1]+=1
                else:
                    take[0]+=1

    return take

########## UPDATES TO LATENT VARIABLES

def UpdateFpop(): # 
    take=[0,0]
    for spop in range(1,NoSubpops+1):

        f=Fpop[spop]
        fp=Fpop[spop]+ 0.02*(random()-0.5) # use prior ratio below
        if fp>0.0001 and fp<(1-0.0001):
            PriorRatio= scipy.stats.beta.pdf(fp,alphaFpop, betaFpop)/scipy.stats.beta.pdf(f,alphaFpop, betaFpop)

            LL0=0.0
            LL1=0.0
            good=0
            for snpID in Q[0]:
                try:
                    LL0+=log( scipy.stats.beta.pdf(Q[spop][snpID],Q[0][snpID]*(1.0-f)/f,  (1-Q[0][snpID])*(1.0-f)/f ) )
                    LL1+=log( scipy.stats.beta.pdf(Q[spop][snpID],Q[0][snpID]*(1.0-fp)/fp, (1-Q[0][snpID])*(1.0-fp)/fp ) )
                except ValueError:
                    take[0]+=1
                    good=1
            if good==0:
                if LL1-LL0+log(PriorRatio)>0:
                    Fpop[spop]=fp
                    take[1]+=1
                else:
                    if random() < PriorRatio*exp(LL1-LL0):
                        Fpop[spop]=fp
                        take[1]+=1
                    else:
                        take[0]+=1
    return take



def UpdateIH(): # sample proposed IH from prior (given current t)
    take=[0,0]
    cdf=[0.0 for k in range(8)]

    if IHPriorModel==1:
        cdf[0]=0.5
        for k in range(1,7):
            pk = 1.0/float(2.0+3.0*k+k*k)        
            cdf[k]=pk+cdf[k-1]



    cdf[7]=1.0


    for j in range(1,NoFamilies+1):
        ih = IH[j] 

        if IHPriorModel==0: # IH sampled based on current t
            cdf[0]=t[ SubPop[j] ]
            for k in range(1,7):
                cdf[k]= cdf[k-1] + t[ SubPop[j] ]*(1.0-t[ SubPop[j] ])**float(k)

        val =random()
        for k in range(8):
            if val<cdf[k]:
                ihp = k
                break

        if ihp != IH[j]:

            LL0=allsnps1(j)
            IH[j]=ihp # temp change
            LL1=allsnps1(j)
            IH[j]=ih # revert

            if LL1-LL0>0:
                IH[j]=ihp
                take[1]+=1
            elif random()< exp(LL1-LL0):
                IH[j]=ihp
                take[1]+=1
            else:
                take[0]+=1

    return take


def UpdateDelta(): 
    take=[0,0]
    for j in range(1,NoFamilies+1):
        for offspring in range(1,len(FullData[NamedSNP][j])):
            deltp= copy.copy(delta[j])

            Sires = max(deltp) # d in notes
            nj=[0 for x in range(Sires+1)]
            ko =0 # the number of offspring in the clutch that were produced via outcrossing
            for z in range(len(deltp)):
                nj[ deltp[z] ] +=1
                if deltp[z]>0:
                    ko+=1

            if deltp[offspring-1]==0: # current state is selfed (case 1 of table)
                priors=[1.0-t[ SubPop[j] ]]
                for x in range(1,Sires+1):
                    priors.append( t[ SubPop[j] ]*float(nj[x])/(ALF[0]+float(ko)) ) # join a current sire family
                priors.append( t[ SubPop[j] ]*ALF[0]/(ALF[0]+float(ko)) ) # make a new sire family

                LL=[]
                for x in range(len(priors)):
                    if priors[x]>0.0:
                        LL.append( log(priors[x]) )
                    else:
                        LL.append( Bigbad )
                    deltp[offspring-1]=x
                    escape=0
                    for snpID in Q[0]:
                        L1=fam1(FullData[snpID][j],j,snpID,deltp,0)
                        if L1>0.0:
                            LL[x]+=log(L1)
                        else:
                            LL[x]+=Bigbad

                pdata = 0.0
                standardizeprob= max(LL)
                for x in range(len(priors)):
                    LL[x]-=standardizeprob
                    pdata+=exp(LL[x])                
                pmod=[]
                cdf=[ 0.0 for x in range(len(priors)) ]
                for x in range(len(priors)):
                    pmod.append( exp(LL[x])/pdata )
                    if x==0:
                        cdf[0]=pmod[0]
                    else:
                        cdf[x]=cdf[x-1]+pmod[x]
                cdf[len(priors)-1]=1.0

                u9 = random() # sample a value (GIBBS sampler)
                for x in range(len(priors)):
                    if u9<cdf[x]:
                        if x != delta[j][offspring-1]:
                            delta[j][offspring-1]=x
                            take[1]+=1
                            nj[0]-=1
                            ko +=1 # the number of offspring in the clutch that were produced via outcrossing
                            if x<=Sires:    
                                nj[x]+=1    
                            else:
                                Sires+=1
                                nj.append(1)
                        else:
                            take[0]+=1
                        break

            elif deltp[offspring-1]>0 and nj[ deltp[offspring-1] ] == 1: # current state is outcrossed (case 2 of table)
                priors=[1.0-t[ SubPop[j] ]]
                for x in range(1,Sires+1):
                    if x==deltp[offspring-1]:
                        priors.append( t[ SubPop[j] ]*ALF[0]/(ALF[0]+float(ko)-1.0) )
                    else:
                        priors.append( t[ SubPop[j] ]*float(nj[x])/(ALF[0]+float(ko)-1.0) ) # join a current sire family
                LL=[]
                for x in range(len(priors)):
                    
                    if priors[x]>0.0:
                        LL.append( log(priors[x]) )
                    else:
                        LL.append( Bigbad )
                    deltp[offspring-1]=x
                    escape=0
                    for snpID in Q[0]:
                        L1=fam1(FullData[snpID][j],j,snpID,deltp,0)
                        if L1>0.0:
                            LL[x]+=log(L1)
                        else:
                            LL[x]+=Bigbad


                pdata = 0.0
                standardizeprob= max(LL)
                for x in range(len(priors)):
                    LL[x]-=standardizeprob
                    pdata+=exp(LL[x])                
                pmod=[]
                cdf=[ 0.0 for x in range(len(priors)) ]
                for x in range(len(priors)):
                    pmod.append( exp(LL[x])/pdata )
                    if x==0:
                        cdf[0]=pmod[0]
                    else:
                        cdf[x]=cdf[x-1]+pmod[x]
                cdf[len(priors)-1]=1.0

                u9 = random()
                for x in range(len(priors)):
                    if u9<cdf[x]:
                        if x != delta[j][offspring-1]:
                            whatiwas = delta[j][offspring-1]
                            delta[j][offspring-1]=x
                            take[1]+=1    
                            for z in range(len(delta[j])): # reindex
                                if delta[j][z]>whatiwas:
                                    delta[j][z]-=1

                            Sires = max(delta[j]) # d in notes
                            nj=[0 for x in range(Sires+1)]
                            ko =0 # the number of offspring in the clutch that were produced via outcrossing
                            for z in range(len(delta[j])):
                                nj[ delta[j][z] ] +=1
                                if delta[j][z]>0:
                                    ko+=1
                        else:
                            take[0]+=1
                        break


            elif deltp[offspring-1]>0 and nj[ deltp[offspring-1] ] > 1: # current state is outcrossed (case 3 of table)
                priors=[1.0-t[ SubPop[j] ]]
                for x in range(1,Sires+1):
                    if x==deltp[offspring-1]:
                        priors.append( t[ SubPop[j] ]*float(nj[x]-1)/(ALF[0]-1.0+float(ko)) )
                    else:
                        priors.append( t[ SubPop[j] ]*float(nj[x])/(ALF[0]-1.0+float(ko)) ) # join a current sire family
                priors.append( t[ SubPop[j] ]*ALF[0]/(ALF[0]-1.0+float(ko)) ) # make a new sire family

                LL=[]
                standardizeprob = 0.0
                for x in range(len(priors)):
                    if priors[x]>0.0:
                        LL.append( log(priors[x]) )
                    else:
                        LL.append( Bigbad )
                    deltp[offspring-1]=x
                    escape=0
                    for snpID in Q[0]:
                        L1=fam1(FullData[snpID][j],j,snpID,deltp,0)
                        if L1>0.0:
                            LL[x]+=log(L1)
                        else:
                            LL[x]+=Bigbad

                pdata = 0.0
                standardizeprob= max(LL)
                for x in range(len(priors)):
                    LL[x]-=standardizeprob
                    pdata+=exp(LL[x])                
                pmod=[]
                cdf=[ 0.0 for x in range(len(priors)) ]
                for x in range(len(priors)):
                    pmod.append( exp(LL[x])/pdata )
                    if x==0:
                        cdf[0]=pmod[0]
                    else:
                        cdf[x]=cdf[x-1]+pmod[x]
                cdf[len(priors)-1]=1.0

                u9 = random()
                for x in range(len(priors)):
                    if u9<cdf[x]:
                        if x != delta[j][offspring-1]:
                            nj[ delta[j][offspring-1] ]-=1
                            delta[j][offspring-1]=x
                            
                            take[1]+=1
                            if x==0:
                                ko -=1 # the number of offspring in the clutch that were produced via outcrossing
                            if x<=Sires:    
                                nj[ delta[j][offspring-1] ]+=1    
                            else:
                                Sires+=1
                                nj.append(1)
                        else:
                            take[0]+=1
                        break
            else:
                print("Should not be here, Error in UpdateDelta()")
    return take


def GenoPP(): # Posterior probabilities for genotypes
    for j in range(1,NoFamilies+1):
        famID=j
        for snpID in Q[0]:
            Data=copy.copy(FullData[snpID][j])
            delt= copy.copy(delta[j])
            Famsize=len(Data)
            PostMom={}
            PostMom[0],PostMom[1],PostMom[2]=fam1(Data,famID,snpID,delt,1)

            if PostMom[0]==-9:
                pass # return 0 # no data for family

            else:
                rx=PostMom[0]+PostMom[1]+PostMom[2]
                for ma in range(3):
                    PostMom[ma]=PostMom[ma]/rx
                    MatGPP[snpID][j][ma]+= PostMom[ma] 
                    
                    

    return 0

 
###########################

# Main program

src  =open(FILEPREFIX+".genotypes.txt", "r") # FILEPREFIX determined in Control.txt

out1a =open(rep+"."+FILEPREFIX+".Q.MAP.txt", "w")
out1b =open(rep+"."+FILEPREFIX+".Dxy.pp.txt", "w")
out2 =open(rep+"."+FILEPREFIX+".Chain.delta.txt", "w")
out2a =open(rep+"."+FILEPREFIX+".OutSelf.pp.txt", "w")
out2b =open(rep+"."+FILEPREFIX+".sibships.pp.txt", "w")
out3 =open(rep+"."+FILEPREFIX+".Chain.IH.txt", "w")
out3a =open(rep+"."+FILEPREFIX+".IH.pp.txt", "w")
out4 =open(rep+"."+FILEPREFIX+".Chain.t.txt", "w")
out4a =open(rep+"."+FILEPREFIX+".t.pp.txt", "w")
out5 =open(rep+"."+FILEPREFIX+".Chain.Fpop.txt", "w")
out5a =open(rep+"."+FILEPREFIX+".Fpop.pp.txt", "w")
out7 =open(rep+"."+FILEPREFIX+".step.history.txt", "w")


if FullOutput == 1:
    out9 =open(rep+"."+FILEPREFIX+"mat.G.pp.txt", "w")
    MatGPP={}
    PatGPP={}


ALF={0:1.0} # the concentration parameter of the Dirichlet process prior.


Fadultmean=0.0
IH={}
SubPop={}
delta={}
for j in range(1,NoFamilies+1): # set initial values for delta[j], IH of moms
    IH[j]=0    # every parent is outbred
    delta[j]=[]

inx=open(popfile, "r")
for line_idx, line in enumerate(inx):
    cols = line.replace('\n', '').split('\t') 
    if line_idx>0:
        SubPop[int(cols[0])]=int(cols[1]) # input fam, output subpop

Q={}
Fpop = {}
t = {} # initial value for population outcrossing rate
for j in range(NoSubpops+1): # Q[0] is the "Ancestral Population" 
    Q[j]={}
    if j>0:
        Fpop[j]=meanFst # initial value for among subpop fst
        t[j] =0.9 
NoSnps=0
FullData={}

SNPList=[]
vcx=0
SNPbyPop={}
for j in range(1,NoSubpops+1):
    SNPbyPop[j]={}

for line_idx, line in enumerate(src):
    cols = line.replace('\n', '').split('\t') 

    if line_idx==0:
        NamedSNP = cols[0]+"_"+cols[1]
# SNP    0    A    G
# 1    par    1.0    1.0    1.0
# 1    off    0.005    1.0    0.005
    if cols[0][0]==NamedSNP[0]:
        NoSnps+=1
        snpID=cols[0]+"_"+cols[1]    
        FullData[snpID]={}
        for j in range(1,NoSubpops+1):
            SNPbyPop[j][snpID]=[0,0.0]

        if (RUNTYPE=="0" and snpID==NamedSNP) or RUNTYPE != "0":
            if FullOutput == 1:
                MatGPP[snpID]={}
                PatGPP[snpID]={}
            SNPList.append(snpID)

            for j in range(1,NoFamilies+1):
                FullData[snpID][j]=[] 
                if FullOutput == 1:
                    MatGPP[snpID][j]=[0.0,0.0,0.0] 
    else:

        if RUNTYPE=="0" and snpID==NamedSNP:
            FullData[snpID][int(cols[0])].append([1.0,1.0,1.0]) # first list element is mom
        elif RUNTYPE != "0":
            FullData[snpID][int(cols[0])].append([float(cols[2]),float(cols[3]),float(cols[4])]) # first list element is mom (mother must precede offspring in each family)
            if float(cols[2])*float(cols[3])*float(cols[4])<0.99: # data at SNP
                denom = 0.25*float(cols[2])+ 0.5*float(cols[3])+ 0.25*float(cols[4]) # uniform prior
                pr = ( 0.25*float(cols[2])+ 0.25*float(cols[3]) )/denom
                SNPbyPop[ SubPop[int(cols[0])] ][snpID][0]+=1
                SNPbyPop[ SubPop[int(cols[0])] ][snpID][1]+=pr

        if snpID==NamedSNP:
            if cols[1]=="off":
                sirenum=len(delta[int(cols[0])])+1
                delta[int(cols[0])].append(sirenum) # everybody is outcrossed half sib at beginning of run
if RUNTYPE=="0":
    NoSnps=1


for k in range(len(SNPList)):
    snpID=SNPList[k]
    qlist=[]
    for j in range(1,NoSubpops+1):
        if SNPbyPop[ j ][snpID][0]>0:
            Q[j][snpID]=SNPbyPop[ j ][snpID][1]/float(SNPbyPop[ j ][snpID][0])
            if Q[j][snpID]<AFboundary:
                Q[j][snpID]=AFboundary
            elif 1.0-Q[j][snpID]<AFboundary:
                Q[j][snpID]=1.0-AFboundary
            qlist.append(Q[j][snpID])
        else:
            Q[j][snpID]=0.5
    Q[0][snpID]=sum(qlist)/float(len(qlist))

print("Number of Snps ",NoSnps)
if FullOutput == 1:
    for snpID in Q[0]:
        for j in range(1,NoFamilies+1):
            Famsize=len(FullData[snpID][j])
            if Famsize>1:
                PatGPP[snpID][j]={}
                for offspring in range(1,Famsize):
                    PatGPP[snpID][j][offspring-1]=[0.0,0.0]
Chainsibs={}
ChainDelta={}
ChainIH={}
Chaint={}
ChainQ={}
ChainDxy={}
ChainFpop={}

for j in range(NoSubpops+1):
    ChainQ[j]={}
    ChainDxy[j]={}
    for k in range(j,NoSubpops+1):
        ChainDxy[j][k]={}
        for k2 in range(1001):
            ChainDxy[j][k][k2]=0

    if j>0:
        ChainFpop[j]={}
        Chaint[j]={}
        for k in range(101):
            ChainFpop[j][k]=0
            Chaint[j][k]=0

    for s in Q[0]: # cycle through all snps
        ChainQ[j][s]={}
        for k in range(101):
            ChainQ[j][s][k]=0


for j in range(1,NoFamilies+1):
    ChainIH[j]=[0 for x in range(8)]
    for offspring in range(1,len(FullData[NamedSNP][j])):
        ChainDelta[str(j)+"_"+str(offspring)]=[0,0]
        for zoff in range(offspring+1,len(FullData[NamedSNP][j])):
            Chainsibs[str(j)+"_"+str(offspring)+"_"+str(zoff)]=[0,0]


for steps in range(TuningSteps): # tune in allele freqs
    print("Tuning allele freqs step ",steps,"Current LL",round(FLL(),2))
    stay,go= UpdateFpop()
    stay,go= UpdateQA()
    stay,go= UpdateQsubpops()


movehistory=[[0,0] for zy in range(10)]
for steps in range(ChainLength):

    if steps % thinningfreq == 0:

# moved updates here
        stayD,goD = UpdateDelta()
        movehistory[2][0]+=stayD
        movehistory[2][1]+=goD
        stay,go= UpdateIH()
        movehistory[3][0]+=stay
        movehistory[3][1]+=go
        Fadultmean=0.0
        for j in range(1,NoFamilies+1):
            Fadultmean+= Fmom[ IH[j] ]/float(NoFamilies)

        t[0]=0.0
        for j in range(1,NoSubpops+1):
            Chaint[j][int(t[j]*100+0.0049)]+=1
            t[0]+= t[j]/float(NoSubpops)        
            cLL=FLL()
            out4.write(str(steps)+'\t'+str(t[0])+'\t'+str(ALF[0])+'\t'+str(cLL)+'\n')

        if steps>=burnin:
            dxytemp={}
            for j in range(NoSubpops+1):
                dxytemp[j]={}
                for k in range(j,NoSubpops+1):
                    dxytemp[j][k]=[]
            for snpid in Q[0]:
                for j in range(1,NoSubpops+1):
                    ChainQ[j][snpid][int(Q[j][snpid]*100+0.0049)]+=1
                    for k in range(j,NoSubpops+1):
                        if SNPbyPop[ j ][snpid][0]>=MinCdxy and SNPbyPop[ k ][snpid][0]>=MinCdxy:
                            dxytemp[j][k].append( Q[j][snpid]+Q[k][snpid]-2.0*Q[j][snpid]*Q[k][snpid] )

            for j in range(1,NoSubpops+1):
                for k in range(j,NoSubpops+1):
                    val = sum(dxytemp[j][k])/float(len(dxytemp[j][k]))
                    # print j,k,len(dxytemp[j][k]),val
                    ChainDxy[j][k][int(val*1000+0.00049)]+=1

                    
            for j in range(1,NoFamilies+1):
                out3.write(str(steps)+'\t'+str(j)+'\t'+str(IH[j])+'\n')
                ChainIH[j][ IH[j] ]+=1
                for offspring in range(1,len(FullData[NamedSNP][j])):
                    out2.write(str(steps)+'\t'+str(j)+'\t'+str(offspring)+'\t'+str(delta[j][offspring-1])+'\n')
                    if delta[j][offspring-1]>0:
                        ChainDelta[str(j)+"_"+str(offspring)][ 0 ]+=1
                        for zoff in range(offspring+1,len(FullData[NamedSNP][j])):
                            if delta[j][zoff-1]>0:
                                if delta[j][offspring-1]==delta[j][zoff-1]: # fullsibs
                                    Chainsibs[str(j)+"_"+str(offspring)+"_"+str(zoff)][0]+=1
                                else: # halfsibs
                                    Chainsibs[str(j)+"_"+str(offspring)+"_"+str(zoff)][1]+=1
                    else:
                        ChainDelta[str(j)+"_"+str(offspring)][ 1 ]+=1
            out5.write(str(steps))
            for j in range(1,NoSubpops+1):
                out5.write('\t'+str(Fpop[j]))
                ChainFpop[j][int(Fpop[j]*100+0.0049)]+=1
            out5.write('\n')

            if FullOutput == 1:
                nothingnumber = GenoPP()

        print("Step ",steps," current LL ",round(cLL,2)," t0,etc ",t)

#####    UPDATE parameters
    stay,go= UpdateFpop()
    movehistory[0][0]+=stay
    movehistory[0][1]+=go
    stay,go= Updatet()
    movehistory[1][0]+=stay
    movehistory[1][1]+=go

    stay,go= UpdateQA()
    movehistory[4][0]+=stay
    movehistory[4][1]+=go
    stay,go= UpdateQsubpops()
    movehistory[5][0]+=stay
    movehistory[5][1]+=go
    stay,go= UpdateALPHA()
    movehistory[6][0]+=stay
    movehistory[6][1]+=go

#output summaries
out2a.write('pop\tfamily\toffspring\tout\tself\n')
out2b.write('pop\tfamily\toff1\toff2\tfull\thalf\n')
out3a.write('pop\tfamily\tIH=0\tIH=1\tIH=2\tIH=3\tIH=4\tIH=5\tIH=6\tIH=7+\n')
for j in range(1,NoFamilies+1):
    out3a.write(str(SubPop[j])+'\t'+str(j))
    for z in range(8):
        out3a.write('\t'+str(ChainIH[j][z]))
    out3a.write('\n')
    for offspring in range(1,len(FullData[NamedSNP][j])):
        out2a.write(str(SubPop[j])+'\t'+str(j)+'\t'+str(offspring)+'\t'+str(ChainDelta[str(j)+"_"+str(offspring)][0])+'\t'+str(ChainDelta[str(j)+"_"+str(offspring)][1])+'\n')
        for zoff in range(offspring+1,len(FullData[NamedSNP][j])):
            out2b.write(str(SubPop[j])+'\t'+str(j)+'\t'+str(offspring)+'\t'+str(zoff)+'\t'+str(Chainsibs[str(j)+"_"+str(offspring)+"_"+str(zoff)][0])+'\t'+str(Chainsibs[str(j)+"_"+str(offspring)+"_"+str(zoff)][1])+'\n')

out4a.write('t value\tfrequency\n')
out5a.write('Subpop\tvalue\tfrequency\n')
for k in range(1, NoSubpops+1):
    for j in range(101):
        out5a.write(str(k)+'\t'+str(float(j)/100 + 0.005)+'\t'+str(ChainFpop[k][j])+'\n')
        out4a.write(str(k)+'\t'+str(float(j)/100 + 0.005)+'\t'+str(Chaint[k][j])+'\n')

out7.write("parameter    rejected    accepted    ARate\n")
# stx=["Fpop","t","delta","IH","Qancestral","Qsubpop","BetaGu"] # stx[j]+'\t'+
for j in range(7):
    out7.write(str(movehistory[j][0])+'\t'+str(movehistory[j][1])+'\t'+str(float(movehistory[j][1])/float(movehistory[j][1]+movehistory[j][0]))+'\n') 

for snpid in Q[0]:
    out1a.write(snpid)
    for j in range(NoSubpops+1):
        mx=0.0
        for k in range(101):
            if ChainQ[j][snpid][k]>mx:
                mx = ChainQ[j][snpid][k]
                val = k
        out1a.write('\t'+str(float(val)/100 + 0.005) )

    out1a.write('\n')


for j in range(1,NoSubpops+1):
    for k in range(j,NoSubpops+1):

        meanx=[0.0,0.0]    
        for k2 in range(1001):
            meanx[0]+=float(ChainDxy[j][k][k2])    
            meanx[1]+=float( ChainDxy[j][k][k2]*float(k2)/1000.0 )
        out1b.write(str(j)+'\t'+str(k)+'\t'+str(meanx[1]/meanx[0])+'\n' )    


if FullOutput == 1:
    for snpID in Q[0]:
        for j in range(1,NoFamilies+1):
            rx=MatGPP[snpID][j][0]+MatGPP[snpID][j][1]+MatGPP[snpID][j][2]
            if rx>0.0:
                out9.write(snpID+'\t'+str(j)+'\t'+str(MatGPP[snpID][j][0]/rx)+'\t'+str(MatGPP[snpID][j][1]/rx)+'\t'+str(MatGPP[snpID][j][2]/rx)+'\n')



