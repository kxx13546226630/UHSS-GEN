      SUBROUTINE UMAT(STRESS,STATEV,DDSDDE,SSE,SPD,SCD,
     1 RPL,DDSDDT,DRPLDE,DRPLDT,
     2 STRAN,DSTRAN,TIME,DTIME,TEMP,DTEMP,PREDEF,DPRED,CMNAME,
     3 NDI,NSHR,NTENS,NSTATV,PROPS,NPROPS,COORDS,DROT,PNEWDT,
     4 CELENT,DFGRD0,DFGRD1,NOEL,NPT,LAYER,KSPT,JSTEP,KINC)
C
      INCLUDE 'ABA_PARAM.INC'
C
      CHARACTER*80 CMNAME
      DIMENSION STRESS(NTENS),STATEV(NSTATV),
     1 DDSDDE(NTENS,NTENS),DDSDDT(NTENS),DRPLDE(NTENS),
     2 STRAN(NTENS),DSTRAN(NTENS),TIME(2),PREDEF(1),DPRED(1),
     3 PROPS(NPROPS),COORDS(3),DROT(3,3),DFGRD0(3,3),DFGRD1(3,3),
     4 JSTEP(4)
      
      DIMENSION DE1(6,6),DE2(6,6),DE(6,6),S(6),S1(6),S2(6)
      DIMENSION DSTRESS(6),DSTRESS1(6),DSTRESS2(6),STRESS1(6),STRESS2(6)
      DIMENSION FG1(6),FG2(6),FG(6),PS(3),SS(6),SS1(6),SS2(6)
      DIMENSION DEP1(6,6),DEP2(6,6)
      DIMENSION DDSTRAN(6),ESTRESS(6)
    
      REAL R1,MF,I1,I2,I3,I11,I22,I33
	  
      SSTOL=1E-3
       
      FLAMA=PROPS(1)   
      FKAPA=PROPS(2)   
      FU=PROPS(3)      
      FM=PROPS(4)      
	  OCR=PROPS(5)     
      FVOID0=PROPS(6)  
	  FGref=PROPS(7)   
      FGamaref=PROPS(8)
	  FGama0=PROPS(9)
	  FAFA=PROPS(10)    
	  FN=PROPS(11)     
	  FS0=PROPS(12)    
 
      FTIME=0.0
      FDTIME=1.0   

 888   CONTINUE

	  FPC0=STATEV(1)
      EVP=STATEV(2) 
	  P=STATEV(3)


      CALL SINV(STRESS,SINV1,SINV2,NDI,NSHR) 
      FP1=SINV1                              
      FQ1=SINV2                               
	  IF(FQ1 .LE. 1.0e-5) FQ1=1.0e-5          

       IF(TIME(2).EQ.0.0) then
	  P=abs(FP1)
	  FPC0=abs(OCR*FP1)	  
	  endif	
	  
	  CALL SPRINC(STRAN,PS,2,NDI,NSHR)
	  GAMAS=SQRT(0.5*((PS(1)-PS(2))**2+(PS(1)-PS(3))**2+(PS(2)-PS(3))**2))
	  

	  IF (GAMAS.LT.FGama0) THEN
	     FGMOD=FGref*OCR**0.2*(P/100000.0)**0.8/(1+3.0/7.0*GAMAS/FGamaref)
	     FKMOD=FGMOD*2.0*(1.0+FU)/3.0/(1.0-2.0*FU)
         CALL GETDE(FKMOD,FGMOD,DDSDDE)    
       CALL GETDSTRESS(DDSDDE,DSTRAN,DSTRESS) 
       DO I=1,6
           STRESS(I)=STRESS(I)+DSTRESS(I) 
      END DO	  
	  ELSE


     
      FP3=FP1                              
	  a=1
	  FKMOD1=(1.0+FVOID0)*abs(FP3)/FKAPA           
      FGMOD1=FKMOD1*3.0*(1.0-2.0*FU)/2.0/(1.0+FU) 
      CALL GETDE(FKMOD1,FGMOD1,DE1)           
	  
		DO I=1,3
		S(I)=STRESS(I)-(-FP3+FS0)**FN/(100000.0)**(FN-1)-FP3
        END DO
		DO I=4,6
		S(I)=STRESS(I)
        END DO	
	
      CALL SINV(S,SINV1,SINV2,NDI,NSHR) 
      FP4=SINV1                               
      FQ4=SINV2                               
	  IF(FQ4.LE.1.0e-5) FQ4=1.0e-5             
 	  CALL SPRINC(S,PS,1,NDI,NSHR)
	  I1=PS(1)+PS(2)+PS(3)
	  I2=PS(1)*PS(2)+PS(1)*PS(3)+PS(2)*PS(3)	  
	  I3=PS(1)*PS(2)*PS(3)	  
	  

	  FQSTAR1=FAFA*(ABS(I1**2-3*I2))**0.5-(1-FAFA)*(2*I1)
     1/(3*(ABS((I1*I2-I3)/(I1*I2-9*I3)))**0.5-1)
		DO I=1,3
		SS(I)=FP4+FQSTAR1/FQ4*(S(I)-FP4)
        END DO
		DO I=4,6
		SS(I)=FQSTAR1/FQ4*S(I)
        END DO
      CALL SINV(SS,SINV1,SINV2,NDI,NSHR)
      FP1=SINV1
      FQ1=SINV2
      IF(FQ1 .LE. 1.0e-5) FQ1=1.0e-5               
	  
      FATA=FQ1/FP1 
	  CHI=1.0*FM*FM/12.0/(3.0-FM)
      R1=1.0*abs(FP1)/FPC0*(1.0+FATA**2/FM**2)
     1 *exp(EVP*(1.0+FVOID0)/(FLAMA-FKAPA))

	  MF=6.0*(SQRT(CHI/R1*(1.0+CHI/R1))-CHI/R1)
	  
      DFDP=(FM*FM-FATA*FATA)/FP1/(FM*FM+FATA*FATA)     
      DFDQ=2.0*FATA/FP1/(FM*FM+FATA*FATA)              
      DGDP=(FM*FM-FATA*FATA)/FP1/(FM*FM+FATA*FATA)      
      DGDQ=2.0*FATA/FP1/(FM*FM+FATA*FATA)                 
	  FPFH=-(1.0+FVOID0)/(FLAMA-FKAPA)
	  FKP1=FPFH*(MF**4-FATA**4)/(FM**2+FATA**2)**2/FP1
   
 999   CONTINUE    
      
      DO I=1,6
          DDSTRAN(I)=DSTRAN(I)*FDTIME  
      END DO

	 
      CALL GETDEP(I1,I2,I3,FAFA,S,FN,FS0,FP3,FP4,FQ4,FQ1,FQSTAR1,
     1 SS,DGDP,DGDQ,DFDP,DFDQ,DE1,FKP1,DEP1,DDSTRAN,DEVP1,FF,FG,FR)                 
	 
	  EVP1=EVP+DEVP1 
	 
      CALL GETDSTRESS(DEP1,DDSTRAN,DSTRESS1) 
 
      DO I=1,6
          STRESS1(I)=STRESS(I)+DSTRESS1(I) 
      END DO


      CALL SINV(STRESS1,SINV1,SINV2,NDI,NSHR) 
      FP33=SINV1                              
      FKMOD2=(1+FVOID0)*abs(FP33)/FKAPA
      FGMOD2=FKMOD2*3.0*(1.0-2.0*FU)/2.0/(1.0+FU)
      CALL GETDE(FKMOD2,FGMOD2,DE2)		  
		DO I=1,3
		S1(I)=STRESS1(I)-(-FP33+FS0)**FN/(100000.0)**(FN-1)-FP33
        END DO
		DO I=4,6
		S1(I)=STRESS1(I)
        END DO	
	
      CALL SINV(S1,SINV1,SINV2,NDI,NSHR) 
      FP44=SINV1                              
      FQ44=SINV2                           
	  IF(FQ44 .LE. 1.0e-5) FQ44=1.0e-5         
 	  CALL SPRINC(S1,PS,1,NDI,NSHR)
	  I11=PS(1)+PS(2)+PS(3)
	  I22=PS(1)*PS(2)+PS(1)*PS(3)+PS(2)*PS(3)	  
	  I33=PS(1)*PS(2)*PS(3)	
	  

	  FQSTAR2=FAFA*(ABS(I11**2-3*I22))**0.5-(1-FAFA)*(2*I11)
     1/(3*(ABS((I11*I22-I33)/(I11*I22-9*I33)))**0.5-1)
		DO I=1,3
		SS1(I)=FP44+FQSTAR2/FQ44*(S1(I)-FP44)
        END DO
		DO I=4,6
		SS1(I)=FQSTAR2/FQ44*S1(I)
        END DO



  
      CALL SINV(SS1,SINV1,SINV2,NDI,NSHR)
      FP2=SINV1
      FQ2=SINV2
      IF(FQ2 .LE. 1.0e-5) FQ2=1.0e-5         
      FATA=FQ2/FP2   

      R1=1.0*abs(FP2)/FPC0*(1.0+FATA**2/FM**2)
     1 *exp(EVP1*(1.0+FVOID0)/(FLAMA-FKAPA))
	  MF=6.0*(SQRT(CHI/R1*(1.0+CHI/R1))-CHI/R1)
	  
     
      DFDP2=(FM*FM-FATA*FATA)/FP2/(FM*FM+FATA*FATA)     
      DFDQ2=2.0*FATA/FP2/(FM*FM+FATA*FATA)                 
      DGDP2=(FM*FM-FATA*FATA)/FP2/(FM*FM+FATA*FATA)       
      DGDQ2=2.0*FATA/FP2/(FM*FM+FATA*FATA)                
  
	  FPFH=-(1.0+FVOID0)/(FLAMA-FKAPA)
	  FKP2=FPFH*(MF**4-FATA**4)/(FM**2+FATA**2)**2/FP2


      CALL GETDEP(I11,I22,I33,FAFA,S1,FN,FS0,FP33,FP44,FQ44,FQ2,FQSTAR2,
     1 SS1,DGDP2,DGDQ2,DFDP2,DFDQ2,DE2,FKP2,DEP2,DDSTRAN,DEVP2,FF,FG,FR)               

      CALL GETDSTRESS(DEP2,DDSTRAN,DSTRESS2)
    
      DO I=1,6
          ESTRESS(I)=0.5*(DSTRESS2(I)-DSTRESS1(I))
          STRESS2(I)=STRESS(I)+0.5*DSTRESS1(I)+0.5*DSTRESS2(I)
      END DO
     
      FEIE=0.0
      FEIS=0.0
      FERR=0.0
      DO I=1,6
          FEIE=FEIE+ESTRESS(I)*ESTRESS(I)
          FEIS=FEIS+STRESS2(I)*STRESS2(I)
      END DO
      FERR=SQRT(abs(FEIE/FEIS))
      IF(FERR.LE.1E-8)FERR=1E-8
      FBETA=0.8*SQRT(abs(SSTOL/FERR))

      IF(FERR.GT.SSTOL)THEN
          
          IF(FBETA.LE.0.1)FBETA=0.1
          FDTIME=FBETA*FDTIME
          GOTO 999
      ELSE
          FTIME=FTIME+FDTIME
 
          IF(FBETA.GE.2.0)FBETA=2.0
          FDTIME=FBETA*FDTIME
      END IF
C-------误差标准满足后计算应力--------------
      
      DO I=1,6
          STRESS(I)=STRESS(I)+0.5*DSTRESS1(I)+0.5*DSTRESS2(I)
      END DO

      DEVP=0.5*(DEVP1+DEVP2)
	  EVP=EVP+DEVP
	  
	  STATEV(1)=FPC0
      STATEV(2)=EVP
      STATEV(3)=P


	  
      IF(FTIME.LT.1.0)THEN
        IF(FDTIME.GT.(1.0-FTIME))THEN
          FDTIME=1.0-FTIME
        END IF
      GOTO 888
      END IF


C过渡应力空间
      CALL SINV(STRESS,SINV1,SINV2,NDI,NSHR) 
        FP3=SINV1                              
      FKMOD=(1.0+FVOID0)*abs(FP3)/FKAPA
      FGMOD=FKMOD*3.0*(1.0-2.0*FU)/2.0/(1.0+FU)
      CALL GETDE(FKMOD,FGMOD,DE)
		DO I=1,3
		S(I)=STRESS(I)-(-FP3+FS0)**FN/(100000.0)**(FN-1)-FP3
        END DO
		DO I=4,6
		S(I)=STRESS(I)
        END DO	
C转换应力空间	
      CALL SINV(S,SINV1,SINV2,NDI,NSHR)
      FP4=SINV1                               
      FQ4=SINV2                               
	  IF(FQ4 .LE. 1.0e-5) FQ4=1.0e-5          
 	  CALL SPRINC(S,PS,1,NDI,NSHR)
	  I1=PS(1)+PS(2)+PS(3)
	  I2=PS(1)*PS(2)+PS(1)*PS(3)+PS(2)*PS(3)	  
	  I3=PS(1)*PS(2)*PS(3)	
	  

	  FQSTAR=FAFA*(ABS(I1**2-3*I2))**0.5-(1-FAFA)*(2*I1)
     1/(3*(ABS((I1*I2-I3)/(I1*I2-9*I3)))**0.5-1)
		DO I=1,3
		SS(I)=FP4+FQSTAR/FQ4*(S(I)-FP4)
        END DO
		DO I=4,6
		SS(I)=FQSTAR/FQ4*S(I)
        END DO


      CALL SINV(SS,SINV1,SINV2,NDI,NSHR)
      FP=SINV1
      FQ=SINV2
      IF(FQ .LE. 1.0e-5) FQ=1.0e-5              
  
      FATA=FQ/FP
      R1=1.0*abs(FP)/FPC0*(1.0+FATA**2/FM**2)
     1*exp(EVP*(1.0+FVOID0)/(FLAMA-FKAPA))
	  MF=6.0*(SQRT(CHI/R1*(1.0+CHI/R1))-CHI/R1)
	  
      DFDP=(FM*FM-FATA*FATA)/FP/(FM*FM+FATA*FATA)      
      DFDQ=2.0*FATA/FP/(FM*FM+FATA*FATA)                
      DGDP=(FM*FM-FATA*FATA)/FP/(FM*FM+FATA*FATA)       
      DGDQ=2.0*FATA/FP/(FM*FM+FATA*FATA)                
	  
	  FPFH=-(1.0+FVOID0)/(FLAMA-FKAPA)
	  FKP=FPFH*(MF**4-FATA**4)/(FM**2+FATA**2)**2/FP
 

      CALL GETDEP(I1,I2,I3,FAFA,S,FN,FS0,FP3,FP4,FQ4,FQ,FQSTAR,
     1 SS,DGDP,DGDQ,DFDP,DFDQ,DE,FKP,DDSDDE,DSTRAN,DEVP,FF,FG,FR)                  
	  END IF

  
	  STATEV(1)=FPC0
      STATEV(2)=EVP	  
      STATEV(3)=P
      STATEV(4)=MF	  
      STATEV(5)=R1
      STATEV(6)=a
	  
      RETURN
      END
	  
     
      SUBROUTINE GETDE(FKMOD,FGMOD,FDE) 
      INCLUDE 'ABA_PARAM.INC'
      DIMENSION FDE(6,6)
      DO I=1,6
          DO J=1,6
              FDE(I,J)=0.0
          END DO
      END DO
 
      FDE(1,1)=FKMOD+4.0/3.0*FGMOD
      FDE(2,2)=FKMOD+4.0/3.0*FGMOD
      FDE(3,3)=FKMOD+4.0/3.0*FGMOD
      FDE(4,4)=FGMOD
      FDE(5,5)=FGMOD
      FDE(6,6)=FGMOD
      FDE(1,2)=FKMOD-2.0/3.0*FGMOD
      FDE(1,3)=FKMOD-2.0/3.0*FGMOD
      FDE(2,1)=FKMOD-2.0/3.0*FGMOD
      FDE(2,3)=FKMOD-2.0/3.0*FGMOD
      FDE(3,1)=FKMOD-2.0/3.0*FGMOD
      FDE(3,2)=FKMOD-2.0/3.0*FGMOD
      END 
      
      SUBROUTINE GETDSTRESS(FDE,DER,DSTRESS)
      INCLUDE 'ABA_PARAM.INC'
      DIMENSION FDE(6,6),DER(6),DSTRESS(6)
     
      DO I=1,6
          DSTRESS(I)=0.0
      END DO
      DO I=1,6
          DO J=1,6
              DSTRESS(I)=DSTRESS(I)+FDE(I,J)*DER(J)
          END DO
      END DO
      END
      
      SUBROUTINE GETDEP(I1,I2,I3,FAFA,FS,FN,FS0,FP3,FP4,FQ4,FQ,FQSTAR,
     1 FSTRESS,DGDP,DGDQ,DFDP,DFDQ,FDE,FKP,DEP,DER,DEVP,FF,FG,FR)       
      INCLUDE 'ABA_PARAM.INC'
      DIMENSION FSTRESS(6),FDE(6,6),DEP(6,6),DER(6)
      DIMENSION FA(6),FB(6),FC(6),FD(6),FE(6),FG(6),FH(6,6)
      DIMENSION DS(6,6),DI1(6),DI2(6),DI3(6)   
      DIMENSION FD1(6),FE1(6),DQSTARDS(6),FB1(6)
      DIMENSION DSS(6,6),FS(6)
      REAL I1,I2,I3	  
      FA(1)=1.0/3.0 
      FA(2)=1.0/3.0
      FA(3)=1.0/3.0
      FA(4)=0.0
      FA(5)=0.0
      FA(6)=0.0

      FB(1)=1.0/2.0/FQ*(2*FSTRESS(1)-FSTRESS(2)-FSTRESS(3))  
      FB(2)=1.0/2.0/FQ*(2*FSTRESS(2)-FSTRESS(1)-FSTRESS(3))
      FB(3)=1.0/2.0/FQ*(2*FSTRESS(3)-FSTRESS(1)-FSTRESS(2))
      FB(4)=3.0/2.0/FQ*2*FSTRESS(4)
      FB(5)=3.0/2.0/FQ*2*FSTRESS(5)
      FB(6)=3.0/2.0/FQ*2*FSTRESS(6)

      FB1(1)=1.0/2.0/FQ4*(2*FS(1)-FS(2)-FS(3)) 
      FB1(2)=1.0/2.0/FQ4*(2*FS(2)-FS(1)-FS(3))
      FB1(3)=1.0/2.0/FQ4*(2*FS(3)-FS(1)-FS(2))
      FB1(4)=3.0/2.0/FQ4*2*FS(4)
      FB1(5)=3.0/2.0/FQ4*2*FS(5)
      FB1(6)=3.0/2.0/FQ4*2*FS(6)
  
      DS(1,1)=2.0/3.0+FN*(-FP3+FS0)**(FN-1)/3/(100000.0)**(FN-1)  
      DS(1,2)=-1.0/3.0+FN*(-FP3+FS0)**(FN-1)/3/(100000.0)**(FN-1)
      DS(1,3)=-1.0/3.0+FN*(-FP3+FS0)**(FN-1)/3/(100000.0)**(FN-1)
      DS(1,4)=0.0
      DS(1,5)=0.0
      DS(1,6)=0.0	  

      DS(2,1)=-1.0/3.0+FN*(-FP3+FS0)**(FN-1)/3/(100000.0)**(FN-1) 
      DS(2,2)=2.0/3.0+FN*(-FP3+FS0)**(FN-1)/3/(100000.0)**(FN-1)
      DS(2,3)=-1.0/3.0+FN*(-FP3+FS0)**(FN-1)/3/(100000.0)**(FN-1)
      DS(2,4)=0.0
      DS(2,5)=0.0
      DS(2,6)=0.0

      DS(3,1)=-1.0/3.0+FN*(-FP3+FS0)**(FN-1)/3/(100000.0)**(FN-1) 
      DS(3,2)=-1.0/3.0+FN*(-FP3+FS0)**(FN-1)/3/(100000.0)**(FN-1)
      DS(3,3)=2.0/3.0+FN*(-FP3+FS0)**(FN-1)/3/(100000.0)**(FN-1)
      DS(3,4)=0.0
      DS(3,5)=0.0
      DS(3,6)=0.0

      DS(4,1)=0.0  
      DS(4,2)=0.0
      DS(4,3)=0.0
      DS(4,4)=1.0
      DS(4,5)=0.0
      DS(4,6)=0.0

      DS(5,1)=0.0
      DS(5,2)=0.0
      DS(5,3)=0.0
      DS(5,4)=0.0
      DS(5,5)=1.0
      DS(5,6)=0.0

      DS(6,1)=0.0 
      DS(6,2)=0.0
      DS(6,3)=0.0
      DS(6,4)=0.0
      DS(6,5)=0.0
      DS(6,6)=1.0
	  
      DI1(1)=1.0  
      DI1(2)=1.0
      DI1(3)=1.0
      DI1(4)=0.0
      DI1(5)=0.0
      DI1(6)=0.0

      DI2(1)=FS(2)+FS(3) 
      DI2(2)=FS(1)+FS(3)
      DI2(3)=FS(2)+FS(1)
      DI2(4)=-2.0*FS(4)
      DI2(5)=-2.0*FS(5)
      DI2(6)=-2.0*FS(6)

      DI3(1)=FS(2)*FS(3)-FS(5)**2  
      DI3(2)=FS(3)*FS(1)-FS(6)**2
      DI3(3)=FS(1)*FS(2)-FS(4)**2
      DI3(4)=2.0*FS(5)*FS(6)-2.0*FS(3)*FS(4)
      DI3(5)=2.0*FS(4)*FS(6)-2.0*FS(1)*FS(5)
      DI3(6)=2.0*FS(5)*FS(4)-2.0*FS(2)*FS(6)

	  DQSTARDI1=FAFA*I1/SQRT(ABS(I1**2-3.0*I2))-(1.0-FAFA)
     1*((6.0*(ABS((I1*I2-I3)/(I1*I2-9*I3)))**0.5-2)+24.0*I1*I2*I3/
     2(I1*I2-9*I3)**2/(ABS((I1*I2-I3)/(I1*I2-9*I3)))**0.5)
     1/(3.0*((I1*I2-I3)/(I1*I2-9*I3))**0.5-1)**2

	  DQSTARDI2=-3.0/2.0*FAFA/SQRT(ABS(I1**2-3.0*I2))
     1   -(1.0-FAFA)*24.0*I1*I1*I3	 
     1/(I1*I2-9*I3)**2/(ABS((I1*I2-I3)/(I1*I2-9*I3)))**0.5	  
     2/(3.0*(ABS((I1*I2-I3)/(I1*I2-9*I3)))**0.5-1)**2
	 
	  DQSTARDI3=(1.0-FAFA)*24.0*I1*I1*I2	 
     1/(I1*I2-9*I3)**2/(ABS((I1*I2-I3)/(I1*I2-9*I3)))**0.5	  
     2/(3.0*(ABS((I1*I2-I3)/(I1*I2-9*I3)))**0.5-1)**2


      DO I=1,6	 
	  DQSTARDS(I)=DQSTARDI1*DI1(I)+DQSTARDI2*DI2(I)+DQSTARDI3*DI3(I)  
      END DO
	  
      DSS(1,1)=1.0/3.0+2.0/3.0*FQSTAR/FQ4+(FS(1)-FP4)
     1*(DQSTARDS(1)*FQ4-FB1(1)*FQSTAR)/FQ4**2  
      DSS(1,2)=1.0/3.0-1.0/3.0*FQSTAR/FQ4+(FS(1)-FP4)
     1*(DQSTARDS(2)*FQ4-FB1(2)*FQSTAR)/FQ4**2  
      DSS(1,3)=1.0/3.0-1.0/3.0*FQSTAR/FQ4+(FS(1)-FP4)
     1*(DQSTARDS(3)*FQ4-FB1(3)*FQSTAR)/FQ4**2 
      DSS(1,4)=(FS(1)-FP4)*(DQSTARDS(4)*FQ4-FB1(4)*FQSTAR)
     1/FQ4**2 
      DSS(1,5)=(FS(1)-FP4)*(DQSTARDS(5)*FQ4-FB1(5)*FQSTAR)
     1/FQ4**2  
      DSS(1,6)=(FS(1)-FP4)*(DQSTARDS(6)*FQ4-FB1(6)*FQSTAR)
     1/FQ4**2 

      DSS(2,1)=1.0/3.0-1.0/3.0*FQSTAR/FQ4+(FS(2)-FP4)
     1*(DQSTARDS(1)*FQ4-FB1(1)*FQSTAR)/FQ4**2  
      DSS(2,2)=1.0/3.0+2.0/3.0*FQSTAR/FQ4+(FS(2)-FP4)
     1*(DQSTARDS(2)*FQ4-FB1(2)*FQSTAR)/FQ4**2 
      DSS(2,3)=1.0/3.0-1.0/3.0*FQSTAR/FQ4+(FS(2)-FP4)
     1*(DQSTARDS(3)*FQ4-FB1(3)*FQSTAR)/FQ4**2  
      DSS(2,4)=(FS(2)-FP4)*(DQSTARDS(4)*FQ4-FB1(4)*FQSTAR)
     1/FQ4**2  
      DSS(2,5)=(FS(2)-FP4)*(DQSTARDS(5)*FQ4-FB1(5)*FQSTAR)
     1/FQ4**2 
      DSS(2,6)=(FS(2)-FP4)*(DQSTARDS(6)*FQ4-FB1(6)*FQSTAR)
     1/FQ4**2 

      DSS(3,1)=1.0/3.0-1.0/3.0*FQSTAR/FQ4+(FS(3)-FP4)
     1*(DQSTARDS(1)*FQ4-FB1(1)*FQSTAR)/FQ4**2  
      DSS(3,2)=1.0/3.0-1.0/3.0*FQSTAR/FQ4+(FS(3)-FP4)
     1*(DQSTARDS(2)*FQ4-FB1(2)*FQSTAR)/FQ4**2  
      DSS(3,3)=1.0/3.0+2.0/3.0*FQSTAR/FQ4+(FS(3)-FP4)
     1*(DQSTARDS(3)*FQ4-FB1(3)*FQSTAR)/FQ4**2  
      DSS(3,4)=(FS(3)-FP4)*(DQSTARDS(4)*FQ4-FB1(4)*FQSTAR)
     1/FQ4**2  
      DSS(3,5)=(FS(3)-FP4)*(DQSTARDS(5)*FQ4-FB1(5)*FQSTAR)
     1/FQ4**2  
      DSS(3,6)=(FS(3)-FP4)*(DQSTARDS(6)*FQ4-FB1(6)*FQSTAR)
     1/FQ4**2 

      DSS(4,1)=FS(4)*(DQSTARDS(1)*FQ4-FB1(1)*FQSTAR)/FQ4**2  
      DSS(4,2)=FS(4)*(DQSTARDS(2)*FQ4-FB1(2)*FQSTAR)/FQ4**2  
      DSS(4,3)=FS(4)*(DQSTARDS(3)*FQ4-FB1(3)*FQSTAR)/FQ4**2 
      DSS(4,4)=FQSTAR/FQ4+FS(4)*(DQSTARDS(4)*FQ4-FB1(4)*FQSTAR)/FQ4**2  
      DSS(4,5)=FS(4)*(DQSTARDS(5)*FQ4-FB1(5)*FQSTAR)/FQ4**2 
      DSS(4,6)=FS(4)*(DQSTARDS(6)*FQ4-FB1(6)*FQSTAR)/FQ4**2  

      DSS(5,1)=FS(5)*(DQSTARDS(1)*FQ4-FB1(1)*FQSTAR)/FQ4**2  
      DSS(5,2)=FS(5)*(DQSTARDS(2)*FQ4-FB1(2)*FQSTAR)/FQ4**2  
      DSS(5,3)=FS(5)*(DQSTARDS(3)*FQ4-FB1(3)*FQSTAR)/FQ4**2 
      DSS(5,4)=FS(5)*(DQSTARDS(4)*FQ4-FB1(4)*FQSTAR)/FQ4**2  
      DSS(5,5)=FQSTAR/FQ4+FS(5)*(DQSTARDS(5)*FQ4-FB1(5)*FQSTAR)/FQ4**2 
      DSS(5,6)=FS(5)*(DQSTARDS(6)*FQ4-FB1(6)*FQSTAR)/FQ4**2 

      DSS(6,1)=FS(6)*(DQSTARDS(1)*FQ4-FB1(1)*FQSTAR)/FQ4**2  
      DSS(6,2)=FS(6)*(DQSTARDS(2)*FQ4-FB1(2)*FQSTAR)/FQ4**2  
      DSS(6,3)=FS(6)*(DQSTARDS(3)*FQ4-FB1(3)*FQSTAR)/FQ4**2  
      DSS(6,4)=FS(6)*(DQSTARDS(4)*FQ4-FB1(4)*FQSTAR)/FQ4**2  
      DSS(6,5)=FS(6)*(DQSTARDS(5)*FQ4-FB1(5)*FQSTAR)/FQ4**2  
      DSS(6,6)=FQSTAR/FQ4+FS(6)*(DQSTARDS(6)*FQ4-FB1(6)*FQSTAR)/FQ4**2  
      DO I=1,6
          FD(I)=DFDP*FA(I)+DFDQ*FB(I)
      END DO

      DO I=1,6
	  FD1(I)=0
	  DO J=1,6
	  DO K=1,6
          FD1(I)=FD1(I)+FD(J)*DSS(J,K)*DS(K,I)
      END DO
      END DO	  
      END DO
	  
      DO I=1,6
		  FC(I)=DGDP*FA(I)+DGDQ*FB(I)
      END DO	  
     
      DO I=1,6
          FE(I)=0.0
          FE1(I)=0.0		  
          DO J=1,6
              FE(I)=FE(I)+FD(J)*FDE(I,J)
              FE1(I)=FE1(I)+FD1(J)*FDE(I,J)		  
          END DO
      END DO

      FF=0.0
      DO I=1,6
      FF=FF+FE1(I)*FC(I)
      END DO

      DO I=1,6
          FG(I)=0.0
          DO J=1,6
              FG(I)=FG(I)+FDE(I,J)*FC(J)
          END DO
      END DO
      DO I=1,6
          DO J=1,6
              FH(I,J)=FG(I)*FE(J)
          END DO
      END DO
	  FR1=0
	  FR2=0	  
      DO I=1,6
	  FR1=FR1+FE(I)*DER(I)
	  FR2=FR2+FE1(I)*DER(I)	  
      END DO	  
	  
	  FR=FR1/FR2
	  
      DO I=1,6
          DO J=1,6
              DEP(I,J)=FDE(I,J)-FH(I,J)/(FKP+FF)/FR 
          END DO
      END DO
	   RAMA=0
      DO I=1,6
          RAMA=RAMA+FE(I)*DER(I)/(FKP+FF)/FR
      END DO
          DEVP=RAMA*DGDP	  
	  
      FUNLOAD=0.0
      DO I=1,6
          FUNLOAD=FUNLOAD+FD(I)*DER(I)
      END DO
      IF(FUNLOAD.LE.0)THEN
          DO I=1,6
              DO J=1,6
                  DEP(I,J)=FDE(I,J)
              END DO
          END DO
          
		  DEVP=0
      END IF

      END
      
      
      
      
