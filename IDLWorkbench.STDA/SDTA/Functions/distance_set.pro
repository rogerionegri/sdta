PRO DISTANCE_SET
   ;Este PRO é apenas para referencia ao
   ;conjunto de distancias disponibilizados abaixo
END


;Adicionado...
FUNCTION BHATTACHARYYA, MuX, MuY, SigmaX, SigmaY
  ;###########################################

  ;Distancia de Bhattacharyya, simplificada para distribuição gaussiana multivariada
  P = (SigmaX + SigmaY)/2.0
  D = (1/8.0)*TRANSPOSE(MuX - MuY)#INVERT(P)#(MuX - MuY) + $
    0.5*ALOG(DETERM(P, /check, /double)/SQRT(DETERM(SigmaX, /check, /double)*DETERM(SigmaY, /check, /double)))

  Return, D
END




;Adicionado...
FUNCTION HELLINGER, MuX, MuY, SigmaX, SigmaY
  ;###########################################

  ;Distancia de Hellinger, simplificada para distribuição gaussiana multivariada

  P = 0.5*(SigmaX + SigmaY)
  A = DETERM(SigmaX, /check, /double)^0.25 * DETERM(SigmaY, /check, /double)^0.25
  B = DETERM( P , /check, /double)^0.5
  D = -0.125 * TRANSPOSE(MuX - MuY)#INVERT(P)#(MuX - MuY)

  ;H = 1.0 - ( (A/B) * EXP(D) )
  H = SQRT(1.0 - ( (A/B) * EXP(D) ))

  Return, H
END




;Adicionado...
FUNCTION KULLBACK_LEIBLER, Mu0, Mu1, Sigma0, Sigma1
  ;###########################################

  ;Distancia de Kullback-Leibler para distribuição gaussiana multivariada

  invSigma0 = INVERT(Sigma0)
  invSigma1 = INVERT(Sigma1)
  k = N_ELEMENTS(Mu0)
  detSigma0 = DETERM(Sigma0, /check, /double)
  detSigma1 = DETERM(Sigma1, /check, /double)

  D01 = 0.5 * ( TRACE(invSigma1 # Sigma0) + TRANSPOSE(Mu1 - Mu0) # invSigma1 # (Mu1 - Mu0) - k + ALOG(detSigma1/detSigma0) )
  D10 = 0.5 * ( TRACE(invSigma0 # Sigma1) + TRANSPOSE(Mu0 - Mu1) # invSigma0 # (Mu0 - Mu1) - k + ALOG(detSigma0/detSigma1) )

  ;cabrito tosco....
  if ~finite(d01) then d01 = randomu(seed,1)
  if ~finite(d10) then d10 = randomu(seed,1)

  D = 0.5*(D01 + D10)

  Return, D
END



;Adicionado...
FUNCTION BOCA_NORMA, Mu0, Mu1, Sigma0, Sigma1
  ;###########################################

  ;experimentacao

  invSigma0 = INVERT(Sigma0)
  invSigma1 = INVERT(Sigma1)
  k = N_ELEMENTS(Mu0)

  AA =  TOTAL(ABS(TRANSPOSE(Mu1 - Mu0) # invSigma1 # (Mu1 - Mu0)))

  BB = invSigma0#Sigma1
  CC = invSigma1#Sigma0
  DD = TOTAL(ABS( BB + CC - 2*IDENTITY(k) ))

  Return, (AA + DD)
END



;###########################################
FUNCTION LP2_METRIC, Mu0, Mu1, Sigma0, Sigma1

  sumD = 0D
  FOR i = 0, N_ELEMENTS(Mu0)-1 DO BEGIN
    A = (SQRT(!PI)*(SQRT(Sigma0[i,i]) + SQRT(Sigma1[i,i]))) / (2 * !PI * SQRT(Sigma0[i,i]) * SQRT(Sigma1[i,i]) )
    B =  SQRT(2*!PI)/(!PI * SQRT(Sigma0[i,i] + Sigma1[i,i]))
    C = EXP(- ( (Mu0[i] - Mu1[i])^2 ) / (2*( Sigma0[i,i] + Sigma1[i,i] )) )

    SumD += SQRT( A - B*C )
  ENDFOR

  Return, SumD
END


;###########################################
FUNCTION MAHALANOBIS_ADAPT, Mu0, Mu1, Sigma0, Sigma1


  Dist = TRANSPOSE(Mu0 - Mu1) # Sigma1 # (Mu0 - Mu1)


  Return, SQRT(Dist)
END