FUNCTION COMPUTE_MULTISCALE_HT, Img1, Img2, Atts, winSize

M = N_ELEMENTS(Atts) + (N_ELEMENTS(Atts) * (N_ELEMENTS(Atts)+1))/2
nei = winSize*winSize

dims = GET_DIMENSIONS(Img1)

;ImgDist = DBLARR(dims[0],dims[1],dims[2])
ImgStat = FLTARR(4,dims[0],dims[1],dims[2])

countTimeBath = 0.0
countTimeJM = 0.0
countTimeKL = 0.0
countTimeHell = 0.0
FOR k = 0, dims[0]-1 DO BEGIN
  FOR i = 0, dims[1]-1 DO BEGIN
    FOR j = 0, dims[2]-1 DO BEGIN
      Par1 = *Img1[k,i,j]
      Par2 = *Img2[k,i,j]
      
      ;CASE distType OF
      
      t0 = systime(/seconds)
         ;Bathacharryya
         ;0: 
         ImgStat[0,k,i,j] = 1.0 - CHISQR_PDF(BHATTACHARYYA(Par1.Mu, Par2.Mu, Par1.Sigma, Par2.Sigma) * 8.0 * (nei/2.0) , M)
         IF ~finite(ImgStat[0,k,i,j]) THEN ImgStat[0,k,i,j] = 0
         t1 = systime(/seconds)
         countTimeBath += t1-t0
             
      t0 = systime(/seconds)
         ;Jeffries-Matusita
         ;1: 
         ImgStat[1,k,i,j] = 1.0 - CHISQR_PDF(2*(1 - EXP(-1*BHATTACHARYYA(Par1.Mu, Par2.Mu, Par1.Sigma, Par2.Sigma))) * 8.0 * (nei/2.0) ,M)    
         IF ~finite(ImgStat[1,k,i,j]) THEN ImgStat[1,k,i,j] = 0
         t1 = systime(/seconds)
         countTimeJM += t1-t0
      
      t0 = systime(/seconds)
         ;Hellinger
         ;2:
         ImgStat[2,k,i,j] = 1.0 - CHISQR_PDF(HELLINGER(Par1.Mu, Par2.Mu, Par1.Sigma, Par2.Sigma) * 8.0 * (nei/2.0), M)
         IF ~finite(ImgStat[2,k,i,j]) THEN ImgStat[2,k,i,j] = 0
         t1 = systime(/seconds)
         countTimeHell += t1-t0
         
         
      t0 = systime(/seconds)
         ;Kullback-Leilber
         ;3:
         ImgStat[3,k,i,j] = 1.0 -  CHISQR_PDF(KULLBACK_LEIBLER(Par1.Mu, Par2.Mu, Par1.Sigma, Par2.Sigma) * 2.0 * alog(10) * (nei/2.0), M)
         IF ~finite(ImgStat[3,k,i,j]) THEN ImgStat[3,k,i,j] = 0
         t1 = systime(/seconds)
         countTimeKL += t1-t0
    
    
      ;ENDCASE
    
    ENDFOR
  ENDFOR
ENDFOR

print, 'Hypothesis test...'
print, 'time bath...',countTimeBath
print, 'time jm...',countTimeJM
print, 'time hell...',countTimeHell
print, 'time kl...',countTimeKL

Return, ImgStat
END