FUNCTION COMPUTE_MULTISCALE_DISTS, Img1, Img2

dims = GET_DIMENSIONS(Img1)

;ImgDist = DBLARR(dims[0],dims[1],dims[2])
ImgDist = FLTARR(4,dims[0],dims[1],dims[2])


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
         ImgDist[0,k,i,j] = BHATTACHARYYA(Par1.Mu, Par2.Mu, Par1.Sigma, Par2.Sigma) 
         IF ~finite(ImgDist[0,k,i,j]) THEN ImgDist[0,k,i,j] = 0
         t1 = systime(/seconds)
         countTimeBath += t1-t0
      
         t0 = systime(/seconds)
         ;Jeffries-Matusita
         ;1: 
         ImgDist[1,k,i,j] = 2*(1 - EXP(-1*BHATTACHARYYA(Par1.Mu, Par2.Mu, Par1.Sigma, Par2.Sigma)))      
         IF ~finite(ImgDist[1,k,i,j]) THEN ImgDist[1,k,i,j] = 0
         t1 = systime(/seconds)
          countTimeJM += t1-t0
      
        t0 = systime(/seconds)
         ;Hellinger
         ;2:
         ImgDist[2,k,i,j] = HELLINGER(Par1.Mu, Par2.Mu, Par1.Sigma, Par2.Sigma)
         IF ~finite(ImgDist[2,k,i,j]) THEN ImgDist[2,k,i,j] = 0
         t1 = systime(/seconds)
         countTimeHell += t1-t0
         
         t0 = systime(/seconds)
         ;Kullback-Leilber
         ;3:
         ImgDist[3,k,i,j] = KULLBACK_LEIBLER(Par1.Mu, Par2.Mu, Par1.Sigma, Par2.Sigma)
         IF ~finite(ImgDist[3,k,i,j]) THEN ImgDist[3,k,i,j] = 0
         t1 = systime(/seconds)
          countTimeKL += t1-t0
    
      ;ENDCASE
    
    ENDFOR
  ENDFOR
ENDFOR
print, 'Distância Estocástica...'
print, 'time bath...',countTimeBath
print, 'time jm...',countTimeJM
print, 'time hell...',countTimeHell
print, 'time kl...',countTimeKL

Return, ImgDist
END