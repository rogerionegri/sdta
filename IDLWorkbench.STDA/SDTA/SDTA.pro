@BocaLib.pro

@compute_parameters.pro
@Functions_ChangeDetection.pro

@kitter_illingworth.pro
@stoch_distances_set.pro

@assess_chandet_report.pro
@compute_multiscale_ht.pro

PRO PROJETO_GODOY_V0

;-------------------------------------------------------------
  PATH_IMG1 = 'img1.tif'                    ;Insert image 1
  PATH_IMG2 = 'img2.tif'                    ;Insert image 2
  PATH_OUTPUT = 'ResultDistancePath'        ;Insert path to save distance map
  PATH_OUTPUT_THRESHOLD = 'ResultMapPath'   ;Insert path to save change map
  PATH_REPORT = 'ReportPath'                ;Insert path to save evaluation of results
  PATH_ROI = 'groundTruth.txt'              ;Inserir amostras de referência
  PATH_SAVE_VAR = 'Save.sav'                ;Insert variables backup path
  Atts = [0,1]                              ;image bands
  sizes = 1                                 ;neighborhood radii
  ;distType = 0   ;Stochastic Distance: 0 - Bhattacharyya, 1 - Jeffries-Matusita, 2 - Hellinger, 3 - Kullback-Leibler
;-------------------------------------------------------------

  winSizes = (2*indgen(sizes)+3)   ;<<<neighborhood windows size (2i+1)x(2i+1); i=1,2...,7

  Img1 = OPEN_IMAGE(PATH_IMG1, Atts)
  Img2 = OPEN_IMAGE(PATH_IMG2, Atts)

  t0 = systime(/seconds)
  ParRegs1 = GET_NEIGH_PARS_MULTI(Img1,winSizes)
  t1 = systime(/seconds)
  print, 'time...', t1-t0
  ParRegs2 = GET_NEIGH_PARS_MULTI(Img2,winSizes)
  t2 = systime(/seconds)
  print, 'time...', t2-t1
  
 ;stop;
  
  ;this function only calculates the distances between the two images, for each of the neighborhood sizes
  SAVE, Filename=PATH_SAVE_VAR
  ImgDists = COMPUTE_MULTISCALE_DISTS(ParRegs1,ParRegs2)

  PTR_FREE, ParRegs1, ParRegs2
  HEAP_GC

  ;stop
  ;here it just saves the result of the previous step
  temp = FLTARR(N_ELEMENTS(Imgdists[0,*,0,0]), N_ELEMENTS(Imgdists[0,0,*,0]),N_ELEMENTS(Imgdists[0,0,0,*]))
  FOR i = 0, 3 DO BEGIN
     temp[*,*,*] = ImgDists[i,*,*,*]
     savePath = PATH_OUTPUT + '/DE/DistImage__distType_'+STRTRIM(STRING(i),1)+'__'+STRTRIM(STRING(winSizes),1)+'.tif'
     WRITE_TIFF, savePath, temp, /DOUBLE
  ENDFOR
  
  ;stop
  temp = FLTARR(N_ELEMENTS(Imgdists[0,0,*,0]),N_ELEMENTS(Imgdists[0,0,0,*]))
  FOR op = 0, 1 DO BEGIN
     FOR di = 0, 3 DO BEGIN
  
        FOR i = 0, N_ELEMENTS(Imgdists[0,*,0,0])-1 DO BEGIN
           temp[*,*] = Imgdists[di,i,*,*]
           
           if di eq 3 then temp[*,*] = alog(temp[*,*]) - min(alog(temp))
           
           t0 = systime(/seconds)
           kiwRes = KIW_THRESHOLD(temp,op)
           t1 = systime(/seconds)
           
           savePath = PATH_OUTPUT_THRESHOLD + 'Thres__distType_'+STRTRIM(STRING(di),1) +'__'+STRTRIM(STRING(winSizes),1)+'__opt_'+STRTRIM(STRING(op),1) + '__win_'+STRTRIM(STRING(winSizes[i]),1) +'.tif'
           WRITE_TIFF, savePath, kiwRes    
           
           ASSESS_CHANDET_REPORT, PATH_REPORT + '.txt', kiwRes, PATH_ROI, [di,op,winSizes[i]], [t1-t0]
        ENDFOR
  
     ENDFOR
  ENDFOR  
  
  t3 = systime(/seconds)
  print, 'Total Thresholding Time...', t3-t2 
  
stop  
  
END
