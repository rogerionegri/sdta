PRO ASSESS_CHANDET_REPORT, PATH_REPORT, img, PATH_ROI, pars, time 

   PtrROI = ASCII_READ_ROI(PATH_ROI)

   roiChange = *PtrROI[0]
   roiNonChange = *PtrROI[1]

   ;calculo da acuracia
   posTP = WHERE(img[roiChange.Roilex] EQ 1)
   posFN = WHERE(img[roiChange.Roilex] NE 1)
   TP = N_ELEMENTS(posTP)
   FN = N_ELEMENTS(posFN)
   
   posTN = WHERE(img[roiNonChange.Roilex] EQ 0)
   posFP = WHERE(img[roiNonChange.Roilex] NE 0)
   TN = N_ELEMENTS(posTN)
   FP = N_ELEMENTS(posFP)

   Ac = (FLOAT(TP) + TN)/N_ELEMENTS([roiChange.Roilex , roiNonChange.Roilex])
   Pr = FLOAT(TP)/(TP + FP)
   Re = FLOAT(TP)/(TP + FN)
   F1 = (2*Pr * Re)/(Pr + Re)
   MCC = ( FLOAT(TP)*TN - FLOAT(FP)*FN )/SQRT( (FLOAT(TP)+FP)*(FLOAT(FN)+TN)*(FLOAT(FP)+TN)*(FLOAT(TP)+FN))
  
   infoAssess = ''
  
   FOR i = 0, N_ELEMENTS(pars)-1 DO infoAssess += STRTRIM(pars[i],1) +' ; '   
  
   infoAssess += STRTRIM(STRING(AC),1) +' ; '+ STRTRIM(STRING(Pr),1) +' ; '+ $
                STRTRIM(STRING(Re),1) +' ; '+ STRTRIM(STRING(F1),1) +' ; '+ STRTRIM(STRING(MCC),1) +' ; '+ STRTRIM(STRING(time),1) 

  
   OPENW, Arq, PATH_REPORT, /APPEND, /GET_LUN
   PRINTF, Arq, infoAssess
   CLOSE, Arq
   FREE_LUN, Arq 
END



;################################################################
PRO PUT_HEADER, PATH_REPORT, parNames
   infoAssess = ''
   FOR i = 0, N_ELEMENTS(parNames)-1 DO infoAssess += STRTRIM(parNames[i],1) +' ; '

   infoAssess += 'Accuracy ; Precision ; Recall ; F1-Score ; MCC ; time(sec.)'

   OPENW, Arq, PATH_REPORT, /APPEND, /GET_LUN
   PRINTF, Arq, infoAssess
   FREE_LUN, Arq
END
