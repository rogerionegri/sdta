PRO Functions_ChangeDetection
   ;caller
END


;#####################################
FUNCTION GET_QUICK_NEIGH_VALS,pi,pj,WinX,WinY,NB,NC,NL,Image

Neighs = MAKE_ARRAY(NB,NC*NL,TYPE = SIZE(Image,/TYPE))
count = 0
FOR i = LONG(-WinX)/2, LONG(+WinX)/2 DO BEGIN
   FOR j = LONG(-WinY)/2, LONG(+WinY)/2 DO BEGIN
      IF ((((pi + i) GE 0) AND ((pi + i) LT NC)) AND (((pj + j) GE 0) AND ((pj + j) LT NL)) ) THEN BEGIN
          Neighs[*,count] = Image[*,pi+i,pj+j]
          count++
      ENDIF    
   ENDFOR
ENDFOR

Return, Neighs[*,0:count-1]
END