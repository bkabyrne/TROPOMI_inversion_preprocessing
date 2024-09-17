;;+
;####################### MOPITT-v8 Data Processor IDL Code ###############################
;;
;;                 Multi-Mission Observation Operator Project
;;                          
;;                             (M2O2)
;;
;;                       L2-Sharp (L2#) data generation
;;
;;                         Copyright 2016,
;;              California Institute of Technology.
;;                       ALL RIGHTS RESERVED.
;;           U.S. Government Sponsorship acknowledged.
;;
;;
;; ROUTINE NAME
;;    MOPITT_xCO
;;
;; FILE NAME 
;;    MOPITT_v8.pro
;;
;; Purpose
;;    Generates L2 sharp data for column CO for data assimilation applications 
;;     
;; AUTHOR
;;    Meemong Lee, July/01/2020
;     Updated by Brendan Byrne, 2023 (there were errors in treatment
;                                     of log retrievals)
;
;; DESCRIPTION
;#################### MOPITT-v8 Data Processor IDL Code ############################

 Pro MOPITT_xCO, year, sM, eM
    
    ; cuf off pressure levels above 68.1295 hPa
    LSave = 0
    dataPath = '/npbackupp19/bbyrne1/MOPITT_L2v19/'
    rootPath = '/nobackup/bbyrne1/MOPITT_v9_xCO_Log_20221216/'
    
    Record = create_struct( $
      'nLon', 144, $
      'nLat', 91, $
      'dLon', 2.5, $
      'dLat', 2.0, $
      'xCO',fltarr(144, 91), $
      'xCO_molC',fltarr(144, 91), $
      'err', fltarr(144, 91), $
      'UnitErr', fltarr(144, 91), $
      'count',fltarr(144, 91) $
    )
   
   Daily = Record
   Monthly = Record 
    
   yearS = string(year, format='(i4,"/")')
    
    for month = sM, eM do begin
      print, 'month = ', month

      Monthly.count(*,*) = 0
      Monthly.xCO(*,*) = 0
      Monthly.xCO_molC(*,*) = 0
      Monthly.err(*,*) = 0
      Monthly.UnitErr(*,*) = 0

      monthS = string(month, format='(i2.2,"/")')
      outPath = rootPath + yearS + monthS
      FILE_MKDIR, outPath

      for day = 1, 31 do begin  
        Daily.count(*,*) = 0
        Daily.xCO(*,*) = 0
        Daily.err(*,*) = 0
        Daily.xCO_molC(*,*) = 0
        Daily.unitErr(*,*) = 0
       
        dateS = string(year,month,day,format='(i4,i2.2,i2.2)') 
                                ;filename = dataPath +
                                ;yearS+monthS+"MOP02N-" + dateS
                                ;+"-L2V19.9.2.he5"
        filename = dataPath + "MOP02J-" + dateS +"-L2V19.9.3.he5"

        output1 = outPath + string(day, format='(i2.2, ".nc")')
        output2 = outPath + string(day, format='("map-", i2.2, ".nc")')
        if (FILE_TEST(filename) eq 1) then begin
          Read_MOPITT_xCO, day, filename, output1, output2, Daily, Lsave
       
          for iy = 0, 90 do begin
          for ix = 0, 143 do begin
            if (Daily.count(ix,iy) gt 0) then begin 
              Monthly.xCO(ix,iy) = Monthly.xCO(ix,iy) + Daily.xCO(ix,iy)
              Monthly.xCO_molC(ix,iy) = Monthly.xCO_molC(ix,iy) + Daily.xCO_molC(ix,iy)
              Monthly.err(ix,iy) = Monthly.err(ix,iy) + Daily.err(ix,iy)
              Monthly.UNitErr(ix,iy) = Monthly.UnitErr(ix,iy) + Daily.UnitErr(ix,iy)
              Monthly.count(ix,iy) = Monthly.count(ix,iy) + Daily.count(ix,iy)
            endif
          endfor
          endfor
        endif else begin
          print, filename, ' does not exist.'
        endelse
      endfor
      
      output = rootPath + yearS + string(month,FORMAT='("summary-",i2.2,".nc")')
      SaveMap,output, Monthly
    endfor 
    end
    
    Pro Read_MOPITT_xCO, day, filename, output, output2, Daily, LSave
    
    LCheck = 0 
    hourlySamples = fltarr(24)
    hourlySamples(*) = 0
        
    ;fill value
    fillV = -9999.0
    ;The retrieved items include estimated value and its uncertainty
    CH1 = 0
    CH2 = 1

    fileId = H5F_OPEN(filename)
    group = '/HDFEOS/SWATHS/MOP02/Geolocation Fields'
    gid = H5G_OPEN(fileId, group)   
    dsInfo = H5G_GET_OBJINFO(gid,'Latitude')    
    
    dsId = H5D_OPEN(gid,'Latitude')
    latitude = H5D_Read(dsId)
    H5D_CLOSE, dsId
    
    dsId = H5D_OPEN(gid, 'Longitude')
    longitude = H5D_Read(dsId)
    H5D_CLOSE, dsId
    
    dsId = H5D_OPEN(gid,  'SecondsinDay')
    hour = H5D_Read(dsId)
    H5D_CLOSE, dsId
    
    dsId = H5D_OPEN(gid, 'Pressure')
    pressureX = H5D_Read(dsId)
    H5D_CLOSE, dsId
    H5G_CLOSE, gid

    nLevels = n_elements(pressureX(*))+1
    nTargets = n_elements(hour(*))
     
    ; allocate memory for the variables    
    pressure = fltarr(nLevels, nTargets)
    Xa = fltarr(nLevels, nTargets)
    colXa = fltarr(nTargets)
    colX = fltarr(nTargets)
    colAK = fltarr(nLevels, nTargets)  
    colErr = fltarr(nTargets)
    weight = fltarr(nLevels, nTargets)
    rLevels = fltarr(nTargets)

    group = '/HDFEOS/SWATHS/MOP02/Data Fields/'
    gid = H5G_OPEN(fileId, group)  
    
    ; data quality related  
    dsId = H5D_Open(gid, 'SurfaceIndex')
    surfaceIndex = H5D_Read(dsId);
    H5D_CLOSE, dsId
    
    dsId = H5D_Open(gid, 'CloudDescription')
    cloudDesc = H5D_Read(dsId);
    H5D_CLOSE, dsId

    dsId = H5D_OPEN(gid,'RetrievalAnomalyDiagnostic')
    anomaly = H5D_Read(dsId)
    H5D_CLOSE, dsId

    dsId = H5D_Open(gid, 'SurfacePressure')
    pSurf = H5D_Read(dsId);
    H5D_CLOSE, dsId

    dsId = H5D_Open(gid, 'APrioriCOMixingRatioProfile')
    XaTemp = H5D_Read(dsId);	
    H5D_CLOSE, dsId

    dsId = H5D_Open(gid, 'APrioriCOSurfaceMixingRatio')
    xaSurf = H5D_Read(dsId);
    H5D_CLOSE, dsId
    
    dsId = H5D_Open(gid, 'RetrievedCOTotalColumn')
    totalCol = H5D_Read(dsId);
    H5D_CLOSE, dsId

    dsId = H5D_Open(gid, 'APrioriCOTotalColumn')
    TotalColXa = H5D_Read(dsId);
    H5D_CLOSE, dsId

    dsId = H5D_Open(gid, 'WaterVaporColumn')
    waterVapor = H5D_Read(dsId);
    H5D_CLOSE, dsId

    dsId = H5D_Open(gid, 'TotalColumnAveragingKernelDimless')
    totalColAK = H5D_Read(dsId);
    H5D_CLOSE, dsId

    dsId = H5D_Open(gid, 'RetrievedCOTotalColumnDiagnostics')
    totalColErr = H5D_Read(dsId);
    H5D_CLOSE, dsId

    dsId = H5D_Open(gid, 'DryAirColumn')
    dryAir = H5D_Read(dsId);
    H5D_CLOSE, dsId

    H5G_CLOSE,gid
    H5F_CLOSE, fileid 
    
    nX = 0L
   
    profile = create_struct ( $
      'PS', fltarr(nLevels), $
      'Xa', fltarr(nLevels) $
       )

    column = create_struct( $
      'PW', fltarr(nLevels), $
      'X', 0.0, $
      'Xa', 0.0, $
      'Err', 0.0, $
      'AK', fltarr(nLevels) $
      )
    
    quality = create_struct( $
      'totalColumn', 0L, $
      'cloud', 0L, $
      'error', 0L $
    )

  negAKCount = 0L
  surfNegAKCount = 0L
    
  for k1 = 0, nTargets-1 do begin
      usability =1
      
      if (totalCol(CH1,k1) lt 5.0e17) then begin
        quality.totalColumn = quality.totalColumn + 1
        ;usability = 0
      endif
      ; cloundDesc
      ; 1 - MOPCLD only clear. thermal only
      ; 2 - MOPCLD and MODIS cloud mask agree on clear (chosen)
      ; 3 = MODIS cloud mask only clear
      ; 4 - MOPCLD overriding MODIS cloud mas over low clouds
      ; 5 - MODIS cloud mask only, clear over polar regions
      
      if (cloudDesc(k1) ne 2.0 or anomaly(k1) ne 0) then begin
        quality.cloud = quality.cloud + 1
        usability = 0
      endif 

      if (usability eq 1) then begin       
        column.AK(*) = 0.0
        column.PW(*) = 0.0
        profile.PS(*) = 0.0
        profile.Xa(*) = 0.0
        
       ; find the valid retrieval levels above the surface pressure
        L0 = 0
        while (L0 lt nLevels-2 and psurf(k1) lt pressureX(L0))do L0 = L0+1    
        rLev = nLevels-L0 
        sLev = L0
        profile.PS(0) = psurf(k1)
        profile.PS(1:rLev-1) = pressureX(L0:nLevels-2)
        profile.Xa(0) = xaSurf(CH1, K1)
        profile.Xa(1:rLev-1) = XaTemp(CH1,L0:nLevels-2, K1)	; in ppb
        column.PW = PressureWeight(profile.PS, rLev)        
        column.AK(0:rLev-1) = totalColAK(L0:nLevels-1, K1)

        ; compute column xCO
        column.Err = TotalColErr(CH1,K1)/TotalCol(CH1,k1)

        ; AK starts from the surface, check for negative value
        negAK = 0
        for k = 0, rLev-1  do begin
          if (column.AK(k) lt 0) then negAK = negAK + 1
          column.AK(k) = column.PW(k) * column.AK(k)
        endfor

        if (negAK gt 0) then begin
	  ; allow AK to be negative at the surface or at the top
          if (negAK eq 1 and (column.AK(0) lt 0 or column.AK(rLev-1) lt 0)) then begin
            negAK = 0
            surfNegAKCount = surfNegAKCount + 1
          endif else begin
            negAKCount = negAKCount + 1
          endelse

        endif
         
        if (column.Err lt 0.3 and negAK eq 0) then begin
          longitude[nX] = longitude[k1]
          latitude[nX] = latitude[k1]
          hour[nX]= hour[k1]/3600.0
          hX = Floor(hour(nX))
          hourlySamples(hX) = hourlySamples(hX)+1
        
          pressure[sLev:nLevels-1,nX] = profile.PS[0:rLev-1]
          colAK[slev:nLevels-1,nX] = column.AK(0:rlev-1)

          column.X = totalCol(CH1,k1)/dryAir(k1) * 1.0e9
          tmpXa = totalColXa(k1)/dryAir(K1) * 1.0e9
          column.Xa = TOTAL(column.PW*profile.Xa)

          weight[*,nX] = column.PW
          ; convert ppb to v/v
          colX[nX] = column.X * 1.0e-9
          colXa[nX] = column.Xa * 1.0e-9
          colErr[nX] = column.Err * colX[nX]
          Xa(slev:nLevels-1,nX) = profile.Xa[0:rLev-1] * 1.0e-9

          ;fill the invalid cells
          if (sLev gt 0) then begin
            colAK[0:sLev-1, nX] = -999
            pressure[0:sLev-1, nX] = -999
            Xa[0:sLev-1,nX] = -999
            ;print, colAK(*,nX)
          endif
          
          ix = long(longitude[nX] + 180) / daily.dLon
          iy = long(latitude[nX] + 90) / daily.dLat
          if (ix ge daily.nLon) then ix = ix-daily.nLon 
             
          Daily.xCO[ix,iy] = Daily.xCO[ix,iy] + column.X
          Daily.xCO_molC[ix,iy] = Daily.xCO_molC[ix,iy] + totalCol(CH1,K1)
          Daily.err[ix,iy] = Daily.err[ix,iy] + column.Err
          Daily.UnitErr[ix,iy] = Daily.UnitErr[ix,iy] + (tmpXa- column.Xa)/column.Xa
          Daily.count[ix,iy] = Daily.count[ix,iy] + 1
        
        ; increment the good sample count
          nX= nX+1
          endif else begin
            quality.error = quality.error+1
          endelse
        endif ;usability loop
    endfor  ;sample loop

    print, day, '  nSamples = ', nTargets, '  nX =', nX, ' negAK = ', negAKCount, surfNegAKCount

    if (nX lt 1000) then begin
      print, '< 1000 good samples found'
      return
    endif
   
    if (LSave) then SaveMap, output2, Daily.xCO, Daily.err, Daily.count


    ; Create a daily NetCDF dile
     id = NCDF_CREATE(output,/CLOBBER)
     NCDF_ATTPUT,id,'TITLE','MOPITT-xCO-L2#',/GLOBAL
     NCDF_ATTPUT,id,'Source','MOPITT-L2_V8',/GLOBAL
     NCDF_ATTPUT,id,'Contact', 'Brendan Byrne <brendan.k.byrne@jpl.nasa.gov>',/GLOBAL
     NCDF_ATTPUT,id,'Created on',  sysTime(), /GLOBAL
     NCDF_ATTPUT,id,'AK Type', 'VARIABLE', /GLOBAL
     NCDF_ATTPUT,id,'AK Space', 'LogVMR', /Global
     NCDF_ATTPUT,id,'Pressure Level Type', 'VARIABLE', /GLOBAL
     NCDF_ATTPUT,id,'Uncertainty Type', 'RMS', /GLOBAL
     NCDF_ATTPUT,id,'_FillValue', -999.0, /GLOBAL
  
     dim0 = NCDF_DIMDEF(id, 'nSamples', nX)
     dim1 = NCDF_DIMDEF(id, 'maxLevels', nLevels) 
     dim2 = NCDF_DIMDEF(id, 'nHours', 24)
            
     vid1 = NCDF_VARDEF(id, 'latitude', dim0, /FLOAT)
            NCDF_ATTPUT,id, vid1, 'long_name', 'latitude'
            NCDF_ATTPUT,id, vid1, 'unit', 'degrees_north'
            
     vid2 = NCDF_VARDEF(id, 'longitude', dim0, /FLOAT)
            NCDF_ATTPUT,id, vid2, 'long_name', 'longitude'
            NCDF_ATTPUT,id, vid2, 'unit', 'degrees_east'
            
     vid3 = NCDF_VARDEF(id, 'samplesPerHour', dim2, /FLOAT)
            NCDF_ATTPUT,id,vid3,'long_name','hourly sample count'
            NCDF_ATTPUT,id,vid3,'unit','none'

     vid4 = NCDF_VARDEF(id, 'time', dim0, /FLOAT)
            NCDF_ATTPUT,id,vid4,'long_name','UTC time of day at location'
            NCDF_ATTPUT,id,vid4,'unit','hrs'
            
     vid5 = NCDF_VARDEF(id, 'pressure', [dim1,dim0], /FLOAT)
            NCDF_ATTPUT, id, vid5, 'long_name', 'sounding pressure profile (from surface pressure)'
            NCDF_ATTPUT,id,vid5,'unit','hPa'
            NCDF_ATTPUT, id, vid5, '_FillValue', -999.0

     vid6 = NCDF_VARDEF(id, 'xCO', dim0, /FLOAT)
            NCDF_ATTPUT, id, vid6,'long_name', 'column CO (AK(CO-COprior) + xCO-prior)'
            NCDF_ATTPUT, id, vid6,'unit','log(v/v)'
            NCDF_ATTPUT, id, vid6,'_FillValue', -999.0
            
     vid7 = NCDF_VARDEF(id, 'xCO-apriori', dim0, /FLOAT)
            NCDF_ATTPUT, id, vid7, 'long_name', 'apriori column CO (pressure weigted sum of CO apriori)'
            NCDF_ATTPUT, id, vid7,'unit','log(v/v)'
            NCDF_ATTPUT, id, vid7, '_FillValue', -999.0
                   
     vid8 = NCDF_VARDEF(id, 'xCO-averagingKernel', [dim1, dim0], /FLOAT)
            NCDF_ATTPUT,id, vid8,'long_name', 'averaging kernel (pressure weighted averaging kernel)'
            NCDF_ATTPUT, id, vid8, 'unit', 'none'
            NCDF_ATTPUT, id, vid8, '_FillValue', -999.0
            
     vid9 = NCDF_VARDEF(id, 'xCO-uncertainty', dim0, /FLOAT)
            NCDF_ATTPUT, id, vid9, 'long_name', 'retrieval error (pressure weighted sum of the retrieval uncertainty)'
            NCDF_ATTPUT, id, vid9, 'unit', 'log(v/v)'
            NCDF_ATTPUT, id, vid9, '_FillValue', -999.0   
                 
     vid10 = NCDF_VARDEF(id, 'xCO-pressureWeight', [dim1,dim0], /FLOAT)
            NCDF_ATTPUT, id, vid10, 'long_name', 'pressure weight vector'
            NCDF_ATTPUT, id, vid10,'unit','none'
            NCDF_ATTPUT, id, vid10, '_FillValue', -999.0

     vid11 = NCDF_VARDEF(id, 'CO-apriori', [dim1,dim0], /FLOAT)
            NCDF_ATTPUT, id, vid11, 'long_name', 'apriori CO profile'
            NCDF_ATTPUT,id,vid11,'unit','v/v'
            NCDF_ATTPUT, id, vid11, '_FillValue', -999.0
            

     NCDF_CONTROL, id, /ENDEF
     
     NCDF_VARPUT, id, vid1, latitude(0:nX-1)
     NCDF_VARPUT, id, vid2, longitude(0:nX-1)
     NCDF_VARPUT, id, vid3, hourlySamples(*)            
     NCDF_VARPUT, id, vid4, hour(0:nX-1)    
     NCDF_VARPUT, id, vid5, pressure(*,0:nX-1)
     NCDF_VARPUT, id, vid11, ALOG(Xa(*, 0:nX-1))
             
     NCDF_VARPUT, id, vid6, ALOG(colX(0:nX-1))
     NCDF_VARPUT, id, vid7, ALOG(colXa(0:nX-1))
     NCDF_VARPUT, id, vid8, colAK(*,0:nX-1)
     NCDF_VARPUT, id, vid9, ALOG(colErr(0:nX-1)) ; THIS IS WRONG!!!!! SHOULD BE ALOG( 1 + COLERR/COLX)
     NCDF_VARPUT, id, vid10, weight(*,0:nX-1)
     NCDF_CLOSE, id      
     return   
    end
    
    Pro SaveMap, outPath, map

      nLons = n_elements(map.xCO(*, 0))
      nLats = n_elements(map.xCO(0,*))
        
      lonGrid = fltarr(nLons)
      latGrid = fltarr(nLats)
      dLon = 360.0/nLons
      dLat = 180.0/nLats
      
      for i=0, nLons-1 do begin
        lonGrid(i) = dLon*i - 180
      endfor
      for i=0, nLats-1 do begin
        latGrid(i) = dLat*i - 90
      endfor
    
      count = fltarr(nLons, nLats)
      meanxCO = fltarr(nLons, nLats)
      meanxCO_molC = fltarr(nLons, nLats)
      meanErr = fltarr(nLons, nLats)
      meanUnitErr = fltarr(nLons, nLats)
      
      for iy = 0, nLats-1 do begin
      for ix = 0, nLons-1 do begin
        tmp = map.count[ix,iy]
        if (tmp gt 0) then begin
          count[ix,iy] = tmp
          meanxCO[ix,iy] = map.xCO[ix,iy]/tmp
          meanxCO_molC[ix,iy] = map.xCO_molC[ix,iy]/tmp
          meanUnitErr[ix,iy] = map.unitErr[ix,iy]/tmp * 100.0
          meanErr[ix,iy] = map.err[ix,iy]/tmp * 100.0
        endif else begin
          count[ix,iy] = -999.0
          meanxCO[ix,iy] = -999.0
          meanxCO_molC[ix,iy] = -999.0
          meanErr[ix,iy] = -999.0
          meanUnitErr[ix,iy] = -999.0
        endelse
      endfor
      endfor
    
      id = NCDF_CREATE(outPath,/CLOBBER)
      NCDF_ATTPUT,id,'TITLE','mean MOPITT-v7-xCO',/GLOBAL
      NCDF_ATTPUT,id,'Author','Meemong Lee/JPL',/GLOBAL
      NCDF_ATTPUT,id,'Date',systime(),/GLOBAL
      NCDF_ATTPUT,id,'_FillValue', -999.0, /GLOBAL
        
      dim0 = NCDF_DIMDEF(id, 'lon', nLons)
      dim1 = NCDF_DIMDEF(id, 'lat', nLats)    

      vid0 = NCDF_VARDEF(id, 'lon', dim0, /FLOAT)       
      vid1 = NCDF_VARDEF(id, 'lat', dim1, /FLOAT)    
      
      vid3 = NCDF_VARDEF(id, 'xCO', [dim0,dim1],/FLOAT)
      NCDF_ATTPUT, id, vid3, '_FillValue', -999.0
      NCDF_ATTPUT, id,vid3,'unit','ppbv'
      
      vid4 = NCDF_VARDEF(id, 'sampleCount', [dim0,dim1],/FLOAT)
      NCDF_ATTPUT, id, vid4, '_FillValue', -999.0
      NCDF_ATTPUT, id,vid4,'unit','none'    
      
      vid5 = NCDF_VARDEF(id, 'err', [dim0,dim1],/FLOAT)
      NCDF_ATTPUT, id, vid5, '_FillValue', -999.0
      NCDF_ATTPUT, id,vid5,'unit','%'

      vid6 = NCDF_VARDEF(id, 'xCO_molC', [dim0,dim1],/FLOAT)
      NCDF_ATTPUT, id, vid6, '_FillValue', -999.0
      NCDF_ATTPUT, id,vid6,'unit','mole count / cm^2'

      vid7 = NCDF_VARDEF(id, 'unitErr', [dim0,dim1],/FLOAT)
      NCDF_ATTPUT, id, vid7, '_FillValue', -999.0
      NCDF_ATTPUT, id,vid7,'unit','%'
      
      NCDF_CONTROL, id, /ENDEF
      
      NCDF_VARPUT, id, vid0, lonGrid(*)
      NCDF_VARPUT, id, vid1, latGrid(*)
      NCDF_VARPUT, id, vid3, meanxCO(*,*)
      NCDF_VARPUT, id, vid4, count(*,*)
      NCDF_VARPUT, id, vid5, meanErr(*,*)
      NCDF_VARPUT, id, vid6, meanxCO_molC(*,*)
      NCDF_VARPUT, id, vid7, meanUnitErr(*,*)
      NCDF_CLOSE, id
    end
    
    Function Convert2MolC_MML, xCO_PPB, PS, waterVapor
      ; xCO is in ppb
      ; CO ppb = CO molC /  air molC * 1.0e9
      ; air molC = molC * air mass - waterVapor
      ; air mass = delP

      airW = 28.966
      molC = 6.022e23/airW

      delP = PS * 100
      airMass = delP / 9.8
      airMolC = molC * 1.0e3 * airMass
      dryAirMolC = airMolC - waterVapor*1.0e4
      xCO_molC = xCO_PPB * 1.0e-9 * dryAirMolC
      xCO_molC = xCO_molC * 1.0e-4      ; molC/cm^2

      ;print, 'dry airMolC    = ', dryAirMolC
      ;print, 'xCO_ppb        = ', xCO_PPB
      ;print, 'xCO_molC       = ', xCO_molC

      return, xCO_molC
    end

    Function Convert2MolC_PPB, xCO_molC, PS, waterVapor
      ; xCO is in molC/cm2
      ; CO molC = CO ppb * 1.0-9 * air molC
      ; air molC = molC/kg * air mass(kg) - waterVapor
      ; air mass(kg) = delP * 100/9.8

      airW = 28.966
      molC = 6.022e23/airW      ; molC/g/m^2
      waterMass = waterVapor/6.022e23 * 18.0 * 1.0e4 * 1.0e-3

      delP = PS * 100
      airMass = delP / 9.8      ; Kg/m^2
      dryAirMass = airMass - waterMass
      dryAirMolC = molC * 1.0e3 * dryAirMass

      xCO_VMR = (xCO_molC * 1.0e4) / dryAirMolC
      xCO_PPB = xCO_VMR * 1.0e9

      ;print, 'dry airMolC    = ', dryAirMolC
      ;print, 'xCO_molC       = ', xCO_molC
      ;print, 'xCO_ppb        = ', xCO_PPB

      return, xCO_PPB
    end

    Function PressureWeight, P, nLev
      ; pressure weight
      ; pressure is from surface to top
      weight = P
      weight[*] = 0.0
      
      for lev=0,nLev-1 do begin
        p1 = 0.0
        p2 = 0.0
        if (lev lt nLev-1) then begin
          v1 = P(lev) - P(lev+1)
          v2 = P(lev)/P(lev+1)
          p1 = P(lev) - v1/ALOG(v2)
        endif
        if (lev gt 0) then begin
          v1 = P(lev-1)- P(lev)
          v2 = P(lev-1)/P(lev)
          p2 = -P(lev) + v1/ALOG(v2)
        endif
        weight[lev] = 1.0/P[0]* ABS(p1+p2)
      endfor
      weight = weight/TOTAL(weight)
      return, weight
    end
