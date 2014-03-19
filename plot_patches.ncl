load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl" ;copyatt

; Plots patches, read from ascii, on Lambert Conformal map projection

begin ; PRESETS
  fnames = (/"xy_01.txt","xy_02.txt","xy_03.txt","xy_04.txt"/)
  headers= (/"xyref_01.txt","xyref_02.txt","xyref_03.txt","xyref_04.txt"/)

  ;fname1 = "lvl01_trace.txt"
  ;fname2 = "lvl02_trace.txt"
  ;fname3 = "lvl03_trace.txt"
  
  ncheader = "/home/adam/data/FLX_MERGED/yearly/nexc_d02_VEI7_200803-200902.nc"

  ;header1 = "lvl01_x0y0_x1y1.txt"
  ;header2 = "lvl02_x0y0_x1y1.txt"
  ;header3 = "lvl03_x0y0_x1y1.txt"
  
  outprefix="dacre_trace_"

  ; The coordinate used to convert pixels to lat-lon:
  xmin_in_lon = todouble(-30.)
  ymin_in_lat = todouble(45.)

  xmax_in_lon = todouble(40.)
  ymax_in_lat = todouble(70.)

  ; Output type
  ;type = "x11"
  ;type = "pdf"
  ;type "ps"
  type = "eps"
  
  ; Map projection settings as in WRF grid (Lambert Conformal is a very special case!)
  res = True
  res@mpProjection = "LambertConformal"
  res@mpLambertMeridianF  = -3.0
  res@mpLambertParallel1F = 50.0
  res@mpLambertParallel2F = 70.0
  ; More map projection under misc plot settings (should be automatic)
end   ;PRESETS

; Function to convert xy-coordinate (on a lat-lon grid) to lat-lon coordinates
function xy2ll(xy,x0,y0,x1,y1,lon0,lat0,lon1,lat1)
local x2lon,y2lat ;variable list
begin
  ; Get factors to convert xy-distance to ll-distance:
  ;printVarSummary(lon1)
  ;printVarSummary(x1)
  dlondx = (lon1-lon0)/(x1-x0)  ; distance
  dlatdy = (lat1-lat0)/(y1-y0)  ;
  ; Convert coordinate to relative (lon,lat)=(0,0) and to lat-lon units
  ll = xy
  ll(:,0) = (xy(:,0)-x0)*dlondx+lon0
  ll(:,1) = (xy(:,1)-y0)*dlatdy+lat0
  return ll
end

begin ; Read lat2d,lon2d from netcdf file
  f    = addfile(ncheader,"r")
  lat2d = f->XLAT  ; Native latitude (2D)
  lon2d = f->XLONG ; Native longitude (2D)
  nxy   = dimsizes(lat2d)
  nlat  = nxy(0)
  nlon  = nxy(1)
end
begin ; Polygon/Polyline resources:
  pgres                  = True
  pgres@gsEdgesOn        = True
  pgres@gsLineColor      = "black"
  pgres@gsLineThicknessF = 2.0
  pgres@gsFillColor      = "tan"
end

begin ; misc plot resources:
  res@gsnDraw                     = False
  res@gsnFrame                    = False

  ; Coast-lines
  res@mpDataBaseVersion           = "Ncarg4_1"
  res@mpDataBaseVersion           = "MediumRes"          ; Finer coast-lines
  ;res@mpGeophysicalLineColor      = "black"
  res@mpGeophysicalLineThicknessF = 0.5
  ;res@mpFillColors                = "black"
  ;res@mpFillOn                    = True                ; turn on map fill
  ;res@mpAreaMaskingOn             = True
  ;res@mpLandFillColor             = "tan"
  res@mpOutlineOn                 = True
  res@mpOutlineDrawOrder          = "PostDraw"
  res@mpFillOn                    = False

  ; Grid and ticks
  res@pmTickMarkDisplayMode       = "Always"             ; turn on tick marks
  res@mpGridAndLimbOn             = True
  res@mpGridSpacingF              = 5.0
  res@mpGridLineThicknessF        = 0.1
  res@mpGridLineColor             = "black"


  ;res@tiMainString                = var2d@description
  ;res@tiMainString                = "Regional Domain on Stereographic Proj"

  ;res@gsnAddCyclic                = False                ; regional data

  ; Zoom in on world map (a bit ugly, could be done better, e.g. manually)
  res@mpLimitMode               = "Corners"           ; method to zoom
  res@mpLeftCornerLatF          = lat2d(0,0)-(lat2d(0,0)-lat2d(0,0))
  res@mpLeftCornerLonF          = lon2d(0,0)-(lon2d(0,0)-lon2d(0,0))
  res@mpRightCornerLatF         = 2.*lat2d(nlat-1,nlon-1)-lat2d(nlat-1,nlon-1)
  res@mpRightCornerLonF         = 2.*lon2d(nlat-1,nlon-1)-lon2d(nlat-1,nlon-1)

  ;res@tfDoNDCOverlay              = True                 ; don't transform 
  res@tfDoNDCOverlay              = False ; Yes, please!
end ;misc plot settings

begin ; MAIN
  ncol = 2
  n   = dimsizes(fnames)
  i = 0
  do while(i.lt.n(0))
    ; Read data from files/headers:
    cxy = readAsciiTable(fnames(i),ncol,"double",0)
    cminmax = readAsciiTable(headers(i),ncol,"double",0)
    cxmin = cminmax(0,0)
    cxmax = cminmax(1,0)
    cymin = cminmax(0,1)
    cymax = cminmax(1,1)
    
    ; Convert xy-coord to ll-coord:
    cll = xy2ll(cxy,cxmin,cymin,cxmax,cymax,xmin_in_lon,ymin_in_lat,xmax_in_lon,ymax_in_lat)

    ; Create work station and define colors:
    outname = outprefix+sprinti("%0.2i",i+1)
    wks  = gsn_open_wks (type, outname)
    gsn_define_colormap(wks,(/"white","black","tan","LightBlue","Blue"/))

    ; Plot polygons on map:
    map = gsn_csm_map(wks,res)
    dum = gsn_add_polygon(wks,map,cll(:,0),cll(:,1),pgres)

    draw(map)
    frame(wks)

    ; Reset variables (next iteration will probably use different sizes)
    delete([/cxy,cll/])

    i=i+1
  end do
end ; MAIN
