; NCEP
begin
  uname = "u_300_NCEP.nc"
  vname = "v_300_NCEP.nc"
      
  ufile = addfile("NCEP/Data/" + uname, "r")
  vfile = addfile("NCEP/Data/" + vname, "r")

  u = ufile->u_300
  v = vfile->v_300

  printVarSummary(u)
  print(dimsizes(u))
  printVarSummary(v)
  print(dimsizes(v))

  uvmsg = 1e+36

  sf = new(dimsizes(u),float, uvmsg)
  vp = new(dimsizes(u),float, uvmsg)

  uv2sfvpf (u,v, sf,vp)
  ; printMinMax(sf,True)
  ; printMinMax(vp,True)

  oh = addfile("NCEP/Data/" + "sf_vp_300_1948-2022.nc", "c")
  oh->SF = sf
  oh->VP = vp
end
exit



