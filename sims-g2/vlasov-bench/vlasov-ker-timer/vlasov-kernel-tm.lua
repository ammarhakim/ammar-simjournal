local Lin = require "Lib.Linalg"
local Time = require "Lib.Time"
local VlasovModDecl = require "Updater.vlasovData.VlasovModDecl"
local ffi = require "ffi"
local Basis = require "Basis"

function timeKernel(bnm, polyOrder, cdim, vdim, nloop)
   local ndim = cdim+vdim

   local confBasis, phaseBasis
   if bnm == "ms" then
      confBasis = Basis.CartModalSerendipity {  ndim = cdim, polyOrder = polyOrder }
      phaseBasis = Basis.CartModalSerendipity {  ndim = cdim+vdim, polyOrder = polyOrder }
   elseif bnm == "mo" then
      confBasis = Basis.CartModalMaxOrder {  ndim = cdim, polyOrder = polyOrder }
      phaseBasis = Basis.CartModalMaxOrder {  ndim = cdim+vdim, polyOrder = polyOrder }
   end

   local nConf, nPhase = confBasis:numBasis(), phaseBasis:numBasis()
   -- allocate memory
   local w, dxv = Lin.Vec(ndim), Lin.Vec(ndim)
   local E = Lin.Vec(8*nConf)
   local fIn, fOut = Lin.Vec(nPhase), Lin.Vec(nPhase)
   local fInL, fOutL = Lin.Vec(nPhase), Lin.Vec(nPhase)
   local fInR, fOutR = Lin.Vec(nPhase), Lin.Vec(nPhase)

   local tmStart, tmEnd
   
   print(string.format("Basis %s polyOrder %d cdim %d vdim %d", bnm, polyOrder, cdim, vdim))
   print(string.format("Num phase-basis %d", nPhase))

   -- streaming updates
   volStreamUpdate = VlasovModDecl.selectVolStream(phaseBasis:id(), cdim, vdim, polyOrder)
   surfStreamUpdate = VlasovModDecl.selectSurfStream(phaseBasis:id(), cdim, vdim, polyOrder)

   tmStart = Time.clock()
   for i = 1, nloop do
      volStreamUpdate(w:data(), dxv:data(), fIn:data(), fOut:data())
   end
   tmEnd = Time.clock()
   local vs = tmEnd-tmStart
   print(string.format("Volume-stream %g secs", tmEnd-tmStart))

   tmStart = Time.clock()   
   for i = 1, nloop do
      surfStreamUpdate[1](w:data(), dxv:data(), fInL:data(), fInR:data(), fOutL:data(), fOutR:data())
   end
   tmEnd = Time.clock()
   local ss = tmEnd-tmStart
   
   print(string.format("Surface-stream %g secs", tmEnd-tmStart))
   print(string.format("Full stream updates should take %g secs", vs+cdim*ss))

   -- force updates
   volForceUpdate = VlasovModDecl.selectVolElcMag(phaseBasis:id(), cdim, vdim, polyOrder)
   surfForceUpdate = VlasovModDecl.selectSurfElcMag(phaseBasis:id(), cdim, vdim, polyOrder)

   tmStart = Time.clock()
   for i = 1, nloop do
      volForceUpdate(w:data(), dxv:data(), E:data(), fIn:data(), fOut:data())
   end
   tmEnd = Time.clock()
   local vf = tmEnd-tmStart
   print(string.format("Volume-force %g secs", tmEnd-tmStart))

   tmStart = Time.clock()
   for i = 1, nloop do
      surfForceUpdate[1](w:data(), dxv:data(), 1.0, E:data(), fInL:data(), fInR:data(), fOutL:data(), fOutR:data())
   end
   tmEnd = Time.clock()
   local sf = tmEnd-tmStart
   
   print(string.format("Surface-force %g secs", tmEnd-tmStart))
   print(string.format("Full force update should take %g secs", vf+vdim*sf))

   print(string.format("Total update should take %g per cell", (vs+cdim*ss + vf+vdim*sf)/nloop))
   print(string.format("Number of updates %d\n\n", nloop))
end

timeKernel("ms", 2, 1, 1, 1e8)
--timeKernel("ms", 2, 1, 2, 1e6)
--timeKernel("ms", 2, 1, 3, 1e6)
--timeKernel("ms", 2, 2, 2, 1e6)
--timeKernel("ms", 2, 2, 3, 1e6)

--timeKernel("mo", 2, 1, 1, 1e7)
--timeKernel("mo", 2, 1, 2, 1e7)
--timeKernel("mo", 2, 1, 3, 1e7)
--timeKernel("mo", 2, 2, 2, 1e7)
--timeKernel("mo", 2, 2, 3, 1e7)


