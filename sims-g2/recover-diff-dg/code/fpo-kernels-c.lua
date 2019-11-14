-- kernels
local updateKernels = {}
local ffi = require "ffi"
local Lin = require "Lib.Linalg"

fpoKernelLib = ffi.load("fpo-kernels")

ffi.cdef [[
    void updateKernelP1(double dxCells[2],
      double *qTL, double *qT, double *qTR, double *qL, double *q, double *qR, double *qBL, double *qB, double *qBR,
      double *gTL, double *gT, double *gTR, double *gL, double *g, double *gR, double *gBL, double *gB, double *gBR,
      double *kerOut);
]]

updateKernels[1] = function(
      dxCells,
      qTL, qT, qTR, qL, q, qR, qBL, qB, qBR,
      gsTL, gsT, gsTR, gsL, gs, gsR, gsBL, gsB, gsBR,
      kerOut)

   local dxVec = Lin.Vec(2)
   dxVec[1], dxVec[2] = dxCells[1], dxCells[2]

   fpoKernelLib.updateKernelP1(
      dxVec:data(),
      qTL:data(), qT:data(), qTR:data(), qL:data(), q:data(), qR:data(), qBL:data(), qB:data(), qBR:data(),
      gsTL:data(), gsT:data(), gsTR:data(), gsL:data(), gs:data(), gsR:data(), gsBL:data(), gsB:data(), gsBR:data(),
      kerOut:data())
end


return updateKernels
