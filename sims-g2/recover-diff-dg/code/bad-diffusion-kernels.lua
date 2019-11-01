-- kernels
local updateKernels = {}

updateKernels[1] = function(diffCoeff, dxCells, q, qL, qR, qT, qB, qTL, qTR, qBL, qBR, kerOut)
   local dx, dy = dxCells[1], dxCells[2]
   local Dxx, Dxy, Dyy, Dyx = diffCoeff.Dxx, diffCoeff.Dxy, diffCoeff.Dyy, diffCoeff.Dyx

kerOut[1] = (-(2.0*qR[4]*Dyx)/(dx*dy))-(2.0*qL[4]*Dyx)/(dx*dy)+(4.0*q[4]*Dyx)/(dx*dy)+(1.732050807568877*qR[3]*Dyx)/(dx*dy)-(1.732050807568877*qL[3]*Dyx)/(dx*dy)-(2.0*qT[4]*Dxy)/(dx*dy)-(2.0*qB[4]*Dxy)/(dx*dy)+(4.0*q[4]*Dxy)/(dx*dy)+(1.732050807568877*qT[2]*Dxy)/(dx*dy)-(1.732050807568877*qB[2]*Dxy)/(dx*dy)-(2.165063509461096*qT[3]*Dyy)/dy^2+(2.165063509461096*qB[3]*Dyy)/dy^2+(2.25*qT[1]*Dyy)/dy^2+(2.25*qB[1]*Dyy)/dy^2-(4.5*q[1]*Dyy)/dy^2-(2.165063509461096*qR[2]*Dxx)/dx^2+(2.165063509461096*qL[2]*Dxx)/dx^2+(2.25*qR[1]*Dxx)/dx^2+(2.25*qL[1]*Dxx)/dx^2-(4.5*q[1]*Dxx)/dx^2 
kerOut[2] = (-(3.464101615137754*qR[4]*Dyx)/(dx*dy))+(3.464101615137754*qL[4]*Dyx)/(dx*dy)+(3.0*qR[3]*Dyx)/(dx*dy)+(3.0*qL[3]*Dyx)/(dx*dy)+(6.0*q[3]*Dyx)/(dx*dy)+(2.0*qT[3]*Dxy)/(dx*dy)+(2.0*qB[3]*Dxy)/(dx*dy)-(4.0*q[3]*Dxy)/(dx*dy)-(1.732050807568877*qT[1]*Dxy)/(dx*dy)+(1.732050807568877*qB[1]*Dxy)/(dx*dy)-(2.165063509461096*qT[4]*Dyy)/dy^2+(2.165063509461096*qB[4]*Dyy)/dy^2+(2.25*qT[2]*Dyy)/dy^2+(2.25*qB[2]*Dyy)/dy^2-(4.5*q[2]*Dyy)/dy^2-(1.75*qR[2]*Dxx)/dx^2-(1.75*qL[2]*Dxx)/dx^2-(11.5*q[2]*Dxx)/dx^2+(2.165063509461095*qR[1]*Dxx)/dx^2-(2.165063509461095*qL[1]*Dxx)/dx^2 
kerOut[3] = (2.0*qR[2]*Dyx)/(dx*dy)+(2.0*qL[2]*Dyx)/(dx*dy)-(4.0*q[2]*Dyx)/(dx*dy)-(1.732050807568877*qR[1]*Dyx)/(dx*dy)+(1.732050807568877*qL[1]*Dyx)/(dx*dy)-(3.464101615137754*qT[4]*Dxy)/(dx*dy)+(3.464101615137754*qB[4]*Dxy)/(dx*dy)+(3.0*qT[2]*Dxy)/(dx*dy)+(3.0*qB[2]*Dxy)/(dx*dy)+(6.0*q[2]*Dxy)/(dx*dy)-(1.75*qT[3]*Dyy)/dy^2-(1.75*qB[3]*Dyy)/dy^2-(11.5*q[3]*Dyy)/dy^2+(2.165063509461095*qT[1]*Dyy)/dy^2-(2.165063509461095*qB[1]*Dyy)/dy^2-(2.165063509461096*qR[4]*Dxx)/dx^2+(2.165063509461096*qL[4]*Dxx)/dx^2+(2.25*qR[3]*Dxx)/dx^2+(2.25*qL[3]*Dxx)/dx^2-(4.5*q[3]*Dxx)/dx^2 
kerOut[4] = (3.464101615137754*qR[2]*Dyx)/(dx*dy)-(3.464101615137754*qL[2]*Dyx)/(dx*dy)-(3.0*qR[1]*Dyx)/(dx*dy)-(3.0*qL[1]*Dyx)/(dx*dy)+(6.0*q[1]*Dyx)/(dx*dy)+(3.464101615137754*qT[3]*Dxy)/(dx*dy)-(3.464101615137754*qB[3]*Dxy)/(dx*dy)-(3.0*qT[1]*Dxy)/(dx*dy)-(3.0*qB[1]*Dxy)/(dx*dy)+(6.0*q[1]*Dxy)/(dx*dy)-(1.75*qT[4]*Dyy)/dy^2-(1.75*qB[4]*Dyy)/dy^2-(11.5*q[4]*Dyy)/dy^2+(2.165063509461095*qT[2]*Dyy)/dy^2-(2.165063509461095*qB[2]*Dyy)/dy^2-(1.75*qR[4]*Dxx)/dx^2-(1.75*qL[4]*Dxx)/dx^2-(11.5*q[4]*Dxx)/dx^2+(2.165063509461095*qR[3]*Dxx)/dx^2-(2.165063509461095*qL[3]*Dxx)/dx^2 
   
end

return updateKernels
