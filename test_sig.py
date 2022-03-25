import class_ABCD
import math

pred=73.817
pred_error=1.689
sig=1

x1, x2 = class_ABCD.calcMinMaxAllowedObs(sig, pred, pred_error)
print(x1,x2)
print(class_ABCD.calcNonClosureSig(pred, pred_error, x1, math.sqrt(x1)), class_ABCD.calcNonClosureSig(pred, pred_error, x2, math.sqrt(x2) ))

c1, c2 = class_ABCD.calcMinMaxAllowedClosure(sig, pred, pred_error)
print(c1, c2)
maxVal=max(abs(1.0 - c1), abs(1.0-c2))
print(maxVal)