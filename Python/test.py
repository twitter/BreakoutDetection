import breakout_detection

Z = []

for i in xrange(64):
    Z.append(0.0)

for i in xrange(64):
    Z.append(1.0)

edm_per = breakout_detection.EdmPer()
edm_per.evaluate(Z)
print edm_per.getLoc()

edm_multi = breakout_detection.EdmMulti()
edm_multi.evaluate(Z)
print edm_multi.getLoc()

edm_tail = breakout_detection.EdmTail()
edm_tail.evaluate(Z)
print edm_tail.getLoc()

edmx = breakout_detection.Edmx()
edmx.evaluate(Z)
print edmx.getLoc()

