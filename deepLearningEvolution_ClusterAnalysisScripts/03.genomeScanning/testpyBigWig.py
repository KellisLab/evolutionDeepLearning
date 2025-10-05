import pyBigWig

bw = pyBigWig.open("/home/tnikitha/data/ancoraUpdated/bigWigOutput/altai.Neanderthal.A/altai.Neanderthal.A.A172.bw")
print(bw.header())
print(bw.isBigWig())
