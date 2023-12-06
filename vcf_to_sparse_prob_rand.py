import sys
import numpy as np
import gzip

inText = sys.argv[1]

if inText.endswith(".gz"):
        inFile = gzip.open(inText, "rt")
        outFile = gzip.open(sys.argv[2],"wt")
else:
        inFile = open(inText,"r")
        outFile = open(sys.argv[2],"w")

vcfOffset = 9

count = 0
for line in inFile:

        count+=1
        if count % 10000 == 0:
                print( "read %s variants" % ( count))

        if line.startswith("##"):
                continue

        if line.startswith("#"):
                header = line.strip().split("\t")
                continue

        line = line.strip()
        spl = line.split("\t")

        inf = spl[0:vcfOffset]
        gt = np.array(spl[vcfOffset:])

        if ":" not in gt[0]:
                # exclude entries without probability score
                continue

        nonZero = np.where(gt != "0|0:.")[0]
        wanted = gt[nonZero]

        inds11 = []
        inds10 = []
        inds01 = []
        inds11Prob = []
        inds10Prob = []
        inds01Prob = []
        for i in range(len(wanted)):
                t = wanted[i].split(":")
                f = t[0]
                try:
                        p = t[1]
                except:
                        print("Unexpected genotype:", wanted[i], inf, "skipping..")
                        continue

                if p == ".":
                        continue
                if f == "0|0":
                        continue

                if f == "1|1":
                        ind = header[nonZero[i] + vcfOffset]
                        inds11.append(ind)
                        inds11Prob.append(p)
                else:
                        rand = np.random.rand(1)
                        if rand > 0.5:
                                ind = header[nonZero[i] + vcfOffset]
                                inds10.append(ind)
                                inds10Prob.append(p)
                        else:
                                ind = header[nonZero[i] + vcfOffset]
                                inds01.append(ind)
                                inds01Prob.append(p)

        if len(inds11) > 0:
                text = "%s\t1|1\t%s\t%s\n" % ("\t".join(inf),";".join(inds11),";".join(inds11Prob))
                outFile.write(text)
        if len(inds10) > 0:
                text = "%s\t1|0\t%s\t%s\n" % ("\t".join(inf),";".join(inds10),";".join(inds10Prob))
                outFile.write(text)
        if len(inds01) > 0:
                text = "%s\t0|1\t%s\t%s\n" % ("\t".join(inf),";".join(inds01),";".join(inds01Prob))
                outFile.write(text)

outFile.close()
print("Finished!")
