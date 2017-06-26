import sys
import random
from subprocess import check_output
def randomSNVAdder(fastaName, nMuts):
  fasta = ""
  with open(fastaName) as f:
    fasta = fasta.join([l.strip() for l in f if '>' not in l])
    fasta = list(fasta)
  mutationList = []

  for i in range(int(nMuts)):
    snv = ()
    while True:
      snvSite = random.randint(0, len(fasta))
      if (snvSite % 2):
        snv = (snvSite, fasta[snvSite], transversion_mutation(fasta[snvSite]))
      else:
        snv = (snvSite, fasta[snvSite], transition_mutation(fasta[snvSite]))
      # call bed to make sure SNV is outside of centromere
      with open('snv.bed', 'w') as f:
        f.write('chr22\t%d\t%d\n' % (snv[0], snv[0]))
      bedRes = check_output(['./bedtools2/bin/bedtools intersect -v -a ./snv.bed -b ./centromeres.bed'],shell=True)
      if bedRes == '':
        print "Failed SNV %d" %(snv[0])
      if snv not in mutationList and snv[2] != 'N' and bedRes != '':
        print "adding %d: %s -> %s" % (snv[0], snv[1], snv[2])
        break

    mutationList.append(snv)
  for snv in mutationList:
    fasta[snv[0]] = snv[2]
  return "".join([c for c in fasta]), sorted(mutationList, key=lambda x:x[0])

def writeFastq(header, fasta, ofname, lineLength):
  of = open(ofname, 'w')
  of.write(header + '\n')
  out = '\n'.join([fasta[i:i+int(lineLength)] for i in range(0, len(fasta), int(lineLength))])
  of.write(out)
  of.close()

def writeSNVs(snvList, ofname):
  of  = open(ofname, 'w')
  of.write('\n'.join(
    ["PointMutation\tchr22\t%d\t%s\t%s" % (site+1, h, t) for site, h, t in snvList]
  ))
  of.close()

def transversion_mutation(base):
    return {
        "A":"T",
        "T":"A",
        "C":"G",
        "G":"C",
        "N":"N"
    }[base]

def transition_mutation(base):
    return {
        "A":"G",
        "G":"A",
        "T":"C",
        "C":"T",
        "N":"N"
    }[base]

def main():
  if len(sys.argv) != 6:
    print "Usage: <exe> <fa> <nmuts> <o_fa_name> <line_len> <o_m_list_name>"
    sys.exit(1)

  fasta, mutationList = randomSNVAdder(sys.argv[1], sys.argv[2])
  f=open(sys.argv[1])
  writeFastq(f.readlines()[0], fasta, sys.argv[3], sys.argv[4])
  f.close()
  writeSNVs(mutationList, sys.argv[5])

if '__main__' == __name__:
  main()
