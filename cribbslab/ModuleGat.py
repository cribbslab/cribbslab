
import cgatcore.iotools as iotools

def gene2trans(infile, outfile, transgenemap):


    mapping = iotools.open_file(transgenemap)

    mapped = {}
    for line in mapping:
        split = line.split()
        key = split[1]
        value = split[0]
        mapped[key] = value
    mapping.close()

    active_gene = iotools.open_file(infile)
    active_trans = iotools.open_file(outfile, "w")
    for line in active_gene:
        line = line.rstrip()
        try:
            trans = mapped[line]
        except KeyError:
            print("Not there")
            pass
        active_trans.write("%s\n" % (trans))
    active_gene.close()
    active_trans.close()


def replace(file, out):
    for line in file:
        chrom,start,end,gene,value, strand = line.split()
        start = int(start) - 2000
        end = int(end) + 2000
    
        final = "\t".join([chrom,str(start),str(end),"CDS",value,strand])
        out.write(final + "\n")


def extend_bed(infile, outfile):
    inf = iotools.open_file(infile)

    outf = iotools.open_file(outfile, "w")

        
    replace(inf, outf)

    outf.close()
