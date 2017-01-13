

import sys
import argparse
import os
from collections import Counter
import gzip

try:
    from StringIO import StringIO # Python 2
except ImportError:
    from io import StringIO # Python 3


from Bio import SeqIO
from intervaltree import IntervalTree
import jellyfish

from search import IgorSuffixTree
from cast import Cast
from utils import *




cd = os.path.dirname(os.path.realpath(__file__))

kmer_len = 3

class Settings(object):
    def __init__(self, **kwargs):
        for name, value in kwargs.items():
            setattr(self, name, value)

    def __str__(self):
        me = ""
        for name, value in self.__dict__.items():
            me += "%s: %s\n" % (name, value)
        return me


class ImReP(object):

    def __init__(self, settings):
        self.__settings = settings
        self._fastq_handle = None

        self.vi_pieces = {}
        self.d_seqs = {}
        self.jay_pieces = {}
        self.pSeq_read_map = {}
        self.just_v = []
        self.just_j = []
        self.just_v_dict = {}
        self.just_j_dict = {}
        self.cdr3_dict = {}
        self.hashV = {}
        self.hashJ = {}
        self.v_chain_type = {}
        self.j_chain_type = {}
        self.__populate_v()
        self.__populate_d()
        self.__populate_j()
        self.__read_reads()

        self.debug_info = {}


    def kmers(self, string, k):
        kmrs = []
        for i in range(len(string) - k + 1):
            kmrs.append(string[i: i + k])
        return kmrs



    def __populate_v(self):
        global cd
        chains_v = map(lambda x: cd + "/db/%s/%sV.faa" % (self.__settings.species, x), self.__settings.chains)
        for ch_v_file in chains_v:
            for record in SeqIO.parse(ch_v_file, "fasta"):
                if "partial in 3'" not in record.description:
                    Vend = str(record.seq)[-20:]
                    kmrs = self.kmers(Vend, kmer_len)
                    for k in kmrs:
                        if k not in self.hashV:
                            self.hashV[k] = set()
                        self.hashV[k].add(record.id)
                    self.v_chain_type[record.id] = getGeneType2(record.id)
                    posC = Vend.rfind("C")
                    if posC != -1:
                        anchor = Vend[:posC]
                        rest = Vend[posC + 1:]
                        self.vi_pieces[record.id] = (anchor, rest)

    def __populate_d(self):
        global cd
        for chain in self.__settings.chains:
            if chain in ["IGH", "TRB", "TRD"]:
                for record in SeqIO.parse(cd + "/db/%s/%sD.faa" % (self.__settings.species, chain), "fasta"):
                    if chain not in self.d_seqs:
                        self.d_seqs[chain] = {}
                    self.d_seqs[chain][record.id] = str(record.seq)


    def __populate_j(self):
        global cd
        chains_j = map(lambda x: cd + "/db/%s/%sJ.faa" % (self.__settings.species, x), self.__settings.chains)
        for ch_j_file in chains_j:
            for record in SeqIO.parse(ch_j_file, "fasta"):
                beginJ = str(record.seq)[:20]
                kmrs = self.kmers(beginJ, kmer_len)
                for k in kmrs:
                    if k not in self.hashJ:
                        self.hashJ[k] = set()
                    self.hashJ[k].add(record.id)
                self.j_chain_type[record.id] = getGeneType2(record.id)
                letter = "F"
                if "IGHJ" in ch_j_file:
                    letter = "W"
                posW = beginJ.find(letter)
                if posW != -1:
                    anchor = beginJ[:posW]
                    rest = beginJ[posW + 1:]
                    self.jay_pieces[record.id] = (anchor, rest)



    def __read_reads(self):
        fastqfile = self.__settings.fastqfile
        formatFile = "fasta"
        if self.__settings.isFastq:
            formatFile = "fastq"
        #if fastqfile.endswith(".fa") or fastqfile.endswith(".fasta") or fastqfile.endswith(".fna") \
        #    or fastqfile.endswith(".fa.gz") or fastqfile.endswith(".fasta.gz") or fastqfile.endswith(".fna.gz"):
        #    formatFile = "fasta"
        #elif fastqfile.endswith(".fq") or fastqfile.endswith(".fastq") \
        #    or fastqfile.endswith(".fq.gz") or fastqfile.endswith(".fastq.gz"):
        #    formatFile = "fastq"
        #else:
        #    raise Exception("Unrecognized format of input file. Please, provide either .fasta or .fastq file")
        if fastqfile.endswith(".gz"):
            with gzip.open(fastqfile, 'rb') as f:
                file_content = f.read()
            self._fastq_handle = SeqIO.parse(StringIO(file_content), formatFile)
        else:
            self._fastq_handle = SeqIO.parse(fastqfile, formatFile)


    def __full_cdr3(self):
        if not self._fastq_handle:
            return []
        vkeys = set(self.hashV.keys())
        jkeys = set(self.hashJ.keys())
        full_cdr3 = []
        for record in self._fastq_handle:
            # If we have paired-end reads,
            # then we have to distinguish them
            if "/1" in record.description:
                record.id += "/1"
            if "/2" in record.description:
                record.id += "/2"
            self.debug_info[record.id] = {}
            pSequences = nucleotide2protein2(str(record.seq))
            if pSequences:
                for pSeq, frame in pSequences:
                    pos1 = pSeq.find("C")
                    pos2 = [pSeq.rfind("F"), pSeq.rfind("W")]
                    v_overlap = "NA"
                    j_overlap = "NA"
                    vtypes = {}
                    jtypes = {}
                    if pos1 != -1:
                        kmrs1 = self.kmers(pSeq[:pos1 + 5], kmer_len)
                        interV = set(kmrs1) & vkeys
                        vlist = []
                        for v in interV:
                            vlist.extend(list(self.hashV[v]))
                        if vlist:
                            vc = [x for x, y in Counter(vlist).items()]
                        else:
                            vc = []
                        v_cl = {}
                        for v in vc:
                            if self.v_chain_type[v] not in v_cl:
                                v_cl[self.v_chain_type[v]] = []
                            v_cl[self.v_chain_type[v]].append(v)
                        f, s = pSeq[:pos1], pSeq[pos1 + 1:]
                        v_overlap = len(f) + len(s) + 1
                        for v1, v2 in v_cl.items():
                            for v3 in v2:
                                if v3 not in self.vi_pieces:
                                    continue
                                v, vv = self.vi_pieces[v3]
                                minlen1 = min(len(f), len(v))
                                minlen2 = min(len(s), len(vv))
                                if minlen1 > 0:
                                    mismatch1 = jellyfish.levenshtein_distance(unicode(f[-minlen1:]), unicode(v[-minlen1:]))
                                else:
                                    mismatch1 = 0
                                if minlen2 > 0:
                                    mismatch2 = jellyfish.levenshtein_distance(unicode(s[:minlen2]), unicode(vv[:minlen2]))
                                else:
                                    mismatch2 = 0
                                if (minlen1 == 0 and mismatch2 <= 1) or (minlen1 > 3 and mismatch1 <= 1 and minlen2 >= 2 and mismatch2 <= 2):
                                    vtypes[v3] = (minlen1 + minlen2 + 1, mismatch1 + mismatch2)
                    if pos2 != [-1, -1]:
                        if pos2[0] != -1:
                            if pos2[1] > 10:
                                offset = pos2[1] - 10
                            else:
                                offset = 0
                            kmrs2 = self.kmers(pSeq[offset:], kmer_len)
                            interJ = set(kmrs2) & jkeys
                            jlist = []
                            for j in interJ:
                                jlist.extend(list(self.hashJ[j]))
                            if jlist:
                                jc = [x for x, y in Counter(jlist).items()]
                            else:
                                jc = []
                            j_cl = {}
                            for j in jc:
                                if self.j_chain_type[j] != "IGHJ" and self.j_chain_type[j] not in j_cl:
                                    j_cl[self.j_chain_type[j]] = []
                                if self.j_chain_type[j] != "IGHJ":
                                    j_cl[self.j_chain_type[j]].append(j)
                            f, s = pSeq[:pos2[0]], pSeq[pos2[0] + 1:]
                            j_overlap = len(f) + len(s) + 1
                            for j1, j2 in j_cl.items():
                                for j3 in j2:
                                    if j3 not in self.jay_pieces:
                                        continue
                                    j, jj = self.jay_pieces[j3]
                                    minlen1 = min(len(f), len(j))
                                    minlen2 = min(len(s), len(jj))
                                    if minlen2 > 0:
                                        mismatch2 = jellyfish.levenshtein_distance(unicode(s[:minlen2]), unicode(jj[:minlen2]))
                                    else:
                                        mismatch2 = 0
                                    if minlen1 > 0:
                                        mismatch1 = jellyfish.levenshtein_distance(unicode(f[-minlen1:]), unicode(j[-minlen1:]))
                                    else:
                                        mismatch1 = 0
                                    if (minlen2 == 0 and mismatch1 <= 1) or (minlen2 > 3 and mismatch2 <= 1 and minlen1 >= 2 and mismatch1 <= 2):
                                        jtypes[j3] = (minlen1 + minlen2 + 1, mismatch1 + mismatch2)
                        if pos2[1] != -1:
                            if pos2[1] > 10:
                                offset = pos2[1] - 10
                            else:
                                offset = 0
                            kmrs2 = self.kmers(pSeq[offset:], kmer_len)
                            interJ = set(kmrs2) & jkeys
                            jlist = []
                            for j in interJ:
                                jlist.extend(list(self.hashJ[j]))
                            if jlist:
                                jc = [x for x, y in Counter(jlist).items()]
                            else:
                                jc = []
                            j_cl = {}
                            for j in jc:
                                if self.j_chain_type[j] == "IGHJ" and self.j_chain_type[j] not in j_cl:
                                    j_cl[self.j_chain_type[j]] = []
                                if self.j_chain_type[j] == "IGHJ":
                                    j_cl[self.j_chain_type[j]].append(j)
                            f, s = pSeq[:pos2[1]], pSeq[pos2[1] + 1:]
                            j_overlap = len(f) + len(s) + 1
                            for j1, j2 in j_cl.items():
                                for j3 in j2:
                                    if j3 not in self.jay_pieces:
                                        continue
                                    j, jj = self.jay_pieces[j3]
                                    minlen1 = min(len(f), len(j))
                                    minlen2 = min(len(s), len(jj))
                                    if minlen2 > 0:
                                        mismatch2 = jellyfish.levenshtein_distance(unicode(s[:minlen2]), unicode(jj[:minlen2]))
                                    else:
                                        mismatch2 = 0
                                    if minlen1 > 0:
                                        mismatch1 = jellyfish.levenshtein_distance(unicode(f[-minlen1:]), unicode(j[-minlen1:]))
                                    else:
                                        mismatch1 = 0
                                    if (minlen2 == 0 and mismatch1 <= 1) or (minlen2 > 3 and mismatch2 <= 1 and minlen1 >= 2 and mismatch1 <= 2):
                                        jtypes[j3] = (minlen1, mismatch1, minlen2, mismatch2)
                    if vtypes or jtypes:
                        vt = {}
                        vscore = {}
                        jt = {}
                        jscore = {}
                        for x in vtypes:
                            chaint = self.v_chain_type[x]
                            if chaint[:3] not in vt:
                                vt[chaint[:3]] = []
                                vscore[chaint[:3]] = []
                            vt[chaint[:3]].append(x)
                            entry = [x] + list(vtypes[x])
                            vscore[chaint[:3]].append(entry)
                        for x in jtypes:
                            chaint = self.j_chain_type[x]
                            if chaint[:3] not in jt:
                                jt[chaint[:3]] = []
                                jscore[chaint[:3]] = []
                            jt[chaint[:3]].append(x)
                            entry = [x] + list(jtypes[x])
                            jscore[chaint[:3]].append(entry)
                        self.debug_info[record.id] = {"vscore": vscore, "jscore": jscore}
                        common = set(vt.keys()) & set(jt.keys())
                        if common:
                            if "IGH" in common:
                                full_cdr3.append(pSeq[pos1: pos2[1] + 1])
                                cdr3 = pSeq[pos1: pos2[1] + 1]
                            else:
                                full_cdr3.append(pSeq[pos1: pos2[0] + 1])
                                cdr3 = pSeq[pos1: pos2[0] + 1]
                            if cdr3 not in self.cdr3_dict:
                                self.cdr3_dict[cdr3] = []
                            self.cdr3_dict[cdr3].append(record.id)
                            if cdr3 not in self.pSeq_read_map or (cdr3 in self.pSeq_read_map and ("v" not in self.pSeq_read_map[cdr3].keys() or "j" not in self.pSeq_read_map[cdr3].keys())):
                                v_t = []
                                j_t = []
                                chtype = {}
                                for key, ch in vt.items():
                                    if key in common:
                                        v_t.extend(ch)
                                        if key not in chtype:
                                            chtype[key] = []
                                        chtype[key].extend(ch)
                                for key, ch in jt.items():
                                    if key in common:
                                        j_t.extend(ch)
                                        if key not in chtype:
                                            chtype[key] = []
                                        chtype[key].extend(ch)
                                self.pSeq_read_map[cdr3] = {"v": map(getGeneType, v_t), "j": map(getGeneType, j_t), "chain_type": chtype}
                        elif vtypes and not jtypes:
                            vi_partial = pSeq[pos1:]
                            if vi_partial not in full_cdr3:
                                self.just_v.append(vi_partial)
                                if vi_partial not in self.just_v_dict:
                                    self.just_v_dict[vi_partial] = []
                                self.just_v_dict[vi_partial].append(record.id)
                            if vi_partial not in self.pSeq_read_map and vi_partial not in full_cdr3:
                                self.pSeq_read_map[vi_partial] = {"v": map(getGeneType, vtypes), "chain_type": vt}
                        elif jtypes and not vtypes:
                            if "IGH" in jt:
                                jay_partial = pSeq[:pos2[1] + 1]
                            else:
                                jay_partial = pSeq[:pos2[0] + 1]
                            if jay_partial not in full_cdr3:
                                self.just_j.append(jay_partial)
                                if jay_partial not in self.just_j_dict:
                                    self.just_j_dict[jay_partial] = []
                                self.just_j_dict[jay_partial].append(record.id)
                            if jay_partial not in self.pSeq_read_map and jay_partial not in full_cdr3:
                                self.pSeq_read_map[jay_partial] = {"j": map(getGeneType, jtypes), "chain_type": jt}
        return full_cdr3



    def __vj_handshakes(self):
        handshakes = []
        just_v = Counter(self.just_v)
        just_j = Counter(self.just_j)

        itree = IntervalTree()

        start = 0
        for v in just_v.keys():
            end = start + len(v) + 1
            itree.addi(start, end, v)
            start = end

        all_v_suf = "|".join(just_v.keys())
        stree = IgorSuffixTree(all_v_suf)

        for j, jj in just_j.items():
            overlap, index, terminal = stree.search_stree(j)
            if terminal and len(j[:overlap]) >= self.__settings.overlapLen:
                overlapping_v = itree.search(index)
                common_chains = set(self.pSeq_read_map[list(overlapping_v)[0].data]["chain_type"].keys()) & set(self.pSeq_read_map[j]["chain_type"].keys())
                if common_chains:
                    v_t = []
                    j_t = []
                    chtype = {}
                    for key, ch in self.pSeq_read_map[list(overlapping_v)[0].data]["chain_type"].items():
                        if key in common_chains:
                            v_t.extend(map(getGeneType, ch))
                            if key not in chtype:
                                chtype[key] = []
                            chtype[key].extend(ch)
                    for key, ch in self.pSeq_read_map[j]["chain_type"].items():
                        if key in common_chains:
                            j_t.extend(map(getGeneType, ch))
                            if key not in chtype:
                                chtype[key] = []
                            chtype[key].extend(ch)
                    newly_born_cdr3 = list(overlapping_v)[0].data + j[overlap:]
                    if newly_born_cdr3 not in self.cdr3_dict:
                        self.cdr3_dict[newly_born_cdr3] = []
                    if list(overlapping_v)[0].data in self.just_v_dict:
                        self.cdr3_dict[newly_born_cdr3].extend(self.just_v_dict[list(overlapping_v)[0].data])
                    if j in self.just_j_dict:
                        self.cdr3_dict[newly_born_cdr3].extend(self.just_j_dict[j])
                    if list(overlapping_v)[0].data in self.just_v_dict:
                        del self.just_v_dict[list(overlapping_v)[0].data]
                    if j in self.just_j_dict:
                        del self.just_j_dict[j]
                    countV = just_v[list(overlapping_v)[0].data]
                    countJ = just_j[j]
                    countVJ = min(countV, countJ)
                    for x in range(countVJ):
                        handshakes.append(newly_born_cdr3)
                    self.pSeq_read_map[newly_born_cdr3] = {"v": v_t, "j": j_t, "chain_type": chtype}
        return handshakes


    def __map_d(self, seq, chain_type):
        d_types = set()
        for d_t, d_seq in self.d_seqs[chain_type].items():
            if seq.find(d_seq) != -1:
                d_types.add(getGeneType(d_t))
        if not d_types:
            return set(["NA"])
        return d_types


    def doComputeClones(self):
        clones = self.__full_cdr3()
        if not self.__settings.noOverlapStep:
            clones2 = self.__vj_handshakes()
            clones.extend(clones2)
        clones = Counter(clones)
        cast_clustering = Cast(clones)
        clustered_clones = cast_clustering.doCast(self.__settings.castThreshold)
        self.clone_dict = {}
        for clone in clustered_clones:
            self.clone_dict[clone[0]] = clone[2]
            del clone[2]
            del clone[1] # remove counts for now
            chain_type = self.pSeq_read_map[clone[0]]["v"][0][:3]
            j_types = None
            if chain_type in ["IGH", "TRB", "TRD"]:
                j_types = self.__map_d(clone[0], chain_type)
            types = [",".join(list(set(self.pSeq_read_map[clone[0]]["v"]))[:3])]
            if j_types:
                types.append(",".join(j_types))
            else:
                types.append("NA")
            types.append(",".join(list(set(self.pSeq_read_map[clone[0]]["j"]))[:3]))
            clone.extend(types)
        return clustered_clones




if __name__ == "__main__":
    ap = argparse.ArgumentParser("python2 imrep.py")

    necessary_arguments = ap.add_argument_group("Necessary Inputs")
    necessary_arguments.add_argument("reads_file", help="unmapped reads in .fasta (default) or .fastq (if flag --fastq is set)  format")
    necessary_arguments.add_argument("output_clones", help="output file with CDR3 clonotypes")

    optional_arguments = ap.add_argument_group("Optional Inputs")
    optional_arguments.add_argument("--fastq", help="a binary flag used to indicate that the input file with unmapped reads is in fastq format", dest="isFastq", action="store_true")
    optional_arguments.add_argument("-s", "--species", help="species (human or mouse, default human)", type=str, dest="species")
    optional_arguments.add_argument("-o", "--overlapLen", help="the minimal length to consider between reads overlapping with a V gene and reads overlapping with a J gene. Default value is 10 amino acids.", type=int)
    optional_arguments.add_argument("--noOverlapStep", help="a binary flag used in case if the user does not want to run the second stage of the ImReP assembly.", dest="noOverlapStep", action="store_true")
    optional_arguments.add_argument("--extendedOutput", help="extended output: write information read by read", dest="extendedOutput", action="store_true")
    optional_arguments.add_argument("-t", "--castThreshold", help="the -t option has been added to control the stringency of CDR3 clustering. The value can be from 0.0 to 1.0. Note that thresholds near 1.0 are more liberal and result in more CDR3 to be reported. The default value is 0.2", type=float)
    optional_arguments.add_argument("-c", "--chains", help="chains: comma separated values from IGH,IGK,IGL,TRA,TRB,TRD,TRG", type=str)

    args = ap.parse_args()

    fastqfile = args.reads_file
    isFastq = args.isFastq
    outFile = args.output_clones

    set_dict = {
        'isFastq': False,
        'species': "human",
        'fastqfile': fastqfile,
        'overlapLen': 10,
        'noOverlapStep': False,
        'extendedOutput': False,
        'castThreshold': 0.2,
        'chains': ['IGH','IGK','IGL','TRA','TRB','TRD','TRG']
    }

    if args.species:
        if args.species in ["human", "mouse"]:
            set_dict["species"] = args.species
        else:
            raise Exception("Species must be either human or mouse")
    if args.overlapLen:
        set_dict["overlapLen"] = args.overlapLen
    if args.noOverlapStep is not None:
        set_dict["noOverlapStep"] = args.noOverlapStep
    if args.isFastq is not None:
        set_dict["isFastq"] = args.isFastq
    if args.extendedOutput is not None:
        set_dict["extendedOutput"] = args.extendedOutput
    if args.castThreshold:
        set_dict["castThreshold"] = args.castThreshold
    if args.chains:
        set_dict["chains"] = args.chains.split(",")

    settings = Settings(**set_dict)

    print "Starting ImReP-0.1"
    imrep = ImReP(settings)
    clones = imrep.doComputeClones()


    final_clones = []
    if set_dict["extendedOutput"]:
        with open("full_cdr3.txt", "w") as f:
            for cl in clones:
                for clon in imrep.clone_dict[cl[0]]:
                    for read in imrep.cdr3_dict[clon]:
                        dinfo_v = imrep.debug_info[read].get("vscore", {})
                        dinfo_j = imrep.debug_info[read].get("jscore", {})
                        di_v = []
                        di_j = []
                        uniq_v, uniq_j = 0, 0
                        for xx, yy in dinfo_v.items():
                            if yy:
                                for u in yy:
                                    geneName = u[0].split("|")[1]
                                    di_v.append(geneName + ":" + ":".join(map(str, u[1:])))
                        for xx, yy in dinfo_j.items():
                            if yy:
                                for u in yy:
                                    geneName = u[0].split("|")[1]
                                    di_j.append(geneName + ":" + ":".join(map(str, u[1:])))
                        if len(di_v) == 1:
                            uniq_v = 1
                        if len(di_j) == 1:
                            uniq_j = 1
                        uniq_vj = uniq_v & uniq_j
                        di_v = ",".join(di_v)
                        if not di_v:
                            di_v = "NA"
                        di_j = ",".join(di_j)
                        if not di_j:
                            di_j = "NA"
                        f.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (read, cl[0], cl[1], cl[2], cl[3], di_v, di_j, uniq_v, uniq_j, uniq_vj))
                        final_clones.append(cl[0] + "\t%s\t" + "%s\t%s\t%s\n" % (cl[1], cl[2], cl[3]))
        with open("partial_cdr3.txt", "w") as f:
            for x, y in imrep.just_v_dict.items():
                for read in y:
                    v = ",".join(list(set(imrep.pSeq_read_map[x].get("v", ["NA"])))[:3])
                    j = ",".join(list(set(imrep.pSeq_read_map[x].get("j", ["NA"])))[:3])
                    dinfo_v = imrep.debug_info[read].get("vscore", {})
                    dinfo_j = imrep.debug_info[read].get("jscore", {})
                    di_v = []
                    di_j = []
                    uniq_v, uniq_j = 0, 0
                    for xx, yy in dinfo_v.items():
                        if yy:
                            for u in yy:
                                geneName = u[0].split("|")[1]
                                di_v.append(geneName + ":" + ":".join(map(str, u[1:])))
                    for xx, yy in dinfo_j.items():
                        if yy:
                            for u in yy:
                                geneName = u[0].split("|")[1]
                                di_j.append(geneName + ":" + ":".join(map(str, u[1:])))
                    if len(di_v) == 1:
                        uniq_v = 1
                    if len(di_j) == 1:
                        uniq_j = 1
                    uniq_vj = uniq_v & uniq_j
                    di_v = ",".join(di_v)
                    if not di_v:
                        di_v = "NA"
                    di_j = ",".join(di_j)
                    if not di_j:
                        di_j = "NA"
                    f.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (read, x, v, "NA", j, di_v, di_j, uniq_v, uniq_j, uniq_vj))
    final_clones = Counter(final_clones)
    print "%s partial-V CDR3 found" % len(imrep.just_v_dict)
    print "%s partial-J CDR3 found" % len(imrep.just_j_dict)
    print "%s full CDR3 found" % len(final_clones)
    clones = []
    for x, y in final_clones.items():
        clones.append(x % y)
    dumpClones2(clones, outFile)
    print "Done. Bye-bye"

