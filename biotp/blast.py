def slice_rgo_by_hgt(input_filename, output_filename, pct):
    qry = {}

    with open(input_filename, 'r') as input_handle:
        for line in input_handle:
            if line.startswith('#'):
                continue

            li = line.strip().split('\t')
            if len(li) != 12:
                continue

            record = {
                'sseqid': sseqid, 
                'pident': pident, 
                'length': length, 
                'mismatch': mismatch, 
                'gapopen': gapopen, 
                'qstart': qstart, 
                'qend': qend, 
                'sstart': sstart, 
                'send': send, 
                'evalue': evalue, 
                'bitscore': bitscore
            }

            qseqid, sseqid, pident, length, mismatch, gapopen, qstart, qend, sstart, send, evalue, bitscore = li

            if qseqid not in qry:
                qry[qseqid] = {
                    'sseqid': [], 
                    'pident': [], 
                    'length': [], 
                    'mismatch': [], 
                    'gapopen': [], 
                    'qstart': [], 
                    'qend': [], 
                    'sstart': [], 
                    'send': [], 
                    'evalue': [], 
                    'bitscore': [], 
                    'rec_bits': [], 
                    'grp_bits': [], 
                    'ogp_bits': []
                }

            for key, value in record.items():
                qry[qseqid][key].append(value)

            if sseqid.startswith('rec'):
                qry[qseqid]['rec_bits'].append(bitscore)
            elif sseqid.startswith('grp'):
                qry[qseqid]['grp_bits'].append(bitscore)
            elif sseqid.startswith('ogp'):
                qry[qseqid]['ogp_bits'].append(bitscore)

    with open(output_filename, 'w') as output_handle:
        for q in qry:
            ##-- Calculate alien-index --##
            if not qry[q]['ogp_bits']:
                continue

            ## Bitscore of besthit in outgroup
            bbh_ogp = max(qry[q]['ogp_bits'])

            ## Bitscore of besthit in recipient
            if not qry[q]['rec_bits']:
                bbh_rec = max(qry[q]['bitscore'])
            else:
                bbh_rec = max(qry[q]['rec_bits'])

            ## Bitscore of besthit in group
            if not qry[q]['grp_bits']:
                bbh_grp = min(qry[q]['bitscore'])
            else:
                bbh_grp = max(qry[q]['grp_bits'])

            ai = (float(bbh_ogp) / float(bbh_rec)) - (float(bbh_grp) / float(bbh_rec))

            ##-- Calculate percentage of outgroup species --##
            num_all = len(qry[q]['bitscore']) ## Number of all hits
            num_ogp = len(qry[q]['ogp_bits']) ## Number of ogp hits
            ogp_pct = int(num_ogp) / int(num_all) * 100

            if float(ai) > 0 and float(ogp_pct) >= float(pct):
                for sseqid, pident, length, mismatch, gapopen, qstart, qend, sstart, send, evalue, bitscore in zip(
                        qry[q]['sseqid'], qry[q]['pident'], qry[q]['length'], qry[q]['mismatch'], qry[q]['gapopen'], 
                        qry[q]['qstart'], qry[q]['qend'], qry[q]['sstart'], qry[q]['send'], 
                        qry[q]['evalue'], qry[q]['bitscore']):
                    output_handle.write(
                        f"{q}\t{sseqid}\t{pident}\t{length}\t{mismatch}\t{gapopen}\t{qstart}\t{qend}\t{sstart}\t{send}\t{evalue}\t{bitscore}\n")
                    
def slice_results_by_rgo(input_filename, output_filename, pct):
    qry = {}

    with open(input_filename, 'r') as input_handle:
        for line in input_handle:
            if line.startswith('#'):
                continue

            li = line.strip().split('\t')
            if len(li) != 12:
                continue

            record = {
                'sseqid': sseqid, 
                'pident': pident, 
                'length': length, 
                'mismatch': mismatch, 
                'gapopen': gapopen, 
                'qstart': qstart, 
                'qend': qend, 
                'sstart': sstart, 
                'send': send, 
                'evalue': evalue, 
                'bitscore': bitscore
            }

            qseqid, sseqid, pident, length, mismatch, gapopen, qstart, qend, sstart, send, evalue, bitscore = li

            if qseqid not in qry:
                qry[qseqid] = {
                    'sseqid': [], 
                    'pident': [], 
                    'length': [], 
                    'mismatch': [], 
                    'gapopen': [], 
                    'qstart': [], 
                    'qend': [], 
                    'sstart': [], 
                    'send': [], 
                    'evalue': [], 
                    'bitscore': [], 
                    'rec_bits': [], 
                    'grp_bits': [], 
                    'ogp_bits': []
                }

            for key, value in record.items():
                qry[qseqid][key].append(value)

            if sseqid.startswith('rec'):
                qry[qseqid]['rec_bits'].append(bitscore)
            elif sseqid.startswith('grp'):
                qry[qseqid]['grp_bits'].append(bitscore)
            elif sseqid.startswith('ogp'):
                qry[qseqid]['ogp_bits'].append(bitscore)

    with open(output_filename, 'w') as output_handle:
        for q in qry:
            ##-- Calculate alien-index --##
            if not qry[q]['ogp_bits']:
                continue

            ## Bitscore of besthit in outgroup
            bbh_ogp = max(qry[q]['ogp_bits'])

            ## Bitscore of besthit in recipient
            if not qry[q]['rec_bits']:
                bbh_rec = max(qry[q]['bitscore'])
            else:
                bbh_rec = max(qry[q]['rec_bits'])

            ## Bitscore of besthit in group
            if not qry[q]['grp_bits']:
                bbh_grp = min(qry[q]['bitscore'])
            else:
                bbh_grp = max(qry[q]['grp_bits'])

            ai = (float(bbh_ogp) / float(bbh_rec)) - (float(bbh_grp) / float(bbh_rec))

            ##-- Calculate percentage of outgroup species --##
            num_all = len(qry[q]['bitscore']) ## Number of all hits
            num_ogp = len(qry[q]['ogp_bits']) ## Number of ogp hits
            ogp_pct = int(num_ogp) / int(num_all) * 100

            if float(ai) > 0 and float(ogp_pct) >= float(pct):
                for sseqid, pident, length, mismatch, gapopen, qstart, qend, sstart, send, evalue, bitscore in zip(
                        qry[q]['sseqid'], qry[q]['pident'], qry[q]['length'], qry[q]['mismatch'], qry[q]['gapopen'], 
                        qry[q]['qstart'], qry[q]['qend'], qry[q]['sstart'], qry[q]['send'], 
                        qry[q]['evalue'], qry[q]['bitscore']):
                    output_handle.write(
                        f"{q}\t{sseqid}\t{pident}\t{length}\t{mismatch}\t{gapopen}\t{qstart}\t{qend}\t{sstart}\t{send}\t{evalue}\t{bitscore}\n")


def detect_hgt(input_filename, output_filename):
    qry = {}

    with open(input_filename, 'r') as input_handle:
        for line in input_handle:
            if line.startswith('#'):
                continue
            li = line.strip().split('\t')
            if len(li) != 12:
                continue

            qseqid, sseqid, pident, length, mismatch, gapopen, qstart, qend, sstart, send, evalue, bitscore = li

            if qseqid not in qry:
                qry[qseqid] = {
                    'sseqid': [],
                    'pident': [],
                    'length': [],
                    'mismatch': [],
                    'gapopen': [],
                    'qstart': [],
                    'qend': [],
                    'sstart': [],
                    'send': [],
                    'evalue': [],
                    'bitscore': [],
                    'rec_bits': [],
                    'grp_bits': [],
                    'ogp_bits': []
                }

            qry[qseqid]['sseqid'].append(sseqid)
            qry[qseqid]['pident'].append(pident)
            qry[qseqid]['length'].append(length)
            qry[qseqid]['mismatch'].append(mismatch)
            qry[qseqid]['gapopen'].append(gapopen)
            qry[qseqid]['qstart'].append(qstart)
            qry[qseqid]['qend'].append(qend)
            qry[qseqid]['sstart'].append(sstart)
            qry[qseqid]['send'].append(send)
            qry[qseqid]['evalue'].append(evalue)
            qry[qseqid]['bitscore'].append(bitscore)

            if sseqid.startswith('rec'):
                qry[qseqid]['rec_bits'].append(bitscore)
            elif sseqid.startswith('grp'):
                qry[qseqid]['grp_bits'].append(bitscore)
            elif sseqid.startswith('ogp'):
                qry[qseqid]['ogp_bits'].append(bitscore)

    with open(output_filename, 'w') as output_handle:
        for q in qry:
            if qry[q]['grp_bits']:
                continue
            if qry[q]['ogp_bits']:
                for sseqid, pident, length, mismatch, gapopen, qstart, qend, sstart, send, evalue, bitscore in zip(
                        qry[q]['sseqid'], qry[q]['pident'], qry[q]['length'], qry[q]['mismatch'],
                        qry[q]['gapopen'], qry[q]['qstart'], qry[q]['qend'], qry[q]['sstart'],
                        qry[q]['send'], qry[q]['evalue'], qry[q]['bitscore']):
                    output_handle.write(
                        f"{q}\t{sseqid}\t{pident}\t{length}\t{mismatch}\t{gapopen}\t{qstart}\t{qend}\t{sstart}\t{send}\t{evalue}\t{bitscore}\n")


def slice_results_by_rgo_2(input_filename, output_filename, pct):
    qry = {}

    with open(input_filename, 'r') as input_handle:
        for line in input_handle:
            if line.startswith('#'):
                continue

            li = line.strip().split('\t')
            if len(li) != 12:
                continue

            record = {
                'sseqid': sseqid, 
                'pident': pident, 
                'length': length, 
                'mismatch': mismatch, 
                'gapopen': gapopen, 
                'qstart': qstart, 
                'qend': qend, 
                'sstart': sstart, 
                'send': send, 
                'evalue': evalue, 
                'bitscore': bitscore
            }

            qseqid, sseqid, pident, length, mismatch, gapopen, qstart, qend, sstart, send, evalue, bitscore = li

            if qseqid not in qry:
                qry[qseqid] = {
                    'sseqid': [], 
                    'pident': [], 
                    'length': [], 
                    'mismatch': [], 
                    'gapopen': [], 
                    'qstart': [], 
                    'qend': [], 
                    'sstart': [], 
                    'send': [], 
                    'evalue': [], 
                    'bitscore': [], 
                    'rec_bits': [], 
                    'grp_bits': [], 
                    'ogp_bits': []
                }

            for key, value in record.items():
                qry[qseqid][key].append(value)

            if sseqid.startswith('rec'):
                qry[qseqid]['rec_bits'].append(bitscore)
            elif sseqid.startswith('grp'):
                qry[qseqid]['grp_bits'].append(bitscore)
            elif sseqid.startswith('ogp'):
                qry[qseqid]['ogp_bits'].append(bitscore)

    with open(output_filename, 'w') as output_handle:
        for q in qry:
            ##-- Calculate alien-index --##
            if not qry[q]['ogp_bits']:
                continue

            ## Bitscore of besthit in outgroup
            bbh_ogp = max(qry[q]['ogp_bits'])

            ## Bitscore of besthit in recipient
            if not qry[q]['rec_bits']:
                bbh_rec = max(qry[q]['bitscore'])
            else:
                bbh_rec = max(qry[q]['rec_bits'])

            ## Bitscore of besthit in group
            if not qry[q]['grp_bits']:
                bbh_grp = min(qry[q]['bitscore'])
            else:
                bbh_grp = max(qry[q]['grp_bits'])

            ai = (float(bbh_ogp) / float(bbh_rec)) - (float(bbh_grp) / float(bbh_rec))

            ##-- Calculate percentage of outgroup species --##
            num_all = len(qry[q]['bitscore']) ## Number of all hits
            num_ogp = len(qry[q]['ogp_bits']) ## Number of ogp hits
            ogp_pct = int(num_ogp) / int(num_all) * 100

            if float(ai) > 0 and float(ogp_pct) >= float(pct):
                for sseqid, pident, length, mismatch, gapopen, qstart, qend, sstart, send, evalue, bitscore in zip(
                        qry[q]['sseqid'], qry[q]['pident'], qry[q]['length'], qry[q]['mismatch'], qry[q]['gapopen'], 
                        qry[q]['qstart'], qry[q]['qend'], qry[q]['sstart'], qry[q]['send'], 
                        qry[q]['evalue'], qry[q]['bitscore']):
                    output_handle.write(
                        f"{q}\t{sseqid}\t{pident}\t{length}\t{mismatch}\t{gapopen}\t{qstart}\t{qend}\t{sstart}\t{send}\t{evalue}\t{bitscore}\n")
                    

def slice_hits_by_crossover_group(input_filename: str, output_filename: str) -> None:
    """Slice BLAST hits by crossover group

    Args
    ----
    input_filename : str
        Input filename.
    output_filename : str
        Output filename.
    
    """
    hits = {}

    with open(input_filename, "r") as input_handle:
        for line in input_handle:
            if line.startswith("#"):
                continue
            li = line.strip().split("\t")
            if len(li) != 12:
                continue
            qseqid, sseqid, pident, length, mismatch, gapopen, qstart, qend, sstart, send, evalue, bitscore = li
            if qseqid not in hits:
                hits[qseqid] = []
            hits[qseqid].append(li)

    with open(output_filename, "w") as output_handle:
        for key in hits:
            if key.startswith("sgp"):
                bbh_sgp = bbh_grp = bbh_ogp = None
                for li in hits[key]:
                    qseqid, sseqid, pident, length, mismatch, gapopen, qstart, qend, sstart, send, evalue, bitscore = li
                    bitscore = float(bitscore)
                    if sseqid.startswith("sgp"):
                        if bbh_sgp is None or bitscore > bbh_sgp:
                            bbh_sgp = bitscore
                    elif sseqid.startswith("grp"):
                        if bbh_grp is None or bitscore > bbh_grp:
                            bbh_grp = bitscore
                    elif sseqid.startswith("ogp"):
                        if bbh_ogp is None or bitscore > bbh_ogp:
                            bbh_ogp = bitscore
                if bbh_sgp is None:
                    continue
                if bbh_ogp is None:
                    continue
                if bbh_grp is None:
                    bbh_grp = 0
                score = (bbh_ogp - bbh_grp) / bbh_sgp
                if score < 0:
                    continue
                for li in hits[key]:
                    output_handle.write("\t".join(li) + "\n")

def output_besthit_for_subgroup(input_filename: str, output_filename: str) -> None:
    """Slice BLAST hits by crossover group

    Args
    ----
    input_filename : str
        Input filename.
    output_filename : str
        Output filename.
    
    """
    hits = {}

    with open(input_filename, "r") as input_handle:
        for line in input_handle:
            if line.startswith("#"):
                continue
            li = line.strip().split("\t")
            if len(li) != 12:
                continue
            qseqid, sseqid, pident, length, mismatch, gapopen, qstart, qend, sstart, send, evalue, bitscore = li
            if qseqid not in hits:
                hits[qseqid] = []
            hits[qseqid].append(li)

    with open(output_filename, "w") as output_handle:
        for key in hits:
            if key.startswith("sgp"):
                bbh_sgp = bbh_grp = bbh_ogp = None
                for li in hits[key]:
                    qseqid, sseqid, pident, length, mismatch, gapopen, qstart, qend, sstart, send, evalue, bitscore = li
                    bitscore = float(bitscore)
                    if sseqid.startswith("sgp"):
                        if bbh_sgp is None or bitscore > bbh_sgp:
                            bbh_sgp = bitscore
                    elif sseqid.startswith("grp"):
                        if bbh_grp is None or bitscore > bbh_grp:
                            bbh_grp = bitscore
                    elif sseqid.startswith("ogp"):
                        if bbh_ogp is None or bitscore > bbh_ogp:
                            bbh_ogp = bitscore
                if bbh_sgp is None:
                    continue
                if bbh_ogp is None:
                    continue
                if bbh_grp is None:
                    bbh_grp = 0
            score = (bbh_ogp - bbh_grp) / bbh_sgp
            if score < 0:
                continue
            for li in hits[key]:
                output_handle.write("\t".join(li) + "\n")



