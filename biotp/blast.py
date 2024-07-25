def slice_rgo_by_hgt(input_filename, output_filename, pct):
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
            else:
                qry[qseqid]['ogp_bits'].append(bitscore)

    with open(output_filename, 'w') as output_handle:
        for q in qry:
            ## Alien Index
            if not qry[q]['ogp_bits']:
                continue
            bbh_rec = max(qry[q]['rec_bits'])
            bbh_ogp = max(qry[q]['ogp_bits'])

            if not qry[q]['grp_bits']:
                bbh_grp = 0
            else:
                bbh_grp = max(qry[q]['grp_bits'])
                
            ai = (float(bbh_ogp) / float(bbh_rec)) - (float(bbh_grp) / float(bbh_rec))

            ## Percentage of outgroup species
            num = len(qry[q]['bitscore'])
            num_ogp = len(qry[q]['ogp_bits'])
            ogp_pct = int(num_ogp) / int(num) * 100

            if float(ai) > 0 and float(ogp_pct) >= float(pct):
                for sseqid, pident, length, mismatch, gapopen, qstart, qend, sstart, send, evalue, bitscore in zip(qry[q]['sseqid'], qry[q]['pident'], qry[q]['length'], qry[q]['mismatch'], qry[q]['gapopen'], qry[q]['qstart'], qry[q]['qend'], qry[q]['sstart'], qry[q]['send'], qry[q]['evalue'], qry[q]['bitscore']):
                    output_handle.write(f"{q}\t{sseqid}\t{pident}\t{length}\t{mismatch}\t{gapopen}\t{qstart}\t{qend}\t{sstart}\t{send}\t{evalue}\t{bitscore}\n")


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
            else:
                qry[qseqid]['ogp_bits'].append(bitscore)

    with open(output_filename, 'w') as output_handle:
        for q in qry:
            if qry[q]['ogp_bits']:
                for sseqid, pident, length, mismatch, gapopen, qstart, qend, sstart, send, evalue, bitscore in zip(
                        qry[q]['sseqid'], qry[q]['pident'], qry[q]['length'], qry[q]['mismatch'],
                        qry[q]['gapopen'], qry[q]['qstart'], qry[q]['qend'], qry[q]['sstart'],
                        qry[q]['send'], qry[q]['evalue'], qry[q]['bitscore']):
                    output_handle.write(
                        f"{q}\t{sseqid}\t{pident}\t{length}\t{mismatch}\t{gapopen}\t{qstart}\t{qend}\t{sstart}\t{send}\t{evalue}\t{bitscore}\n")

