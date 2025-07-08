import argparse
import math
import sys
import re
import fileinput
from collections import defaultdict

parser = argparse.ArgumentParser(description='Parse domain table output from HMMer (`hmmsearch`) to calculate HMM model coverages from metagenomic reads.')
parser.add_argument('-i', "--input_fp", type=str,
                    default='-', required=False,
                    help="Input fasta/fastq filepath. Use '-' for STDIN (default).")
parser.add_argument('-o', "--output_fp", type=str,
                    default='-', required=False,
                    help="Output TSV filepath. Use '-' for STDOUT (default).")

args = parser.parse_args()


cols = [
    'target_name',
    'target_accession',
    'tlen',
    'query_name',
    'query_accession',
    'qlen',
    'overall_evalue',
    'overall_score',
    'overall_bias',
    'domain_n',
    'domain_total',
    'domain_c_evalue',
    'domain_i_evalue',
    'domain_score',
    'domain_bias',
    'hmm_coord_from',
    'hmm_coord_to',
    'ali_coord_from',
    'ali_coord_to',
    'env_coord_from',
    'env_coord_to',
    'acc',
    'description_of_target',
]

def extract_from_line(line_dict):
    l, d, b, e = int(line_dict['read_length']), int(line_dict['read_frame']), int(line_dict['read_frame_begin']), int(line_dict['read_frame_end'])
    if d<0:
        b = l-b
        e = l-e

    return {
        'read_frame': (l,d,b,e),
        'query_accession': line_dict['query_accession'],
        'overall_score': float(line_dict['overall_score']),
        'overall_evalue': float(line_dict['overall_evalue']),
        'acc': float(line_dict['acc']),
        'tlen': int(line_dict['tlen']),
        'qlen': int(line_dict['qlen']),
        'ali_coord_from': int(line_dict['ali_coord_from']),
        'ali_coord_to': int(line_dict['ali_coord_to']),
        'hmm_coord_from': int(line_dict['hmm_coord_from']),
        'hmm_coord_to': int(line_dict['hmm_coord_to'])
    }

if __name__ == '__main__':
    # Parse file
    read_hits = defaultdict(list)
    for line in fileinput.input([] if args.input_fp=='-' else args.input_fp):
        if line[0] == '#':
            continue
        line_dict = dict(zip(cols, [v.strip() for v in line.strip().split()]))

        match = re.match(r'^([^_]*)_length=(\d+)(?:_frame=(-?\d+)_begin=(\d+)_end=(\d+))?$', line_dict['target_name'])
        if not match:
            sys.stderr.write(f"[WARN] Skipping unrecognized target name: {line_dict['target_name']}\n")
            continue

        line_dict['read_name'] = match.group(1)
        line_dict['read_length'] = match.group(2)
        line_dict['read_frame'] = match.group(3) if match.group(3) is not None else '0'
        line_dict['read_frame_begin'] = match.group(4) if match.group(4) is not None else '0'
        line_dict['read_frame_end'] = match.group(5) if match.group(5) is not None else line_dict['read_length']


        read_hits[line_dict['read_name']].append(extract_from_line(line_dict))

    top_read_hits = {}
    for k,vs in read_hits.items():
        # greedy resolution of overlaps
        deoverlapped = []
        ali_coverage = {i:False for i in range(vs[0]['read_frame'][0])}
        for d in sorted(vs, key=lambda x:x['overall_evalue']):
            phase = d['read_frame'][1]
            direction = -1 if phase<0 else 1
            start, end = d['read_frame'][2:4]

            f_sort = lambda a,b: (b,a) if a>b else (a,b)
            m = lambda x: (start-1)+direction*(x-1)*3 + ((phase*direction)-1)
            nt_base_idxs = list(range(*f_sort(m(d['ali_coord_from']), m(d['ali_coord_to']))))
            nt_base_idxs = [i for i in nt_base_idxs if 0 <= i < len(ali_coverage)]

            if not any([ali_coverage[i] for i in nt_base_idxs]):
                deoverlapped.append(d)
                for i in nt_base_idxs:
                    ali_coverage[i] = True


        top_read_hits[k] = list(deoverlapped)

    # get hmm coverage and read counts
    hmm_hits_coverage = {}
    hmm_hit_count = defaultdict(int)
    for k,vs in top_read_hits.items():
        for d in vs:
            hmm_hit_count[d['query_accession']] += 1
            if not d['query_accession'] in hmm_hits_coverage:
                hmm_hits_coverage[d['query_accession']] = {i+1:0 for i in range(d['qlen'])}
            for i in range(d['hmm_coord_from'], d['hmm_coord_to']):
                hmm_hits_coverage[d['query_accession']][i] += 1

    # Collect and write
    hmm_hits_coverage_stats = {}
    for k,d in hmm_hits_coverage.items():
        if not len(d)>0:
            continue
        depth = sum(list(d.values()))/len(d)
        breadth = sum([v>0 for _,v in d.items()])/len(d)
        
        depth_ = 709 if depth>709 else depth  # prevents and float overflow with math.exp
        expected = 1-(1/math.log2(1+math.exp(depth_)))

        hmm_hits_coverage_stats[k] = {
            'depth': depth,
            'breadth': breadth,
            'count': hmm_hit_count[k],
            'expected_breadth': expected,
            'ratio': breadth/expected,
        }

    outfile = sys.stdout if args.output_fp=='-' else open(args.output_fp, 'wt')
    outfile.write(f"# Query Accession\tRead Count\tCoverage Depth\tCoverage Breadth\n")
    for k,d in sorted(hmm_hits_coverage_stats.items(), key=lambda x:-x[1]['depth']):
        outfile.write(f"{k}\t{d['count']}\t{d['depth']}\t{d['breadth']}\n")
    outfile.close()

