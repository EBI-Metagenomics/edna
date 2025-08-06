import argparse
import math
import sys
import re
import fileinput
import json
from collections import defaultdict

parser = argparse.ArgumentParser(description='Parse domain table output from HMMer (`hmmsearch`) to calculate HMM model coverages from metagenomic reads.')
parser.add_argument('-i', "--input_fp", type=str,
                    default='-',
                    help="Input fasta/fastq filepath. Use '-' for STDIN (default).")
parser.add_argument('-o', "--output_fp", type=str,
                    default='-',
                    help="Output TSV filepath. Use '-' for STDOUT (default).")
parser.add_argument('-s', "--stats_output_fp", type=str,
                    default='',
                    help="Output JSON filepath for recording mapping stats.")
parser.add_argument('-j', "--json_fp", type=str, 
                    required=False,
                    help="fastp JSON file for total reads.")

args = parser.parse_args()

# --- Read total reads from fastp JSON ---
with open(args.json_fp) as f:
    fastp_data = json.load(f)
total_reads = fastp_data.get("summary", {}).get("after_filtering", {}).get("total_reads", None)

if total_reads is None:
    sys.stderr.write("[ERROR] Could not find after_filtering.total_reads in fastp JSON.\n")
    sys.exit(1)


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
    d, b, e = int(line_dict['read_frame']), int(line_dict['read_frame_begin']), int(line_dict['read_frame_end'])
    if d<0:
        b, e = e, b

    return {
        'read_frame': (d,b,e),
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

        read_header_split = re.findall(r'^(.*?)_frame=(-?\d+)_begin=(\d+)_end=(\d+)\s*$',
                                       line_dict['target_name'])[0]
        line_dict['read_name'] = read_header_split[0]
        line_dict['read_frame'] = read_header_split[1]
        line_dict['read_frame_begin'] = read_header_split[2]
        line_dict['read_frame_end'] = read_header_split[3]

        read_hits[line_dict['read_name']].append(extract_from_line(line_dict))

    top_read_hits = {}
    for k,vs in read_hits.items():
        # greedy resolution of overlaps
        deoverlapped = []
        ali_coverage = set()
        for d in sorted(vs, key=lambda x:x['overall_evalue']):
            phase = int(d['read_frame'][0])
            direction = -1 if phase<0 else 1
            phase *= direction
            start, end = d['read_frame'][1:3]

            m = lambda x: (start-1)+direction*(x-1)*3 + (phase-1)
            nt_base_idxs = list(range(*list(sorted((m(d['ali_coord_from']), m(d['ali_coord_to']))))))

            if not any([i in ali_coverage for i in nt_base_idxs]):
                deoverlapped.append(d)
                for i in nt_base_idxs:
                    ali_coverage.add(i)

        top_read_hits[k] = list(deoverlapped)

    # get hmm coverage and read counts
    hmm_hits_coverage = {}
    hmm_hit_count = defaultdict(set)
    for k,vs in top_read_hits.items():
        for d in vs:
            hmm_hit_count[d['query_accession']].add(k)
            if not d['query_accession'] in hmm_hits_coverage:
                hmm_hits_coverage[d['query_accession']] = {i+1:0 for i in range(d['qlen'])}
            for i in range(d['hmm_coord_from'], d['hmm_coord_to']):
                hmm_hits_coverage[d['query_accession']][i] += 1
    hmm_hit_count = {k:len(vs) for k,vs in hmm_hit_count.items()}

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
            'norm_count' : hmm_hit_count[k] / total_reads if total_reads > 0 else 0
        }

    outfile = sys.stdout if args.output_fp=='-' else open(args.output_fp, 'wt')
    outfile.write("#function\tread_count\tcoverage_depth\tcoverage_breadth\treads_percentage\n")
    for k,d in sorted(hmm_hits_coverage_stats.items(), key=lambda x:-x[1]['depth']):
        outfile.write(f"{k}\t{d['count']}\t{d['depth']}\t{d['breadth']}\t{d['norm_count']}\n")
    outfile.close()

    if args.stats_output_fp:
        stats = {
            'reads_mapped': len(top_read_hits),
            'hmm_count': len(hmm_hit_count),
            'read_hit_count': sum(list(hmm_hit_count.values()))
        }
        with open(args.stats_output_fp, 'wt') as f:
            json.dump(stats, f)

