import argparse
import gzip

parser = argparse.ArgumentParser(description='Add the length of the sequence to the header of the sequence in a fasta or fastq file.')
parser.add_argument('-i', "--input_fp", type=str,
                    required=True,
                    help="Input fasta/fastq filepath.")
parser.add_argument('-o', "--output_fp", type=str,
                    required=True,
                    help="Output fasta/fastq filepath.")
parser.add_argument('-z', "--output_gzip", action='store_true',
                    help="Gzip the output. If not specified then will only gzip when inputs are gzipped.")
parser.add_argument('-c', "--no_output_gzip", action='store_false', dest='output_gzip',
                    help="Don't gzip the output. If not specified then will only gzip when inputs are gzipped.")
parser.add_argument('-p', "--length_prefix", type=str,
                    default='_length=',
                    help="Prefix for the length text to be embedded in the header.")
args = parser.parse_args()

def get_reads(infile, gz: bool):
    prev_line = None
    lines = []
    header = None
    after_header = None
    counting = False
    count = 0
    for line in infile:
        line_str = line.decode('utf8') if gz else line
        line_str_strip = line_str.strip()
        if line_str[0]=='+':
            counting = False
        if (not prev_line=='+') and (line_str[0] in {'>','@'}):
            if not header is None:
                l = sum([len(v) for v in lines[1:]])
                new_header = header + args.length_prefix + str(count)
                if after_header:
                    new_header += f" {after_header}"
                yield (new_header, '\n'.join(lines[1:]))
            lines = []
            header = line_str.split()[0]
            after_header = ' '.join(line_str.split()[1:])
            counting = True
            count = 0
        else:
            if counting:
                count += len(line_str_strip)
        lines.append(line_str_strip)
        prev_line = line_str_strip
    else:
        new_header = header + args.length_prefix + str(count)
        if after_header:
            new_header += f" {after_header}"
        yield (new_header, '\n'.join(lines[1:]))

if __name__ == '__main__':
    gz = args.input_fp[-3:]=='.gz'
    in_file = gzip.open(args.input_fp, 'rb') if gz else open(args.input_fp, 'rt')

    out_gz = (args.output_gzip) if (args.output_gzip is not None) else (args.output_fp[-3:]=='.gz')
    out_file = gzip.open(args.output_fp, 'wb') if out_gz else open(args.output_fp, 'wt')

    for k, v in get_reads(in_file, gz):
        out_str = f"{k}\n{v}\n\n"
        out_file.write(out_str.encode() if out_gz else out_str)

    in_file.close()
    out_file.close()

