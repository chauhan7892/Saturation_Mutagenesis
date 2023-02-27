import argparse
from datetime import datetime 

def fastqProcess(lines=None):
    ks = ['name', 'sequence', 'optional', 'quality']
    return {k: v for k, v in zip(ks, lines)}

def main( ):
    start_time = datetime.now()
    parameters = args_.config_argument
    input_file = args_.input_argument
    output_file = args_.output_argument

    MIN_QUALITY = int(parameters[0])

    with open(output_file[0],'w') as f_out_fasta, open(output_file[1],'w') as f_out_id:

        n = 4
        with open(input_file[0], 'r') as f_in_fastq:
            lines = []
            for line in f_in_fastq:
                lines.append(line)
                if len(lines) == n:
                    record = fastqProcess(lines)
                    quality = record['quality'].strip('\n')
                    phred_score = 0
                    for q in quality:
                        phred_score += (ord(q) - 33)

                    if (phred_score/len(quality)) >= MIN_QUALITY:
                        f_out_fasta.write('>'+record['name'])
                        f_out_fasta.write(record['sequence'])
                    else:
                        f_out_id.write(record['name'] + '\t' + str(phred_score) +'\n')

                    lines = []


    end_time = datetime.now()
    print(f"PROGRAM IS COMPLETED||Duration||H:M:S||{end_time - start_time}|")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="A script to parse fastq file")
    parser.add_argument('-c', nargs='+', dest="config_argument",default="something", help="input parameters")
    parser.add_argument('-i', nargs='+', dest="input_argument",default="something.fq", help="input a fastq file")
    parser.add_argument('-o', nargs='+', dest="output_argument",default="something.fa", help="write to a fasta file")
    args_ = parser.parse_args()
    main(  )
    print('***** well done *****')
