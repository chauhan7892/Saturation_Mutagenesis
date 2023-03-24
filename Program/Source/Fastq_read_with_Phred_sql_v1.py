import argparse
from datetime import datetime 
import sqlite3

def fastqProcess(lines=None):
    ks = ['seq_id', 'seq', 'optional', 'quality']
    return {k: v for k, v in zip(ks, lines)}


def main( ):
    start_time = datetime.now()
    parameters = args_.config_argument
    input_file = args_.input_argument
    output_file = args_.output_argument

    MIN_QUALITY = int(parameters[0])

    ## connect to sqlite3 database
    db_path = input_file[0] # sql database 'fasta_database.sqlite'
    table_name = input_file[1] ## table name "fasta_table"
    conn = sqlite3.connect(db_path) ## connect to database
    cursor = conn.cursor()

    # create table
    # cursor.execute('DROP TABLE IF EXISTS fasta_table ') ## drop table
    # cursor.execute('CREATE TABLE IF NOT EXISTS fasta_table (sample TEXT, seq_id TEXT, seq TEXT, PRIMARY KEY (sample, seq_id))') ## hard coded
    cursor.execute(f'CREATE TABLE IF NOT EXISTS {table_name} (sample TEXT, seq_id TEXT, seq TEXT, PRIMARY KEY (sample, seq_id))') ## variable table

    with open(output_file[0],'w') as f_out_id:
        SAMPLE = input_file[2]
        n = 4 ## four lines of fastq at a time 
        with open(input_file[3], 'r') as f_in_fastq:
            lines = [] ## initilize fastq data
            for line in f_in_fastq:
                lines.append(line)
                if len(lines) == n: ## if four lines covered 
                    record = fastqProcess(lines) ## map fastq read 
                    quality = record['quality'].strip('\n')
                    phred_score = 0
                    for q in quality:
                        phred_score += (ord(q) - 33)
                    if (phred_score/len(quality)) >= MIN_QUALITY:
                        # write data to sqlite3 table 
                        record_seq_id = '>'+record['seq_id'].strip('\n')
                        record_seq = record['seq'].strip('\n')
                        cursor.execute('INSERT OR IGNORE INTO fasta_table (sample, seq_id, seq) VALUES (?, ?, ?)', (SAMPLE, record_seq_id, record_seq))
                    else:
                        #print(phred_score)
                        f_out_id.write(record['seq_id'].strip('\n')  + '\t' + str(phred_score/len(quality)) +'\n')

                    lines = [] ## empty the fastq data for next read
                    
    # commit changes and close connection
    conn.commit()
    cursor.close()

    # cursor = conn.cursor()    
    # cursor.execute("SELECT seq_id FROM fasta_table WHERE sample = ?", (SAMPLE,))
    # seq_ids_tuple = cursor.fetchall()
    # seq_ids_list = [x for x, in seq_ids_tuple]
    # for seq_id in seq_ids_list:
    #     cursor.execute("SELECT seq FROM fasta_table WHERE (seq_id = ? AND sample = ?)", (str(seq_id), SAMPLE,))
    #     fasta_seq = cursor.fetchone()[0]
    #     print(fasta_seq)

    # cursor.close()
    conn.close()

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
