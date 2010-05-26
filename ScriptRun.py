import subprocess
import sys

bowtie_dir =  '~/bowtie-0.10.1/'
workspace_dir = '~/workspace/EmFusion/'

EM_bin = workspace_dir + 'Release/EmFusion'
bowtie_bin = bowtie_dir + 'bowtie'
bowtie_build_bin = bowtie_dir + 'bowtie_build'

exon_structure_file = ''
exon_seq_file = ''

initial_bowtie_flags = '-l 20 -e 250 -n 3 -y -a -m 150 -p 4 --nomaqround'
pe_bowtie_flags = '-a -p 4 -y -X 800 -m 100 --nomaqround'

offset = 33

def main(fastq1, fastq2, reference):
    
    
    p = subprocess.Popen(' '.join([bowtie_bin, initial_bowtie_flags, '--maxfa 1.repeat',
                                  reference, fastq1, fastq1 + '.bowtie']), shell = True)
    p.wait()
    p = subprocess.Popen(' '.join([bowtie_bin, initial_bowtie_flags, '--maxfa 2.repeat',
                                  reference, fastq2, fastq2 + '.bowtie']), shell = True)
    p.wait()
    
    p = subprocess.Popen(' '.join(['python', workspace_dir + 'PythonSort.py',
                                   'bowtie', fastq1 + '.bowtie']), shell = True)
    p.wait()
    p = subprocess.Popen(' '.join(['python', workspace_dir + 'PythonSort.py',
                                   'bowtie', fastq2 + '.bowtie']), shell = True)
    p.wait()
    
    p = subprocess.Popen(' '.join(['python', workspace_dir + 'PythonSort.py',
                                   'fasta', '1.repeat']), shell = True)
    p.wait()
    p = subprocess.Popen(' '.join(['python', workspace_dir + 'PythonSort.py',
                                   'fasta', '2.repeat']), shell = True)
    p.wait()
    
    p = subprocess.Popen(' '.join(['python', workspace_dir + 'PythonSort.py',
                                   'fastq', fastq1]), shell = True)
    p.wait()
    p = subprocess.Popen(' '.join(['python', workspace_dir + 'PythonSort.py',
                                   'fastq', fastq2]), shell = True)
    p.wait()

    p = subprocess.Popen(' '.join([EM_bin, 'sift', fastq1 + '.bowtie', fastq2 + '.bowtie',
                                   fastq1, fastq2, offset, '1.repeat', '2.repeat']),
                                   shell = True)
    p.wait()
    
    p = subprocess.Popen(' '.join(['sort -k 2,2 -k 3,3', fastq1 + '.discord', '>',
                                   fastq1 + '.discord.sorted']), shell = True)
    p.wait()
    
    p = subprocess.Popen(' '.join(['python', workspace_dir + 'ExtendReference.py',
                                   exon_structure_file, fastq1 + '.discord.sorted']),
                                   shell = True)
    p.wait()
    
    p = subprocess.Popen(' '.join(['sort -k 2, -u', fastq1 + '.discord.sorted.exons', 
                                   '>', fastq1 + '.discord.sorted.exons.uniq']),
                                    shell = True)
    p.wait()
    
    p = subprocess.Popen(' '.join(['python', workspace_dir + 'MakeFusionTranscripts.py',
                                    exon_seq_file, fastq1 + '.discord.sorted.exons.uniq']),
                                    shell = True)
    p.wait()
    
    p = subprocess.Popen(' '.join(['cat', fastq1 + '.discord.sorted.exons.uniq.seqs',
                                   reference + '.fa', '>', 'augmented_ref.fa']),
                                   shell = True)
    p.wait()
    
    p = subprocess.Popen(' '.join([bowtie_build_bin, '-o 0', 'augmented_ref.fa',
                                   'augmented_ref']))
    p.wait()
    
    p = subprocess.Popen(' '.join([bowtie_bin, pe_bowtie_flags, '--unfa pe.unmapped',
                                   'augmented_ref', '-1', fastq1, '-2', fastq2,
                                    fastq1 + '.pe.bowtie']), shell = True)
    p.wait()

    p = subprocess.Popen(' '.join(['cat', 'pe_1.unmapped', 'pe_2.unmapped', '>',
                                   'pe.unmapped']), shell = True)
    p.wait()

    p = subprocess.Popen(' '.join(['sh', workspace_dir + 'sortPEBowtie_Mapping.sh', 
                                  fastq1 + '.pe.bowtie']), shell = True)
    p.wait()
    
    p = subprocess.Popen(' '.join(['sort -n -k 1,1', fastq1 + '.mapdist', '>',
                                   fastq1 + '.mapdist.sorted']), shell = True)
    p.wait()

    p = subprocess.Popen(' '.join(['python', workspace_dir + 'getdistprob.py', 
                                   fastq1 + '.mapdist.sorted', 'dist.prob']),
                                   shell = True)
    p.wait()

    p = subprocess.Popen(' '.join([EM_bin, 'EM', fastq1 + '.pe.bowtie.sorted', 'dist.prob',
                                   'augmented_ref.fa', 'pe.unmapped', '33']), shell = True)
    p.wait()


if __name__ == '__main__':
    main(*sys.argv[1:])
