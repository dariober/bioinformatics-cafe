#!/usr/bin/python

""" Start aligning a FASTQ file with reads right-trimmed (e.g. CAGE_050810_21bp.fq)
against a reference using BWA.
After this alignment is done, collect the multi-mapping reads and add one base to the
read sequence (i.e. match against the un-trimmed FASTQ, CAGE_50810_no_adapter.fq, and make it 22 bp). 
This is done with rescue_multimappers.py.
Re-start the loop until the last multimappers are full length (27 bp).

This strategy should maximize the number of mappable reads because we start without the
poor quality 5' end. But it should minimize the multimappers (because multimappers
are discriminated by adding one base at a time)

"""

import os
import re
import sys
import subprocess
import time

# -----------------------------[ Input/Output ]--------------------------------

## Full length FASTQ file, without adapters:
fastq_27bp= '/exports/work/vet_roslin_nextgen/dario/fastq/20100805_CAGE_AM/CAGE_050810_no_adapter.fq'

## Reference sequence for BWA
ref_db='/exports/work/vet_roslin_nextgen/dario/bwa/indexes/Sus_scrofa.Sscrofa9.56.dna.toplevel.fa'

## This is the 'start' fastq file, the shortest. Put the whole path:
fastq = '/exports/work/vet_roslin_nextgen/dario/bwa/output/20101218_cage_050810_increment/CAGE_050810_21bp.fq'

## Where input/output file will be:
cwd= '/exports/work/vet_roslin_nextgen/dario/bwa/output/20101218_cage_050810_increment'

## Start read length (shortest, e.g. 21) and finish (probably the full length, e.g. 27 for CAGE)
start_len= 21
end_len= 27

## Where bwa executable is:
bwa_exe= '/exports/work/vet_roslin_nextgen/dario/bwa/bwa-0.5.8c/bwa'

## Where rescue_multimappers.py is (get dir and exact name right!):
rescue_multimappers= '/exports/work/vet_roslin_nextgen/dario/bwa/output/20101218_cage_050810_increment/20101218_rescue_multimappers.py'

# -----------------------------------------------------------------------------

os.chdir(cwd)
job_id= None
for i in range(start_len, end_len+1):
    """ Create .sh file """

    """ BWA commands """
    out_bwa= 'aln' + str(i) + '.bwa'
    out_bwa= os.path.join(cwd, out_bwa)
    
    out_sam= 'aln' + str(i) + '.sam'
    out_sam= os.path.join(cwd, out_sam)
    
    job= 'bwa_' + str(i) + '.sh'
    script_bwa= open(job, 'w')
    script_bwa.write("""#!/bin/bash\n""")
    script_bwa.write(bwa_exe + ' aln -t 4 ' + ref_db + ' ' + fastq + ' > ' + out_bwa) ## Change the number of bwa threads here with -t
    script_bwa.write('\n')
    script_bwa.write(bwa_exe + ' samse ' + ref_db + ' ' + out_bwa  + ' ' + fastq + ' > ' + out_sam)
    script_bwa.write('\n\n')
    
    """  Prepare input for rescue_multimappers.py 
    IMPORTANT: Input SAM for rescue_multimappers.py has to be called bwa_output.sam. 
               Therefore the output of BWA is going to be copied to bwa_output.sam
    rescue_multimappers.py returns the reads that have more then one best alignment (X0:i:2 or more).
    """
    bwa_output= os.path.join(cwd, 'bwa_output.sam')
    script_bwa.write('cp ' + out_sam + ' ' + bwa_output)
    script_bwa.write('\n')

    """ Run rescue_multimappers.py """
    script_bwa.write('python ' + rescue_multimappers + '\n\n' )

    """ Get the output of rescue_multimappers.py and make it suitable for next job 
    The output of rescue_multimappers.py (a FASTQ) is called bwa_input.fq. Rename it
    before feeding it to the next bwa round.
    """
    fastq= str(i+1) + '_bwa_input.fq'  ## E.g. new fastq name will be 22_bwa_input.fq where 22 is the length of the reads
    fastq= os.path.join(cwd, fastq)
    bwa_input= os.path.join(cwd, 'bwa_input.fq')
    script_bwa.write('cp ' + bwa_input + ' ' + fastq  + '\n\n')
    script_bwa.close()
    
    """ The actual script should look something like this:
    #!/bin/bash
    /exports/work/vet_roslin_nextgen/dario/bwa/bwa-0.5.8c/bwa aln -t 4 /exports/work/vet_roslin_nextgen/dario/bwa/indexes/Sscrofa9.56_plus_Human_ribosomal_DNA_complete_repeating_unit.fa /exports/work/vet_roslin_nextgen/dario/bwa/output/20101213c_cage_iterate/22_bwa_input.fq > /exports/work/vet_roslin_nextgen/dario/bwa/output/20101213c_cage_iterate/aln22.bwa
    /exports/work/vet_roslin_nextgen/dario/bwa/bwa-0.5.8c/bwa samse /exports/work/vet_roslin_nextgen/dario/bwa/indexes/Sscrofa9.56_plus_Human_ribosomal_DNA_complete_repeating_unit.fa /exports/work/vet_roslin_nextgen/dario/bwa/output/20101213c_cage_iterate/aln22.bwa /exports/work/vet_roslin_nextgen/dario/bwa/output/20101213c_cage_iterate/22_bwa_input.fq > /exports/work/vet_roslin_nextgen/dario/bwa/output/20101213c_cage_iterate/aln22.sam

    cp /exports/work/vet_roslin_nextgen/dario/bwa/output/20101213c_cage_iterate/aln22.sam /exports/work/vet_roslin_nextgen/dario/bwa/output/20101213c_cage_iterate/bwa_output.sam
    python /exports/work/vet_roslin_nextgen/dario/bwa/output/20101213c_cage_iterate/rescue_multimappers.py

    cp /exports/work/vet_roslin_nextgen/dario/bwa/output/20101213c_cage_iterate/bwa_input.fq /exports/work/vet_roslin_nextgen/dario/bwa/output/20101213c_cage_iterate/23_bwa_input.fq
    """    


    """ Send job to eddie. 
    Make sure the number of threads is consistent with bwa -t.
    """
    if job_id is None:
        """ This is for the first job """
        os.system('qsub -P vet_roslin -pe memory 4 -l h_rt=02:00:00 ' + job)
    else:
        """ For the following jobs, you need to make the script wait for the end of the previous one 
        by setting -hold_jid jobid. See below how to get jobid """
        os.system('qsub -P vet_roslin -pe memory 4 -l h_rt=02:00:00 -hold_jid ' +  job_id + ' ' + job)
    
    
    """ Get job_id """
    time.sleep(20) ## This is to give qstat the time to update
    p= subprocess.Popen('qstat', shell= True, stdout= subprocess.PIPE)
    out= p.stdout.readlines()
    for j in out:
        if job in j:
            j= j.split(' ')
            job_id= j[0]
            break
    print('Job: ' + job + ' job_id: ' + job_id)
sys.exit()


    