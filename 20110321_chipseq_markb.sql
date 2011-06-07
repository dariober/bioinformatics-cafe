/*
  Mark Barnett chipseq experiment (H3K4me3 on ctrl, lps, rabbit IgG).
  - Raw sequences aligned to mm9 with BWA, default settings. 
  - BED files generated from reads with MAPQ >= 15 (samtools view -q 15).
  - Peaks detecetd by findPeaks using IgG as control, -size 1000 -minDist 1000.
  - Peaks common to LPS and CTRL identified with intersectBed. 70% reciprocal overlap required (-r -f 0.7)
  
*/

set search_path to 'chipseq';

drop table if exists chipseq.findpeaks;
CREATE TABLE chipseq.findpeaks
(
  dataset_id text,
  peak_id text,
  chr text,
  "start" integer,
  "end" integer,
  strand text,
  score double precision,
  region_size double precision,
  total_tags double precision,
  control_tags double precision,
  control_fc double precision,
  clonal_fc double precision
);
comment on table chipseq.findpeaks is 'Output of HOMER findPeaks. See 20110321_chipseq_markb.sql. LabBook 21/03/2011.';

/*
  findPeaks output
*/

select read_table($$ table: 'chipseq.findpeaks', file: 'F:/data/20110315_markb_chipseq/con2_k4-1.peaks', apply:
"""
def apply(line):
    line= 'con2' + '\t' + line
    return(line)     
""", append: True, skip: 36 $$)

select read_table($$ table: 'chipseq.findpeaks', file: 'F:/data/20110315_markb_chipseq/lps2_k4-1.peaks', apply:
"""
def apply(line):
    line= 'lps2' + '\t' + line
    return(line)     
""", append: True, skip: 36 $$)

/*
  Upload peak coverage as output from coverageBed
*/

select public.read_table($$ file:'D:/Tritume/con2_k4.coverageBed.gtf', table: 'chipseq.peakcoverage', 
  header: ['rname', 'source', 'feature', 'f_start', 'f_end', 'score', 'strand', 'frame', 'attributes', 'read_count', 'no_zero_coverage', 'feature_length', 'no_zero_fract'],
  overwrite: True,
  apply: 
"""def apply(line):
    line= line.split('\t')
    line[1]= 'con2'
    line= '\t'.join(line)
    return(line) """$$);

select public.read_table($$ file:'D:/Tritume/lps2_k4.coverageBed.gtf', table: 'chipseq.peakcoverage', 
  header: ['rname', 'source', 'feature', 'f_start', 'f_end', 'score', 'strand', 'frame', 'attributes', 'read_count', 'no_zero_coverage', 'feature_length', 'no_zero_fract', 'peak_id'],
  append: True,
  apply: 
"""def apply(line):
    line= line.split('\t')
    line[1]= 'lps2'
    line= '\t'.join(line)
    return(line) """$$);
select * from chipseq.peakcoverage order by read_count desc limit 10;
comment on table chipseq.peakcoverage is 'Output of coverageBed to count reads in peaks defined in a GTF file. See labbook 22/03/2011';
comment on column chipseq.peakcoverage.attributes is 'The attribute peak_id comes from mergePeaks. This is NOT the peak id used for edgeR/DEseq';

alter table chipseq.peakcoverage add column peak_id text;
update chipseq.peakcoverage set peak_id = rname || '_' || f_start || '_' || f_end || '_' || strand;
comment on column chipseq.peakcoverage.peak_id is 'Identifier of the peak used for edger/deseq etc...'

select count(distinct peak_id) from chipseq.peakcoverage limit 10;


--
-- edgeR output
--
-- tmp_edger_chipseq comes from 20110323_edger_chipseq_markb.R.
alter table tmp_edger_chipseq set schema chipseq;
alter table chipseq.tmp_edger_chipseq rename to edger_toptable;
comment on table chipseq.edger_toptable is 'Output of edgeR topTable function. Data is Mark Barnett chip-seq H3K4me3. See labbook 21/03/2011 and 20110321_chipseq_markb.sql'
drop table tmp_edger_chipseq;

select * from chipseq.edger_toptable where fdr < 0.01 order by abs(logfc) desc;
select * from chipseq.edger_toptable where peak_id in( 'chr1_3395081_3396335_+', 'chr1_4774332_4777028_+');


/*
  Prepare table for Homer's annotatePeaks.pl
*/
copy (
select distinct 
    chipseq.deseq_nbinomtest.id,
    chipseq.peakcoverage.rname,
    chipseq.peakcoverage.f_start,
    chipseq.peakcoverage.f_end,
    chipseq.peakcoverage.strand,
    chipseq.deseq_nbinomtest.*
from chipseq.peakcoverage inner join chipseq.deseq_nbinomtest on peak_id = id) 
to 'D:/Tritume/deseq_peaks.txt';
-- Send this table to annotatePeaks.pl
-- Get the output back
select read_table($$ file:'D:/Tritume/annotation-1.txt.gz', table:'chipseq.annotate_peaks', select: [0,1,2,3,4, range(7,18)], skip: 1, overwrite:True,
    header:['peak_id', 'rname', 'f_start', 'f_end', 'strand', 'skip1', 'skip2', 'Annotation', 'Detailed Annotation', 'Distance to TSS', 'Nearest PromoterID', 'PromoterID', 'Nearest Unigene', 'Nearest Refseq', 'Nearest Ensembl', 'Gene Name', 'Gene Alias', 'Gene Description']
    $$);
select * from chipseq.annotate_peaks limit 10;

-- Put annotation and precipitation level togheter:
copy (
select * from chipseq.annotate_peaks inner join chipseq.deseq_nbinomtest on 
chipseq.annotate_peaks.peak_id = chipseq.deseq_nbinomtest.id order by padj
) to 'F:/data/20110315_markb_chipseq/annotation_deseq.txt' with csv header delimiter E'\t';

----------------------------------[ TRITUME ]----------------------------------
"intron (NM_012008, intron 1 of 16)"

/*
  Annotation of peaks in peakcoverage.
  Better than annotating findpeaks since peakcoverage has peaks common to both.
*/

create temp table commonpeaks as (
    select distinct peak_id, rname, f_start, f_end, strand from chipseq.peakcoverage
);
create unique index indx_peakposition on commonpeaks(f_start, f_end, rname, strand);
vacuum analyze commonpeaks;

select commonpeaks.*, ensembl.mus_musculus_ncbim37_61_gtf.*
from ensembl.mus_musculus_ncbim37_61_gtf inner join commonpeaks on
    'chr' || ensembl.mus_musculus_ncbim37_61_gtf.rname = commonpeaks.rname and
    commonpeaks.f_start between ensembl.mus_musculus_ncbim37_61_gtf.f_start - 1000 and ensembl.mus_musculus_ncbim37_61_gtf.f_start + 1000
limit 100;

select chipseq.deseq_nbinomtest.*, ensembl.mus_musculus_ncbim37_61_gtf.*
from ensembl.mus_musculus_ncbim37_61_gtf inner join chipseq.deseq_nbinomtest on
    'chr' || ensembl.mus_musculus_ncbim37_61_gtf.rname = commonpeaks.rname and
    commonpeaks.f_start between ensembl.mus_musculus_ncbim37_61_gtf.f_start - 1000 and ensembl.mus_musculus_ncbim37_61_gtf.f_start + 1000
limit 100;

select * from chipseq.peakcoverage where rname like 'chr5' and f_start between 30339701 - 10000 and 30339701 + 10000;

select * from chipseq.edger_toptable inner join chipseq.deseq_nbinomtest on peak_id = id 
where id like 'chr3_96361177_96368399_+'

/*
  Merged peaks (intersectBed output).
*/	
select public.read_table($$ table: 'intersectbed', file: 'D:/Tritume/merge-1.bed', overwrite:True, temp: True,
    header: [ 'rname_a', 'start_a', 'end_a', 'strand_a', 'peak_id_a', 'file_a', 'rname_b', 'start_b', 'end_b', 'strand_b', 'peak_id_b', 'file_b'] $$);

-- Peaks that are most far apart.
select *, 
    (start_a + end_a)/2 as con2_mid, 
    (start_b + end_b)/2 as lps2_mid,
    abs(((start_a + end_a)/2) - ((start_b + end_b)/2)) AS peak_dist,
    ((end_a - start_b)::float / (end_a - start_a)::float) AS perc_overlap
from intersectbed 
order by perc_overlap asc
limit 100;

select * from intersectbed where start_a >= start_b;

chr2	120994442	120998920

select * from findpeaks where chr = 'chr2' and "start" >= 120994442 and "start" <= 120998920 order by start;
and "start" <= 120998920;

select count(*), dataset_id from findpeaks group by dataset_id;

select * from findpeaks where (peak_id = 'chr11-20135' and dataset_id = 'lps2') or (peak_id = 'chr11-15906' and dataset_id = 'con2') 


select con.dataset_id, 
       con.peak_id, con.chr, 
       con.start, 
       con.end, 
       con.strand, 
       con.score, 
       con.region_size, 
       con.total_tags,
       lps.dataset_id, 
       lps.peak_id,
       lps.start, 
       lps.end, 
       lps.score, lps.region_size, lps.total_tags,       
       (con.start + con.end)/2 AS con_mid, 
       (lps.start + lps.end)/2 AS lps_mid
from
   (select * from chipseq.findpeaks where dataset_id like 'con2') AS con 
inner join
   (select * from chipseq.findpeaks where dataset_id like 'lps2') AS lps
   on con.chr = lps.chr and con.strand = lps.strand
where abs( ((con.start + con.end)/2) - ((lps.start + lps.end)/2) ) < 500
limit 200;

select names($$ t:'findpeaks', quote: '', apply: ['"con." + x'] $$);

select * from chipseq.findpeaks limit 100;

select * from chipseq.merged order by v2, v3 limit 20;

select * from chipseq.findpeaks where peak_id like 'chr1-7945' order by "chr", "start", "end" limit 10000;

copy (select chr, "start", "end", strand, peak_id, dataset_id from chipseq.findpeaks where dataset_id like 'con2' order by chr, start) to 'D:/Tritume/con2_k4-1.bed' with csv delimiter E'\t';
copy (select chr, "start", "end", strand, peak_id, dataset_id from chipseq.findpeaks where dataset_id like 'lps2' order by chr, start) to 'D:/Tritume/lps2_k4-1.bed' with csv delimiter E'\t';

select read_table($$ file:'D:/Tritume/merge-1.bed', table: 'intersect_bed', temp:True $$)
select * from intersect_bed limit 10;
select v5 from intersect_bed;