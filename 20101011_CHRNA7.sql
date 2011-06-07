select names('sus_scrofa_sscrofa9_56_gtf');
select read_table($$ file:'D:/Downloads/Homo_sapiens.GRCh37.59.gtf', dest_table:'homo_sapiens_grch37_59_gtf', sep:'\t', header:["rname", "source", "feature", "f_start", "f_end", "score", "strand", "frame", "attributes"], overwrite:True $$)

alter table homo_sapiens_grch37_59_gtf set tablespace hdd_free_agent

select * from homo_sapiens_grch37_59_gtf where attributes like '%ENSG00000175344%' limit 10;

copy 
  (select * from homo_sapiens_grch37_59_gtf where rname = '15' and f_start between 32322701 - 1000000 and 32322701 + 1000000 order by f_start)
to 'D:/Tritume/hs_chrna7.gtf' with csv header delimiter E'\t';

copy(
select distinct rname, source, min(f_start) AS min_start, strand, gtf_attribute(attributes, 'gene_id') AS gene_id, gtf_attribute(attributes, 'transcript_id') AS transcript_id, gtf_attribute(attributes, 'gene_name') AS gene_name
from sus_scrofa_sscrofa9_56_gtf where rname like '1' and f_start between 144186000 and 151041173 
group by rname, source, strand, gtf_attribute(attributes, 'gene_id'), gtf_attribute(attributes, 'transcript_id'), gtf_attribute(attributes, 'gene_name')
order by min_start
) to 'D:/Tritume/chrna7_pig.gtf' with csv header delimiter E'\t';